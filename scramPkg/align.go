// BSD 3-Clause License

// Copyright (c) 2022, Stephen Fletcher
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

package scramPkg

import (
	"fmt"
	"runtime"
	"sort"
	"sync"
)

// MaxAlignLen is the longest read length that can be aligned (a 2-bit packed
// read must fit in a uint64).
const MaxAlignLen = 32

// baseCode maps A/C/G/T to 2-bit codes; every other byte is invalid (0xFF).
// Complement of a valid code c is 3-c.
var baseCode = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 0xFF
	}
	t['A'] = 0
	t['C'] = 1
	t['G'] = 2
	t['T'] = 3
	return t
}()

// ReadIndex holds, for each requested read length, a map from the 2-bit packed
// read sequence to the read string. Packed uint64 keys make the per-position
// lookups while scanning a reference much cheaper than hashing 18-32 byte
// string keys, and the per-length buckets keep each map small.
type ReadIndex map[int]map[uint64]string

// IndexReads packs all reads of the requested lengths into a ReadIndex.
// Reads containing non-ACGT characters cannot be packed and are excluded from
// alignment — such reads could previously only exact-match reference runs of
// the same ambiguity characters, which produced phantom alignments.
func IndexReads(seqMap map[string][]float64, lengths []int) ReadIndex {
	index := make(ReadIndex, len(lengths))
	for _, l := range lengths {
		index[l] = make(map[uint64]string)
	}
	skipped := 0
	for srna := range seqMap {
		bucket, ok := index[len(srna)]
		if !ok {
			continue
		}
		packed, ok := packSeq(srna)
		if !ok {
			skipped++
			continue
		}
		bucket[packed] = srna
	}
	if skipped > 0 {
		fmt.Printf("%v reads containing non-ACGT characters excluded from alignment\n", skipped)
	}
	return index
}

// packSeq 2-bit encodes a sequence of up to MaxAlignLen bases.
func packSeq(seq string) (uint64, bool) {
	if len(seq) > MaxAlignLen {
		return 0, false
	}
	var packed uint64
	for i := 0; i < len(seq); i++ {
		code := baseCode[seq[i]]
		if code == 0xFF {
			return 0, false
		}
		packed = packed<<2 | uint64(code)
	}
	return packed, true
}

// refAlignments carries one reference's alignments for every requested length.
type refAlignments struct {
	header   string
	byLength map[int]map[string][]int
}

// AlignReadsMulti aligns reads of all requested lengths to one or more reference
// sequences in a single pass per reference, with exact matches in forward or
// reverse complement accepted. It returns one ref_header:[srna_seq:[pos,pos,...],...]
// map per requested length. Positions are 1-based coordinates of the leftmost
// aligned base on the forward strand, negated for reverse-strand alignments;
// each position slice is sorted ascending.
func AlignReadsMulti(index ReadIndex, refSlice []*HeaderRef, lengths []int) map[int]map[string]map[string][]int {
	uniqueLengths := dedupeLengths(lengths)
	for _, l := range uniqueLengths {
		if l < 1 || l > MaxAlignLen {
			fmt.Printf("Can't align %v nt reads: lengths must be between 1 and %v\n", l, MaxAlignLen)
			errorShutdown()
		}
	}

	wg := &sync.WaitGroup{}
	refSeqChan := make(chan *HeaderRef, len(refSlice))
	for _, headerRefPair := range refSlice {
		refSeqChan <- headerRefPair
	}
	close(refSeqChan)

	resultChan := make(chan refAlignments, len(refSlice))
	workers := runtime.NumCPU()
	if workers > len(refSlice) {
		workers = len(refSlice)
	}
	if workers < 1 {
		workers = 1
	}
	wg.Add(workers)
	for a := 0; a < workers; a++ {
		go func() {
			defer wg.Done()
			for refSeq := range refSeqChan {
				resultChan <- alignSingleRef(index, refSeq, uniqueLengths)
			}
		}()
	}
	go func() {
		wg.Wait()
		close(resultChan)
	}()

	final := make(map[int]map[string]map[string][]int, len(uniqueLengths))
	for _, l := range uniqueLengths {
		final[l] = make(map[string]map[string][]int)
	}
	for res := range resultChan {
		for l, alignments := range res.byLength {
			if len(alignments) > 0 {
				final[l][res.header] = alignments
			}
		}
	}
	return final
}

// AlignReads aligns reads of a single length nt, kept for API compatibility.
// For multiple lengths, IndexReads + AlignReadsMulti scans each reference once in total.
func AlignReads(seqMap map[string][]float64, refSlice []*HeaderRef, nt int) map[string]map[string][]int {
	return AlignReadsMulti(IndexReads(seqMap, []int{nt}), refSlice, []int{nt})[nt]
}

// alignSingleRef rolls a 2-bit packed window (and its reverse complement) along
// one reference sequence, checking every requested length at each position.
func alignSingleRef(index ReadIndex, refSeq *HeaderRef, lengths []int) refAlignments {
	byLength := make(map[int]map[string][]int, len(lengths))
	buckets := make([]map[uint64]string, len(lengths))
	masks := make([]uint64, len(lengths))
	rcShifts := make([]uint, len(lengths))
	for i, l := range lengths {
		byLength[l] = make(map[string][]int)
		buckets[i] = index[l]
		if l == MaxAlignLen {
			masks[i] = ^uint64(0)
		} else {
			masks[i] = uint64(1)<<(2*uint(l)) - 1
		}
		rcShifts[i] = 64 - 2*uint(l)
	}

	seq := refSeq.Seq
	// packed holds the last 32 bases 2-bit encoded, newest base in the low bits;
	// rc holds their complements in reverse order, newest base in the high bits.
	// A window of length l ending at p is packed&masks[i]; its reverse
	// complement is rc>>rcShifts[i]. Windows overlapping a non-ACGT base
	// (position <= lastInvalid) are skipped.
	var packed, rc uint64
	lastInvalid := -1
	for p := 0; p < len(seq); p++ {
		code := baseCode[seq[p]]
		c := uint64(code)
		if code == 0xFF {
			lastInvalid = p
			c = 0
		}
		packed = packed<<2 | c
		rc = rc>>2 | (3-c)<<62
		for i, l := range lengths {
			start := p - l + 1
			if start < 0 || start <= lastInvalid {
				continue
			}
			bucket := buckets[i]
			if srna, ok := bucket[packed&masks[i]]; ok {
				byLength[l][srna] = append(byLength[l][srna], 1+start)
			}
			if srna, ok := bucket[rc>>rcShifts[i]]; ok {
				byLength[l][srna] = append(byLength[l][srna], -1-start)
			}
		}
	}
	for _, alignments := range byLength {
		for _, positions := range alignments {
			sort.Ints(positions)
		}
	}
	return refAlignments{refSeq.Header, byLength}
}

func dedupeLengths(lengths []int) []int {
	seen := make(map[int]bool, len(lengths))
	var unique []int
	for _, l := range lengths {
		if !seen[l] {
			seen[l] = true
			unique = append(unique, l)
		}
	}
	return unique
}
