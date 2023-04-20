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
	"sync"
)

// AlignReads aligns reads of  length nt to one or more reference sequences, with exact matches in forward or reverse
// complement accepted.
// A map of ref_header:[srna_seq:[pos,pos,...],...] is returned.
func AlignReads(seq_map map[string]interface{}, ref_slice []*HeaderRef, nt int) map[string]map[string][]int {
	wg := &sync.WaitGroup{}
	wg.Add(len(ref_slice))

	ref_seq_chan := make(chan *HeaderRef, len(ref_slice))
	for _, header_ref_pair := range ref_slice {
		ref_seq_chan <- header_ref_pair
	}
	close(ref_seq_chan)

	header_map_chan := make(chan map[string]map[string][]int, 1)
	for a := 0; a < len(ref_slice); a++ {
		go worker_go(seq_map, ref_seq_chan, nt, header_map_chan, wg)
	}
	go func(cs chan map[string]map[string][]int, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(header_map_chan, wg)

	final_alignment_map := compile_alignments(header_map_chan)
	return final_alignment_map
}

// Align reads to individual reference sequences
func worker_go(seq_map map[string]interface{}, ref_seqs chan *HeaderRef, nt int,
	header_map_chan chan map[string]map[string][]int, wg *sync.WaitGroup) {

	ref_seq := <-ref_seqs
	header_mapped := make(map[string][]int)
	mapped := make(map[string]map[string][]int, 1)
	position := 0
	ref_seq_len := len(ref_seq.Seq)
	for position <= ref_seq_len-nt {
		//each alignment position for an srna
		fwd_seq := ref_seq.Seq[position : position+nt]
		rvs_seq := ref_seq.ReverseSeq[position : position+nt]
		if _, ok := seq_map[fwd_seq]; ok {
			header_mapped[fwd_seq] = append(header_mapped[fwd_seq], 1+position)
		}
		if _, ok := seq_map[rvs_seq]; ok {
			header_mapped[rvs_seq] = append(header_mapped[rvs_seq], -1-(ref_seq_len-position-nt))
		}
		position++
	}
	if len(header_mapped) > 0 {
		mapped[ref_seq.Header] = header_mapped
	}
	header_map_chan <- mapped
	wg.Done()
}

// compile_alignments compiles the alignments
func compile_alignments(header_map_chan chan map[string]map[string][]int) map[string]map[string][]int {
	final_alignment_map := make(map[string]map[string][]int)
	for combined_alignment := range header_map_chan {
		for header, alignment := range combined_alignment {
			final_alignment_map[header] = alignment
		}
	}
	return final_alignment_map
}

type mean_se_dup struct {
	mean_se interface{}
	dup     float64
}

type counts_dup struct {
	counts []float64
	dup    float64
}
