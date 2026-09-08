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
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"sync"
)

// singleAlignment holds the details of one alignment of a discrete srna.
type singleAlignment struct {
	Seq          string // Seq is an aligned read sequence
	timesAligned int    // No of times the read has aligned
	// Pos is the 1-based forward-strand coordinate of the leftmost aligned base for
	// both strands. For "+" alignments this is the read's 5' end; for "-" alignments
	// it is the read's 3' end.
	Pos        int
	Strand     string    // Strand
	Alignments []float64 // Per-read-file counts (may be shared between alignments of the same read)
}

// collection for sorting
type singleAlignments []*singleAlignment

// sortAlignments orders alignments by position, then strand ("+" before "-"), then
// read sequence, so output row order is deterministic across runs.
func sortAlignments(alignments singleAlignments) {
	sort.Slice(alignments, func(i, j int) bool {
		a, b := alignments[i], alignments[j]
		if a.Pos != b.Pos {
			return a.Pos < b.Pos
		}
		if a.Strand != b.Strand {
			return a.Strand < b.Strand
		}
		return a.Seq < b.Seq
	})
}

// ProfileNoSplit takes an alignment map and a sequence map as an input.  It returns a map of single alignments
// with a reference header as key and a single alignments struct as value.  The count for each read alignment
// is NOT split by the number of times a read aligns.
func ProfileNoSplit(alignmentMap map[string]map[string][]int, seqMap map[string][]float64) map[string]*singleAlignments {
	srnaAlignmentMap := calcTimesReadAligns(alignmentMap)
	profileAlignmentsMap := make(map[string]*singleAlignments)
	for header, alignments := range alignmentMap {
		var combinedAlignments singleAlignments
		for srna, positions := range alignments {
			counts := seqMap[srna]
			timesAligned := srnaAlignmentMap[srna]
			for _, position := range positions {
				switch {
				case position > 0:
					combinedAlignments = append(combinedAlignments,
						&singleAlignment{srna, timesAligned, position, "+", counts})
				case position < 0:
					combinedAlignments = append(combinedAlignments,
						&singleAlignment{srna, timesAligned, 0 - position, "-", counts})
				}
			}
		}
		sortAlignments(combinedAlignments)
		profileAlignmentsMap[header] = &combinedAlignments
	}
	return profileAlignmentsMap
}

type alignmentStruct struct {
	header     string
	alignments map[string][]int
}

type outputStruct struct {
	header             string
	combinedAlignments *singleAlignments
}

// calcTimesReadAligns calculates the number of times an aligned read aligns across all references
func calcTimesReadAligns(alignmentMap map[string]map[string][]int) map[string]int {
	srnaAlignmentMap := make(map[string]int)
	for _, alignment := range alignmentMap {
		for srna, pos := range alignment {
			srnaAlignmentMap[srna] += len(pos)
		}
	}
	return srnaAlignmentMap
}

// ProfileSplit takes an alignment map and a sequence map as an input.  It returns a map of single alignments
// with a reference header as key and a single alignments struct as value.  The count for each read alignment
// is split by the number of times the read aligns.
func ProfileSplit(alignmentMap map[string]map[string][]int, seqMap map[string][]float64) map[string]*singleAlignments {
	alignmentsNo := len(alignmentMap)
	wg := &sync.WaitGroup{}

	srnaAlignmentMap := calcTimesReadAligns(alignmentMap)

	alignmentsForGoroutine := make(chan alignmentStruct, alignmentsNo)
	outputFromGoroutine := make(chan outputStruct, alignmentsNo)

	for header, alignments := range alignmentMap {
		alignmentsForGoroutine <- alignmentStruct{header, alignments}
	}
	close(alignmentsForGoroutine)

	workers := runtime.NumCPU()
	if workers > alignmentsNo {
		workers = alignmentsNo
	}
	if workers < 1 {
		workers = 1
	}
	wg.Add(workers)
	for a := 0; a < workers; a++ {
		go profileSplitWorker(alignmentsForGoroutine, outputFromGoroutine, seqMap, srnaAlignmentMap, wg)
	}

	go func() {
		wg.Wait()
		close(outputFromGoroutine)
	}()

	profileAlignmentsMap := make(map[string]*singleAlignments)
	for singleOutputStruct := range outputFromGoroutine {
		profileAlignmentsMap[singleOutputStruct.header] = singleOutputStruct.combinedAlignments
	}

	return profileAlignmentsMap
}

func profileSplitWorker(alignmentsForGoroutine chan alignmentStruct, outputFromGoroutine chan outputStruct,
	seqMap map[string][]float64, srnaAlignmentMap map[string]int, wg *sync.WaitGroup) {
	defer wg.Done()
	for singleAlign := range alignmentsForGoroutine {
		var combinedAlignments singleAlignments
		for srna, positions := range singleAlign.alignments {
			counts := seqMap[srna]
			timesAligned := srnaAlignmentMap[srna]
			splitCounts := make([]float64, len(counts))
			for i, count := range counts {
				splitCounts[i] = count / float64(timesAligned)
			}
			for _, position := range positions {
				switch {
				case position > 0:
					combinedAlignments = append(combinedAlignments,
						&singleAlignment{srna, timesAligned, position, "+", splitCounts})
				case position < 0:
					combinedAlignments = append(combinedAlignments,
						&singleAlignment{srna, timesAligned, 0 - position, "-", splitCounts})
				}
			}
		}
		sortAlignments(combinedAlignments)
		outputFromGoroutine <- outputStruct{singleAlign.header, &combinedAlignments}
	}
}

// ProfileToCsv streams the profile results to a csv file (one row per alignment).
func ProfileToCsv(profileAlignmentsMap map[string]*singleAlignments, refSlice []*HeaderRef, nt int,
	outPrefix string, fileOrder []string) {
	outFile := outPrefix + "_" + strconv.Itoa(nt) + ".csv"
	outDir := filepath.Dir(outFile)
	if err := os.MkdirAll(outDir, 0755); err != nil {
		fmt.Println("Can't create output directory " + outDir + ": " + err.Error())
		errorShutdown()
	}
	f, err := os.Create(outFile)
	if err != nil {
		fmt.Println("Can't open csv file for writing: " + err.Error())
		errorShutdown()
	}
	bw := bufio.NewWriterSize(f, 1<<20)
	w := csv.NewWriter(bw)

	headerRow := append([]string{"Header", "len", "sRNA", "Position", "Strand", "Times aligned"}, fileOrder...)
	if err := w.Write(headerRow); err != nil {
		log.Fatalln("error writing csv:", err)
	}
	row := make([]string, 0, len(headerRow))
	for _, ref := range refSlice {
		alignments, ok := profileAlignmentsMap[ref.Header]
		if !ok {
			continue
		}
		refLen := strconv.Itoa(len(ref.Seq))
		for _, alignment := range *alignments {
			row = row[:0]
			row = append(row, ref.Header, refLen, alignment.Seq, strconv.Itoa(alignment.Pos),
				alignment.Strand, strconv.Itoa(alignment.timesAligned))
			for _, v := range alignment.Alignments {
				row = append(row, strconv.FormatFloat(v, 'f', 6, 64))
			}
			if err := w.Write(row); err != nil {
				log.Fatalln("error writing csv:", err)
			}
		}
	}
	w.Flush()
	if err := w.Error(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
	if err := bw.Flush(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
	if err := f.Close(); err != nil {
		log.Fatalln("error closing csv:", err)
	}
}
