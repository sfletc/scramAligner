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
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"sync"
)

// Details of an alignment for a discrete srna - meanSe
type singleAlignment struct {
	Seq          string // Seq is an aligned read sequence
	timesAligned int    // No of times the read has aligned
	Pos          int    // Pos is an aligned read position (from 5' fwd, starting at 1)
	Strand       string // Strand
	Alignments   interface{}
}

// collection for sorting
type singleAlignments []*singleAlignment

func (slice singleAlignments) Len() int {
	return len(slice)
}
func (slice singleAlignments) Less(i, j int) bool {
	return slice[i].Pos < slice[j].Pos
}
func (slice singleAlignments) Swap(i, j int) {
	slice[i], slice[j] = slice[j], slice[i]
}

// ProfileNoSplit takes and alignment map and a sequence map as an input.  It returns a map of single alignments
// with a reference header as key and a single alignments struct as value.  Each single alignments struct is comprised of
// single_alignment structs (read seq, position, count, se).  The count for each read alignment is NOT split by the
// number of times a read aligns.
func ProfileNoSplit(alignmentMap map[string]map[string][]int, seqMap map[string]interface{}) map[string]interface{} {
	srnaAlignmentMap := calcTimesReadAligns(alignmentMap)
	profileAlignmentsMap := make(map[string]interface{})
	for header, alignments := range alignmentMap {
		var combinedAlignments singleAlignments
		for srna, positions := range alignments {
			for _, position := range positions {
				switch {
				case position > 0:
					switch v := seqMap[srna].(type) {
					case *meanSe:
						alignment := singleAlignment{srna, srnaAlignmentMap[srna],
							position, "+", &meanSe{v.Mean, v.Se}}
						combinedAlignments = append(combinedAlignments, &alignment)
					case *[]float64:
						alignment := singleAlignment{srna, srnaAlignmentMap[srna],
							position, "+", v}
						combinedAlignments = append(combinedAlignments, &alignment)
					}
				case position < 0:
					switch v := seqMap[srna].(type) {
					case *meanSe:
						alignment := singleAlignment{srna, srnaAlignmentMap[srna],
							0 - position, "-", &meanSe{v.Mean, v.Se}}
						combinedAlignments = append(combinedAlignments, &alignment)
					case *[]float64:
						alignment := singleAlignment{srna, srnaAlignmentMap[srna],
							0 - position, "-", v}
						combinedAlignments = append(combinedAlignments, &alignment)
					}
				}
			}
		}
		sort.Sort(combinedAlignments)
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
	combinedAlignments interface{}
}

// Calculates the number of times an aligned read aligns
func calcTimesReadAligns(alignmentMap map[string]map[string][]int) map[string]int {
	srnaAlignmentMap := make(map[string]int)
	for _, alignment := range alignmentMap {
		for srna, pos := range alignment {
			if _, ok := srnaAlignmentMap[srna]; ok {
				srnaAlignmentMap[srna] += len(pos)
			} else {
				srnaAlignmentMap[srna] = len(pos)
			}
		}
	}
	return srnaAlignmentMap
}

// ProfileSplit takes and alignment map and a sequence map as an input.  It returns a map of single alignments
// with a reference header as key and a single alignments struct as value.  Each single alignments struct is comprised
// of single_alignment structs (read seq, position, count, se).  The count for each read alignment is split by the
// number of times a read aligns.
func ProfileSplit(alignmentMap map[string]map[string][]int, seqMap map[string]interface{}) map[string]interface{} {
	alignmentsNo := len(alignmentMap)
	wg := &sync.WaitGroup{}
	wg.Add(alignmentsNo)

	srnaAlignmentMap := calcTimesReadAligns(alignmentMap)

	alignmentsForGoroutine := make(chan alignmentStruct, alignmentsNo)
	outputFromGoroutine := make(chan outputStruct, alignmentsNo)

	for header, alignments := range alignmentMap {
		alignmentsForGoroutine <- alignmentStruct{header, alignments}

	}
	for a := 0; a < alignmentsNo; a++ {
		go profileSplitWorker(alignmentsForGoroutine, outputFromGoroutine, seqMap, srnaAlignmentMap, wg)
	}

	go func(cs chan outputStruct, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(outputFromGoroutine, wg)

	profileAlignmentsMap := make(map[string]interface{})
	for singleOutputStruct := range outputFromGoroutine {
		profileAlignmentsMap[singleOutputStruct.header] = singleOutputStruct.combinedAlignments
	}

	return profileAlignmentsMap
}

func profileSplitWorker(alignmentsForGoroutine chan alignmentStruct, outputFromGoroutine chan outputStruct,
	seqMap map[string]interface{}, srnaAlignmentMap map[string]int, wg *sync.WaitGroup) {
	singleAlign := <-alignmentsForGoroutine
	var combinedAlignmentsMeanSe singleAlignments
	for srna, positions := range singleAlign.alignments {
		for _, position := range positions {
			switch v := seqMap[srna].(type) {
			case *meanSe:
				splitCountMean := v.Mean / float64(srnaAlignmentMap[srna])
				splitSe := v.Se / float64(srnaAlignmentMap[srna])
				switch {
				case position > 0:
					alignment := singleAlignment{srna, srnaAlignmentMap[srna],
						position, "+", &meanSe{splitCountMean, splitSe}}
					combinedAlignmentsMeanSe = append(combinedAlignmentsMeanSe, &alignment)
				case position < 0:
					alignment := singleAlignment{srna, srnaAlignmentMap[srna],
						0 - position, "-", &meanSe{splitCountMean, splitSe}}
					combinedAlignmentsMeanSe = append(combinedAlignmentsMeanSe, &alignment)
				}
			case *[]float64:
				var splitCounts []float64
				for _, i := range *v {
					splitCounts = append(splitCounts, i/float64(srnaAlignmentMap[srna]))
				}
				switch {
				case position > 0:
					alignment := singleAlignment{srna, srnaAlignmentMap[srna],
						position, "+", &splitCounts}
					combinedAlignmentsMeanSe = append(combinedAlignmentsMeanSe, &alignment)
				case position < 0:
					alignment := singleAlignment{srna, srnaAlignmentMap[srna],
						0 - position, "-", &splitCounts}
					combinedAlignmentsMeanSe = append(combinedAlignmentsMeanSe, &alignment)
				}
			}
		}
	}

	sort.Sort(combinedAlignmentsMeanSe)
	outputFromGoroutine <- outputStruct{singleAlign.header, &combinedAlignmentsMeanSe}
	wg.Done()
}

// ProfileToCsv writes the  den results to a csv file
func ProfileToCsv(profileAlignmentsMap map[string]interface{}, refSlice []*HeaderRef, nt int, outPrefix string, fileOrder []string) {

	firstRow := true
	var rows [][]string
	for _, ref := range refSlice {
		if alignments, ok := profileAlignmentsMap[ref.Header]; ok {
			for _, alignment := range *alignments.(*singleAlignments) {
				switch v := alignment.Alignments.(type) {
				case *meanSe:
					if firstRow {
						rows = [][]string{
							{"Header", "len", "sRNA", "Position", "Strand", "Count", "Std. Err", "Times aligned"},
						}
						firstRow = false
					}
					row := []string{ref.Header, strconv.Itoa(len(ref.Seq)),
						alignment.Seq, strconv.Itoa(alignment.Pos),
						alignment.Strand,
						strconv.FormatFloat(v.Mean, 'f', 3, 64),
						strconv.FormatFloat(v.Se, 'f', 8, 64),
						strconv.Itoa(alignment.timesAligned)}
					rows = append(rows, row)
				case *[]float64:
					if firstRow {
						row := []string{"Header", "len", "sRNA", "Position", "Strand", "Times aligned"}
						row = append(row, fileOrder...)
						rows = append(rows, row)
						firstRow = false
					}

					row := []string{ref.Header, strconv.Itoa(len(ref.Seq)),
						alignment.Seq, strconv.Itoa(alignment.Pos),
						alignment.Strand, strconv.Itoa(alignment.timesAligned)}
					pos := 0
					for pos < len(*v) {
						row = append(row, strconv.FormatFloat((*v)[pos], 'f', 3, 64))
						pos++
					}
					rows = append(rows, row)
				}
			}
		}
	}
	outFile := outPrefix + "_" + strconv.Itoa(nt) + ".csv"
	outDir := filepath.Dir(outFile)
	os.MkdirAll(outDir, 0777)
	f, err := os.Create(outFile)
	if err != nil {
		fmt.Println("Can't open csv file for writing")
		errorShutdown()
	}
	w := csv.NewWriter(f)
	w.WriteAll(rows)
	if err := w.Error(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
}
