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
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"math"
	"os"
	"path"
	"strconv"
	"strings"
	"sync"

	"github.com/dustin/go-humanize"
	"github.com/montanaflynn/stats"
)

func errorShutdown() {
	fmt.Println("\nExiting scram")
	os.Exit(1)
}

// SeqLoad loads 1 or more small RNA seq. read files.
// It returns a map with a read sequence as key and a meanSe struct (normalised or raw read mean and standard error) as a value.
// Little format checking is  performed.  It is required that the input file is correctly formatted.
func SeqLoad(seqFiles []string, fileType string, adapter string, minLen int, maxLen int,
	minCount float64, noNorm bool) map[string]interface{} {
	noOfFiles, srnaMaps := loadFiles(seqFiles, fileType, minLen, maxLen, minCount, noNorm, adapter)
	seqMapAllCounts, _ := compileCounts(srnaMaps, noOfFiles, minCount)
	seqMap := calcMeanSe(seqMapAllCounts, noOfFiles)
	return seqMap
}

// IndvSeqLoad loads 1 or more small RNA seq. read files.
// It returns a map with a read sequence as key and a slice of normalized or raw individual read counts as a value.
// Little format checking is  performed.  It is required that the input file is correctly formatted.
func IndvSeqLoad(seqFiles []string, fileType string, adapter string, minLen int, maxLen int,
	minCount float64, noNorm bool) (map[string]interface{}, []string) {
	noOfFiles, srnaMaps := loadFiles(seqFiles, fileType, minLen, maxLen, minCount, noNorm, adapter)
	seqMapAllCounts, loadOrder := compileCounts(srnaMaps, noOfFiles, minCount)
	return seqMapAllCounts, loadOrder
}

// loadFiles loads replicate read files into a channel - ref_name map of read / count pairs
func loadFiles(seqFiles []string, fileType string, minLen int, maxLen int, minCount float64, noNorm bool, adapter string) (int, chan map[string]map[string]float64) {
	wg := &sync.WaitGroup{}
	noOfFiles := len(seqFiles)
	wg.Add(noOfFiles)
	fileNames := make(chan string, len(seqFiles))
	for _, fileName := range seqFiles {
		fileNames <- fileName
	}
	close(fileNames)
	srnaMaps := make(chan map[string]map[string]float64, len(seqFiles))
	for a := 0; a < len(seqFiles); a++ {
		switch {
		case fileType == "cfa":
			if a == 0 {
				fmt.Println("\nSCRAM is attempting to load read files in the default collapsed FASTA format")
			}
			go loadCfaFile(fileNames, srnaMaps, minLen, maxLen, minCount, "-", noNorm, wg)
		case fileType == "clean":
			if a == 0 {
				fmt.Println("\nSCRAM is attempting to load read files in BGI clean format")
				go loadCfaFile(fileNames, srnaMaps, minLen, maxLen, minCount, " ", noNorm, wg)
			}
		case fileType == "fa":
			if a == 0 {
				fmt.Println("\nSCRAM is attempting to load read files in FASTA format")
			}
			go loadFastx(fileNames, []byte(">"), adapter, srnaMaps, minLen, maxLen, minCount, noNorm, wg)
		case fileType == "fq":
			if a == 0 {
				fmt.Println("\nSCRAM is attempting to load read files in FASTQ format")
			}
			go loadFastx(fileNames, []byte("@"), adapter, srnaMaps, minLen, maxLen, minCount, noNorm, wg)
		}
	}
	go func(cs chan map[string]map[string]float64, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(srnaMaps, wg)
	return noOfFiles, srnaMaps
}

// loadCfaFile loads a single collapsed read file and return map of read sequence as key and normalised RPMR count as value
func loadCfaFile(fileNames chan string, srnaMaps chan map[string]map[string]float64,
	minLen int, maxLen int, minCount float64, sep string, noNorm bool, wg *sync.WaitGroup) {
	srnaMap := make(map[string]float64)
	var count float64
	var totalCount float64

	fileName := <-fileNames
	f, err := os.Open(fileName)

	if err != nil {
		fmt.Println("\nCan't load collapsed read file " + fileName)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	if fileName[len(fileName)-2:] == "gz" {
		gz, err := gzip.NewReader(f)
		if err != nil {
			fmt.Println("\nCan't decompress read file " + fileName)
			os.Exit(1)
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	}
	seqNext := false
	for scanner.Scan() {
		fastaLine := scanner.Text()
		if seqNext && len(fastaLine) == 0 {
			fmt.Println("Read file format problem - blank line between header and sequence in " + fileName)
			errorShutdown()
		}
		if strings.HasPrefix(fastaLine, ">") {
			headerLine := strings.Split(fastaLine, sep)
			err := checkHeaderError(headerLine, fileName)
			if err != nil {
				fmt.Println("Read file format problem - header error in " + fileName)
				errorShutdown()
			}
			count, err = strconv.ParseFloat(strings.Split(fastaLine, sep)[1], 32)
			seqNext = true
		} else if count >= minCount && len(fastaLine) >= minLen && len(fastaLine) <= maxLen {
			totalCount += count
			srnaMap[strings.ToUpper(fastaLine)] = count
			seqNext = false
		}
	}
	if !noNorm {
		srnaMap = rpmrNormalize(srnaMap, totalCount)
	}
	finalMap := map[string]map[string]float64{fileName: srnaMap}
	srnaMaps <- finalMap

	fmt.Println(fileName + " - " + humanize.Comma(int64(totalCount)) + " reads processed")
	wg.Done()
}

// loadFastx loads a single FASTA or FASTQ file and return map of read sequence as key and normalised RPMR count as
// value. Trim adapter from 3' end using up to 12 nt of 5' end of adapter as seed if required
func loadFastx(fileNames chan string, firstChar []byte, adapter string, srnaMaps chan map[string]map[string]float64,
	minLen int, maxLen int, minCount float64, noNorm bool, wg *sync.WaitGroup) {
	trim := false
	var seed string
	trim, seed = trimAdapter(adapter, trim, seed)

	srnaMap := make(map[string]float64)
	var totalCount float64
	fileName := <-fileNames
	f, err := os.Open(fileName)
	if err != nil {
		fmt.Println("\nCan't load read file " + fileName)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	if fileName[len(fileName)-2:] == "gz" {
		gz, err := gzip.NewReader(f)
		if err != nil {
			fmt.Println("\nCan't decompress read file " + fileName)
			os.Exit(1)
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	}

	seqNext := false
	for scanner.Scan() {
		fasta_line := scanner.Bytes()
		switch {
		case bytes.Equal(fasta_line[:1], firstChar):
			seqNext = true
		case seqNext && trim == true:
			srnaMap, totalCount, seqNext = addTrimmedRead(fasta_line, seed, minLen,
				maxLen, srnaMap, totalCount, seqNext)
		case seqNext && len(fasta_line) >= minLen && len(fasta_line) <= maxLen:
			srnaMap, totalCount, seqNext = addFullLengthRead(srnaMap, fasta_line, totalCount, seqNext)
		}
	}
	srnaMap, totalCount = removeReadsBelowMin(minCount, srnaMap, totalCount)

	if !noNorm {
		srnaMap = rpmrNormalize(srnaMap, totalCount)
	}
	final_map := map[string]map[string]float64{fileName: srnaMap}
	srnaMaps <- final_map
	fmt.Println(fileName + " - " + humanize.Comma(int64(totalCount)) + " reads processed")
	wg.Done()
}

// If adapter present, set seed for trimming (5' up to 12 nt of adapter sequence)
func trimAdapter(adapter string, trim bool, seed string) (bool, string) {
	if adapter != "nil" {
		trim = true
		switch {
		case len(adapter) < 12:
			seed = adapter
		default:
			seed = adapter[:11]
		}
	}
	return trim, seed
}

// Trim read and add it to the srna_map
func addTrimmedRead(fasta_line []byte, seed string, min_len int, max_len int, srna_map map[string]float64,
	total_count float64, seq_next bool) (map[string]float64, float64, bool) {
	read_slice := bytes.Split(fasta_line, []byte(seed))
	if len(read_slice) == 2 && len(read_slice[0]) >= min_len && len(read_slice[0]) <= max_len {
		if srna_count, ok := srna_map[string(read_slice[0])]; ok {
			srna_map[string(read_slice[0])] = srna_count + 1.0
			total_count += 1.0
		} else {
			srna_map[string(read_slice[0])] = 1.0
			total_count += 1.0
		}
	}
	seq_next = false
	return srna_map, total_count, seq_next
}

// Add full-length read to the srna map
func addFullLengthRead(srna_map map[string]float64, fasta_line []byte, total_count float64,
	seq_next bool) (map[string]float64, float64, bool) {
	if srna_count, ok := srna_map[string(fasta_line)]; ok {
		srna_map[string(fasta_line)] = srna_count + 1.0
		total_count += 1.0
	} else {
		srna_map[string(fasta_line)] = 1.0
		total_count += 1.0
	}
	seq_next = false
	return srna_map, total_count, seq_next
}

// Remove reads with count below the stated minimum for the srna_map
func removeReadsBelowMin(minCount float64, srnaMap map[string]float64, totalCount float64) (map[string]float64,
	float64) {
	if minCount > 1 {
		for srna, srnaCount := range srnaMap {
			if srnaCount < minCount {
				delete(srnaMap, srna)
				totalCount -= srnaCount
			}
		}
	}
	return srnaMap, totalCount
}

// Reads per million reads normalization of an input read library
func rpmrNormalize(srnaMap map[string]float64, total_count float64) map[string]float64 {
	for srna, srnaCount := range srnaMap {
		srnaMap[srna] = 1000000 * srnaCount / total_count
	}
	return srnaMap
}

// Checks for error in collapsed fasta header
func checkHeaderError(headerLine []string, file_name string) error {
	if len(headerLine) < 2 || len(headerLine) > 2 {
		return errors.New("\n" + file_name + " is incorrectly formatted")
	}
	return nil
}

// Compile_counts generates a map wit read seq as key and a slice of normalised counts for each read file
func compileCounts(srna_maps chan map[string]map[string]float64, no_of_files int, min_count float64) (map[string]interface{}, []string) {
	// map [srna:[count1,count2....], ...]
	seq_map_all_counts := make(map[string]interface{})
	var load_order []string
	pos := 0
	for singleSeqMap := range srna_maps {
		for file, seqMap := range singleSeqMap {
			load_order = append(load_order, path.Base(file))
			for srna, count := range seqMap {
				if _, ok := seq_map_all_counts[srna]; ok {
					// a:= append(*seq_map_all_counts[srna].(*[]float64), count)
					a := seq_map_all_counts[srna].(*[]float64)
					(*a)[pos] = count
					seq_map_all_counts[srna] = a
				} else {
					a := make([]float64, no_of_files)
					a[pos] = count
					seq_map_all_counts[srna] = &a

				}
			}
			pos++
		}
	}
	if min_count > 1 {
		removeUnderMinCount(seq_map_all_counts)
	}
	return seq_map_all_counts, load_order
}

// Remove read if its count is under the specified minimum
func removeUnderMinCount(seq_map_all_counts map[string]interface{}) {
	for srna := range seq_map_all_counts {
		counts := seq_map_all_counts[srna]
		for _, i := range *counts.(*[]float64) {
			// If using a min_count > 1, unless the srna is present in all libraries, it's removed so as not
			// to generate spurious means and standard errors
			if i == 0.0 {
				delete(seq_map_all_counts, srna)
			}
		}
	}
}

// meanSe is a struct comprising a normalised mean and standard error for a read
type meanSe struct {
	Mean float64
	Se   float64
}

// calcMeanSe calculates the mean and standard error for each slice of counts
func calcMeanSe(seq_map_all_counts map[string]interface{}, no_of_files int) map[string]interface{} {
	seqMap := make(map[string]interface{})
	sqrt := math.Sqrt(float64(no_of_files))
	for srna := range seq_map_all_counts {
		counts := seq_map_all_counts[srna]
		switch {
		case no_of_files > 1:
			countsMean, _ := stats.Mean(*counts.(*[]float64))
			countsStdDev, _ := stats.StandardDeviationSample(*counts.(*[]float64))
			seqMap[srna] = &meanSe{countsMean, countsStdDev / sqrt}
		default:
			seqMap[srna] = &meanSe{(*counts.(*[]float64))[0], 0.0}
		}
	}
	return seqMap
}

// HeaderRef is a struct comprising a reference sequence header, seques and reverse complement
type HeaderRef struct {
	Header     string
	Seq        string
	ReverseSeq string
}

// RefLoad loads a reference sequence DNA file (FASTA format).
// It returns a slice of HeaderRef structs (individual reference header, sequence and reverse complement).
func RefLoad(refFile string) []*HeaderRef {
	var totalLength int
	var refSlice []*HeaderRef
	var singleHeaderRef *HeaderRef
	var header string
	var refSeq bytes.Buffer
	f, err := os.Open(refFile)
	if err != nil {
		fmt.Println("Problem opening fasta reference file " + refFile)
		errorShutdown()
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fastaLine := scanner.Text()
		switch {
		case strings.HasPrefix(fastaLine, ">"):
			seq := refSeq.String()
			singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
			refSlice = append(refSlice, singleHeaderRef)
			header = strings.Split(fastaLine[1:], " ")[0]
			refSeq.Reset()
		case len(fastaLine) != 0:
			refSeq.WriteString(strings.ToUpper(fastaLine))
			totalLength += len(fastaLine)
		}
	}
	seq := refSeq.String()
	singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
	refSlice = append(refSlice, singleHeaderRef)
	refSlice = refSlice[1:]

	fmt.Println("No. of reference sequences: ", len(refSlice))
	fmt.Println("Combined length of reference sequences: " + humanize.Comma(int64(totalLength)) + " nt")
	return refSlice
}

// Reverse complements a DNA sequence
func reverseComplement(seq string) string {
	complement := map[rune]rune{
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'N': 'N',
	}
	runes := []rune(seq)
	var result bytes.Buffer
	for i := len(runes) - 1; i >= 0; i-- {
		result.WriteRune(complement[runes[i]])
	}
	return result.String()
}
