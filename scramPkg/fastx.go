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
	"errors"
	"fmt"
	"io"
	"os"
	"path"
	"strconv"
	"strings"
	"sync"

	"github.com/dustin/go-humanize"
	"github.com/klauspost/compress/gzip"
)

func errorShutdown() {
	fmt.Println("\nExiting scram")
	os.Exit(1)
}

// baseNormalize maps each sequence byte to its canonical form: lowercase to
// uppercase and U/u to T, so RNA and soft-masked input match the reference.
var baseNormalize = func() [256]byte {
	var t [256]byte
	for i := 0; i < 256; i++ {
		t[i] = byte(i)
	}
	for c := byte('a'); c <= 'z'; c++ {
		t[c] = c - 'a' + 'A'
	}
	t['U'] = 'T'
	t['u'] = 'T'
	return t
}()

func normalizeSeqInPlace(seq []byte) {
	for i, c := range seq {
		seq[i] = baseNormalize[c]
	}
}

// fileJob identifies one input read file and its column slot in the output.
type fileJob struct {
	idx  int
	name string
}

// fileResult carries one loaded read file's counts back with its column slot,
// so output column order matches the -f argument order regardless of which
// goroutine finishes first.
type fileResult struct {
	idx    int
	counts map[string]float64
}

// openFastxScanner opens a (possibly gzipped) FASTA/FASTQ file and returns a
// scanner that can handle single-line chromosome-scale sequences, plus a close
// function. The default bufio.Scanner token limit of 64KB would silently
// truncate unwrapped reference FASTA otherwise.
func openFastxScanner(fileName string) (*bufio.Scanner, func(), error) {
	f, err := os.Open(fileName)
	if err != nil {
		return nil, nil, err
	}
	var r io.Reader = bufio.NewReaderSize(f, 1<<20)
	closeFn := func() { f.Close() }
	if strings.HasSuffix(fileName, ".gz") {
		gz, err := gzip.NewReader(r)
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		r = gz
		closeFn = func() {
			gz.Close()
			f.Close()
		}
	}
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 1<<20), 1<<30)
	return scanner, closeFn, nil
}

func checkScannerError(scanner *bufio.Scanner, fileName string) {
	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading " + fileName + ": " + err.Error())
		errorShutdown()
	}
}

// scanFastaRecords joins wrapped sequence lines before visiting each record.
// The sequence buffer is reused after visit returns.
func scanFastaRecords(scanner *bufio.Scanner, visit func(string, []byte) error) error {
	var header string
	var seq []byte
	haveHeader := false
	flush := func() error {
		if !haveHeader {
			return nil
		}
		if len(seq) == 0 {
			return fmt.Errorf("missing sequence for header %q", header)
		}
		return visit(header, seq)
	}
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			if err := flush(); err != nil {
				return err
			}
			header = string(line)
			seq = seq[:0]
			haveHeader = true
		} else {
			if !haveHeader {
				return fmt.Errorf("sequence before the first FASTA header")
			}
			seq = append(seq, line...)
		}
	}
	if err := scanner.Err(); err != nil {
		return err
	}
	return flush()
}

// IndvSeqLoad loads 1 or more small RNA seq. read files.
// It returns a map with a read sequence as key and a slice of normalized or raw individual read counts as a value,
// plus the CSV column names for the read files in -f argument order.
// Little format checking is performed.  It is required that the input file is correctly formatted.
func IndvSeqLoad(seqFiles []string, fileType string, adapter string, minLen int, maxLen int,
	minCount float64, noNorm bool) (map[string][]float64, []string) {
	srnaMaps := loadFiles(seqFiles, fileType, minLen, maxLen, minCount, noNorm, adapter)
	seqMapAllCounts := compileCounts(srnaMaps, len(seqFiles), minCount)
	return seqMapAllCounts, columnNames(seqFiles)
}

// loadFiles loads replicate read files concurrently, one goroutine per file.
func loadFiles(seqFiles []string, fileType string, minLen int, maxLen int, minCount float64, noNorm bool, adapter string) chan fileResult {
	wg := &sync.WaitGroup{}
	wg.Add(len(seqFiles))
	srnaMaps := make(chan fileResult, len(seqFiles))
	switch fileType {
	case "cfa":
		fmt.Println("\nSCRAM is attempting to load read files in the default collapsed FASTA format")
	case "clean":
		fmt.Println("\nSCRAM is attempting to load read files in BGI clean format")
	case "fa":
		fmt.Println("\nSCRAM is attempting to load read files in FASTA format")
	case "fq":
		fmt.Println("\nSCRAM is attempting to load read files in FASTQ format")
	}
	for i, fileName := range seqFiles {
		job := fileJob{i, fileName}
		switch fileType {
		case "cfa":
			go loadCfaFile(job, srnaMaps, minLen, maxLen, minCount, "-", noNorm, wg)
		case "clean":
			go loadCfaFile(job, srnaMaps, minLen, maxLen, minCount, " ", noNorm, wg)
		case "fa":
			go loadFastx(job, '>', adapter, srnaMaps, minLen, maxLen, minCount, noNorm, wg)
		case "fq":
			go loadFastx(job, '@', adapter, srnaMaps, minLen, maxLen, minCount, noNorm, wg)
		}
	}
	go func() {
		wg.Wait()
		close(srnaMaps)
	}()
	return srnaMaps
}

// loadCfaFile loads a single collapsed read file and sends a map of read sequence to normalised RPMR count.
func loadCfaFile(job fileJob, srnaMaps chan fileResult,
	minLen int, maxLen int, minCount float64, sep string, noNorm bool, wg *sync.WaitGroup) {
	defer wg.Done()
	srnaMap := make(map[string]float64)
	var totalCount float64

	scanner, closeFile, err := openFastxScanner(job.name)
	if err != nil {
		fmt.Println("\nCan't load collapsed read file " + job.name + ": " + err.Error())
		os.Exit(1)
	}
	defer closeFile()
	err = scanFastaRecords(scanner, func(header string, seq []byte) error {
		headerLine := strings.Split(header, sep)
		if err := checkHeaderError(headerLine, job.name); err != nil {
			return err
		}
		count, err := strconv.ParseFloat(headerLine[1], 64)
		if err != nil {
			return fmt.Errorf("can't parse read count in header %q: %w", header, err)
		}
		if count > 0 && len(seq) >= minLen && len(seq) <= maxLen {
			normalizeSeqInPlace(seq)
			totalCount += count
			srnaMap[string(seq)] += count
		}
		return nil
	})
	if err != nil {
		fmt.Printf("Read file format problem in %s: %v\n", job.name, err)
		errorShutdown()
	}
	// Apply the threshold to each sequence's combined count, as for FASTA/FASTQ.
	srnaMap, totalCount = removeReadsBelowMin(minCount, srnaMap, totalCount)
	if totalCount == 0 {
		fmt.Println("Warning: no reads in " + job.name + " passed the length/count filters")
	}
	if !noNorm {
		srnaMap = rpmrNormalize(srnaMap, totalCount)
	}
	srnaMaps <- fileResult{job.idx, srnaMap}

	fmt.Println(job.name + " - " + humanize.Comma(int64(totalCount)) + " reads processed")
}

// loadFastx loads a single FASTA or FASTQ file and sends a map of read sequence to normalised RPMR count.
// Trims the adapter from the 3' end using up to 12 nt of the 5' end of the adapter as seed if required.
func loadFastx(job fileJob, firstChar byte, adapter string, srnaMaps chan fileResult,
	minLen int, maxLen int, minCount float64, noNorm bool, wg *sync.WaitGroup) {
	defer wg.Done()
	trim, seed := trimAdapter(adapter)

	srnaMap := make(map[string]float64)
	var totalCount float64
	var untrimmed int64

	scanner, closeFile, err := openFastxScanner(job.name)
	if err != nil {
		fmt.Println("\nCan't load read file " + job.name + ": " + err.Error())
		os.Exit(1)
	}
	defer closeFile()

	addRead := func(line []byte) {
		normalizeSeqInPlace(line)
		if trim {
			i := bytes.Index(line, seed)
			if i < 0 {
				untrimmed++
				return
			}
			line = line[:i]
		}
		if len(line) >= minLen && len(line) <= maxLen {
			srnaMap[string(line)]++
			totalCount++
		}
	}

	if firstChar == '@' {
		// FASTQ: strict 4-line records. Sniffing line types by first character is
		// unsafe ('+' lines may repeat the header, quality lines may start with '@').
		lineNo := 0
		var seq []byte
		for scanner.Scan() {
			line := scanner.Bytes()
			switch lineNo % 4 {
			case 0:
				if len(line) == 0 || line[0] != '@' {
					fmt.Printf("Read file format problem - expected a FASTQ header at line %d in %s\n", lineNo+1, job.name)
					errorShutdown()
				}
			case 1:
				// Scanner's buffer may change before the quality line is read.
				seq = append(seq[:0], line...)
			case 2:
				if len(line) == 0 || line[0] != '+' {
					fmt.Printf("Read file format problem - expected a FASTQ '+' separator at line %d in %s\n", lineNo+1, job.name)
					errorShutdown()
				}
			case 3:
				if len(line) != len(seq) {
					fmt.Printf("Read file format problem - FASTQ quality length %d does not match sequence length %d at line %d in %s\n",
						len(line), len(seq), lineNo+1, job.name)
					errorShutdown()
				}
				addRead(seq)
			}
			lineNo++
		}
		checkScannerError(scanner, job.name)
		if lineNo%4 != 0 {
			fmt.Printf("Read file format problem - incomplete FASTQ record at end of %s after line %d\n", job.name, lineNo)
			errorShutdown()
		}
	} else {
		err := scanFastaRecords(scanner, func(_ string, seq []byte) error {
			addRead(seq)
			return nil
		})
		if err != nil {
			fmt.Printf("Read file format problem in %s: %v\n", job.name, err)
			errorShutdown()
		}
	}
	if untrimmed > 0 {
		fmt.Println(job.name + " - " + humanize.Comma(untrimmed) +
			" reads discarded (adapter seed not found)")
	}
	srnaMap, totalCount = removeReadsBelowMin(minCount, srnaMap, totalCount)
	if totalCount == 0 {
		fmt.Println("Warning: no reads in " + job.name + " passed the length/count filters")
	}

	if !noNorm {
		srnaMap = rpmrNormalize(srnaMap, totalCount)
	}
	srnaMaps <- fileResult{job.idx, srnaMap}
	fmt.Println(job.name + " - " + humanize.Comma(int64(totalCount)) + " reads processed")
}

// trimAdapter returns whether to trim and the seed to trim with
// (5' up to 12 nt of the adapter sequence, normalized like read sequences).
func trimAdapter(adapter string) (bool, []byte) {
	if adapter == "nil" || adapter == "" {
		return false, nil
	}
	seed := []byte(adapter)
	if len(seed) > 12 {
		seed = seed[:12]
	}
	normalizeSeqInPlace(seed)
	return true, seed
}

// removeReadsBelowMin removes reads with a count below the stated minimum from the srnaMap
func removeReadsBelowMin(minCount float64, srnaMap map[string]float64, totalCount float64) (map[string]float64,
	float64) {
	for srna, srnaCount := range srnaMap {
		if srnaCount < minCount {
			delete(srnaMap, srna)
			totalCount -= srnaCount
		}
	}
	return srnaMap, totalCount
}

// rpmrNormalize applies reads per million reads normalization to an input read library
func rpmrNormalize(srnaMap map[string]float64, totalCount float64) map[string]float64 {
	if totalCount == 0 {
		return srnaMap
	}
	for srna, srnaCount := range srnaMap {
		srnaMap[srna] = 1000000 * srnaCount / totalCount
	}
	return srnaMap
}

// checkHeaderError checks for errors in a collapsed fasta header
func checkHeaderError(headerLine []string, fileName string) error {
	if len(headerLine) != 2 {
		return errors.New("\n" + fileName + " is incorrectly formatted")
	}
	return nil
}

// compileCounts generates a map with read seq as key and a slice of normalised counts (one slot per
// read file, in -f argument order) as value.
func compileCounts(srnaMaps chan fileResult, noOfFiles int, minCount float64) map[string][]float64 {
	seqMapAllCounts := make(map[string][]float64)
	for res := range srnaMaps {
		for srna, count := range res.counts {
			counts, ok := seqMapAllCounts[srna]
			if !ok {
				counts = make([]float64, noOfFiles)
				seqMapAllCounts[srna] = counts
			}
			counts[res.idx] = count
		}
	}
	if minCount > 1 {
		removeUnderMinCount(seqMapAllCounts)
	}
	return seqMapAllCounts
}

// columnNames derives the CSV column names from the read file paths: the base name,
// or the full path where base names collide.
func columnNames(seqFiles []string) []string {
	baseCount := make(map[string]int)
	for _, f := range seqFiles {
		baseCount[path.Base(f)]++
	}
	names := make([]string, len(seqFiles))
	for i, f := range seqFiles {
		if baseCount[path.Base(f)] > 1 {
			names[i] = f
		} else {
			names[i] = path.Base(f)
		}
	}
	return names
}

// removeUnderMinCount removes a read unless it is present in all libraries, so as not
// to generate spurious means and standard errors when using a minCount > 1.
func removeUnderMinCount(seqMapAllCounts map[string][]float64) {
	removed := 0
	for srna, counts := range seqMapAllCounts {
		for _, count := range counts {
			if count == 0.0 {
				delete(seqMapAllCounts, srna)
				removed++
				break
			}
		}
	}
	if removed > 0 {
		fmt.Printf("%d reads removed as not present in all libraries (--minCount > 1)\n", removed)
	}
}

// HeaderRef is a struct comprising a reference sequence header and sequence.
// The reverse complement is not stored; reverse-strand matches are computed
// from rolling 2-bit windows during alignment.
type HeaderRef struct {
	Header string
	Seq    string
}

// RefLoad loads a reference sequence DNA file (FASTA format).
// It returns a slice of HeaderRef structs (individual reference header and sequence).
func RefLoad(refFile string) []*HeaderRef {
	var totalLength int
	var refSlice []*HeaderRef
	var header string
	haveHeader := false
	seenHeaders := make(map[string]bool)
	var refSeq bytes.Buffer

	scanner, closeFile, err := openFastxScanner(refFile)
	if err != nil {
		fmt.Println("Problem opening fasta reference file " + refFile + ": " + err.Error())
		errorShutdown()
	}
	defer closeFile()

	flush := func() {
		if !haveHeader {
			return
		}
		refSlice = append(refSlice, &HeaderRef{header, refSeq.String()})
		refSeq.Reset()
	}
	for scanner.Scan() {
		line := scanner.Bytes()
		switch {
		case len(line) > 0 && line[0] == '>':
			flush()
			header = firstField(string(line[1:]))
			if seenHeaders[header] {
				fmt.Println("Duplicate reference header (truncated at whitespace): " + header)
				errorShutdown()
			}
			seenHeaders[header] = true
			haveHeader = true
		case len(line) != 0:
			normalizeSeqInPlace(line)
			refSeq.Write(line)
			totalLength += len(line)
		}
	}
	checkScannerError(scanner, refFile)
	flush()
	if len(refSlice) == 0 {
		fmt.Println("No reference sequences loaded from " + refFile + " - is it FASTA format?")
		errorShutdown()
	}

	fmt.Println("No. of reference sequences: ", len(refSlice))
	fmt.Println("Combined length of reference sequences: " + humanize.Comma(int64(totalLength)) + " nt")
	return refSlice
}

// firstField returns s truncated at the first space or tab.
func firstField(s string) string {
	if i := strings.IndexAny(s, " \t"); i >= 0 {
		return s[:i]
	}
	return s
}
