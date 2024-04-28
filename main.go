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

package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
	"time"

	scramPkg "scramAligner/scramPkg"

	"github.com/spf13/cobra"
)

const version = "2.1.0"

var cfgFile string
var fastaSet1 string
var readFileType string
var length string
var alignTo string
var outFilePrefix string
var noSplit bool
var minLen int
var maxLen int
var minCount float64
var adapter string
var noNorm bool

var rootCmd = &cobra.Command{
	Use:   "scramAligner -r <reference> -f <fileSet> -l <length> -o <outFilePrefix>",
	Short: "Exact match small RNA read alignment v" + version,
	Long:  `Exact match small RNA read alignment v` + version,
	Run:   align,
}

func main() {
	cobra.CheckErr(rootCmd.Execute())
}

func init() {

	rootCmd.PersistentFlags().StringVarP(&alignTo, "alignTo", "r", "", "path/to/FASTA reference file - REQUIRED")
	rootCmd.MarkFlagRequired("alignTo")
	rootCmd.PersistentFlags().StringVarP(&fastaSet1, "fastxSet", "f", "", "comma-separated path/to/read file set 1. GZIPped files must have .gz file extension - REQUIRED")
	rootCmd.MarkFlagRequired("fastxSet")
	rootCmd.PersistentFlags().StringVarP(&readFileType, "readFileType", "t", "fq", "Read file type: cfa (collapsed FASTA), fa (FASTA), fq (FASTQ), clean (BGI clean.fa).")
	rootCmd.PersistentFlags().StringVarP(&length, "length", "l", "", "comma-separated read (sRNA) lengths to align - REQUIRED")
	rootCmd.MarkFlagRequired("length")
	rootCmd.PersistentFlags().StringVarP(&outFilePrefix, "outFilePrefix", "o", "", "path/to/outfile prefix (len.csv will be appended) - REQUIRED")
	rootCmd.MarkFlagRequired("outFilePrefix")
	rootCmd.PersistentFlags().BoolVar(&noSplit, "noSplit", false, "Do not split alignment count for each read by the number of times it aligns")
	rootCmd.PersistentFlags().BoolVar(&noNorm, "noNorm", false, "Do not normalize read counts by library size (i.e. reads per million reads)")
	rootCmd.PersistentFlags().IntVar(&minLen, "minLen", 18, "Minimum read length to include for RPMR normalization")
	rootCmd.PersistentFlags().IntVar(&maxLen, "maxLen", 32, "Maximum read length to include for RPMR normalization")
	rootCmd.PersistentFlags().Float64Var(&minCount, "minCount", 1.0, "Minimum read count for alignment and to include for RPMR normalization")
	rootCmd.PersistentFlags().StringVar(&adapter, "adapter", "nil", "3' adapter sequence to trim - FASTA & FASTQ only")

}

func align(cmd *cobra.Command, args []string) {
	if readFileType != "cfa" && readFileType != "fa" && readFileType != "fq" && readFileType != "clean" {
		fmt.Println("\nCan't parse read file type " + readFileType)
		os.Exit(1)
	}
	t0 := time.Now()
	var a map[string]interface{}
	var fileOrder []string

	fmt.Println("\nLoading individual read counts")
	a, fileOrder = scramPkg.IndvSeqLoad(strings.Split(fastaSet1, ","), readFileType, adapter, minLen, maxLen, minCount, noNorm)

	fmt.Println("\nLoading reference")
	c := scramPkg.RefLoad(alignTo)
	for _, nt := range strings.Split(length, ",") {
		nt, _ := strconv.Atoi(nt)
		fmt.Printf("\nAligning %v nt reads\n", nt)
		d := scramPkg.AlignReads(a, c, nt)
		switch {
		case !noSplit:
			e := scramPkg.ProfileSplit(d, a)
			scramPkg.ProfileToCsv(e, c, nt, outFilePrefix, fileOrder)
		default:
			e := scramPkg.ProfileNoSplit(d, a)
			scramPkg.ProfileToCsv(e, c, nt, outFilePrefix, fileOrder)

		}
	}
	t1 := time.Now()
	fmt.Printf("\nAlignment complete.  Total time taken = %s\n", t1.Sub(t0))
}
