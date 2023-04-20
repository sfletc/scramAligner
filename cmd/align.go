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

package cmd

import (
	"fmt"
	"os"
	scramPkg "scramAligner/scramPkg"
	"strconv"
	"strings"
	"time"

	"github.com/spf13/cobra"
)

// profileCmd represents the align command

var alignCmd = &cobra.Command{
	Use:   "align",
	Short: "Align reads of length l from 1 read file set to all sequences in a reference file",
	Long: `Align reads of length l from 1 read file set to all sequences in a reference file

For example:

scramAligner align -r ref.fa -1 seq1a.fa,seq1b.fa,seq1c.fa -l 21,22,24 -o testAlign

`,
	Run: func(cmd *cobra.Command, args []string) {
		if readFileType != "cfa" && readFileType != "fa" && readFileType != "fq" && readFileType != "clean" {
			fmt.Println("\nCan't parse read file type " + readFileType)
			os.Exit(1)
		}
		t0 := time.Now()
		var a map[string]interface{}
		var fileOrder []string
		switch {
		case !indv:
			fmt.Println("\nLoading mean and standard errors of replicate reads")
			a = scramPkg.SeqLoad(strings.Split(fastaSet1, ","), readFileType, adapter, minLen, maxLen, minCount, noNorm)
		case indv:
			fmt.Println("\nLoading individual read counts")
			a, fileOrder = scramPkg.IndvSeqLoad(strings.Split(fastaSet1, ","), readFileType, adapter, minLen, maxLen, minCount, noNorm)
		}
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
	},
}

func init() {
	RootCmd.AddCommand(alignCmd)

}
