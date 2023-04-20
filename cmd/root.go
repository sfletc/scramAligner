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

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

const version = "2.0.1"

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
var indv bool

var RootCmd = &cobra.Command{
	Use:   "scram2",
	Short: "Exact match small RNA read alignment v" + version,
	Long:  `Exact match small RNA read alignment v` + version,
}

func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

func init() {
	cobra.OnInitialize(initConfig)
	RootCmd.PersistentFlags().StringVarP(&alignTo, "alignTo", "r", "", "path/to/FASTA reference file")
	RootCmd.PersistentFlags().StringVarP(&fastaSet1, "fastxSet1", "1", "", "comma-separated path/to/read file set 1. GZIPped files must have .gz file extension")
	RootCmd.PersistentFlags().StringVarP(&readFileType, "readFileType", "t", "fq", "Read file type: cfa (collapsed FASTA), fa (FASTA), fq (FASTQ), clean (BGI clean.fa).")
	RootCmd.PersistentFlags().StringVarP(&length, "length", "l", "", "comma-separated read (sRNA) lengths to align")
	RootCmd.PersistentFlags().StringVarP(&outFilePrefix, "outFilePrefix", "o", "", "path/to/outfile prefix (len.csv will be appended)")
	RootCmd.PersistentFlags().BoolVar(&noSplit, "noSplit", false, "Do not split alignment count for each read by the number of times it aligns")
	RootCmd.PersistentFlags().BoolVar(&noNorm, "noNorm", false, "Do not normalize read counts by library size (i.e. reads per million reads)")
	RootCmd.PersistentFlags().IntVar(&minLen, "minLen", 18, "Minimum read length to include for RPMR normalization")
	RootCmd.PersistentFlags().IntVar(&maxLen, "maxLen", 32, "Maximum read length to include for RPMR normalization")
	RootCmd.PersistentFlags().Float64Var(&minCount, "minCount", 1.0, "Minimum read count for alignment and to include for RPMR normalization")
	RootCmd.PersistentFlags().StringVar(&adapter, "adapter", "nil", "3' adapter sequence to trim - FASTA & FASTQ only")
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {
	if cfgFile != "" { // enable ability to specify config file via flag
		viper.SetConfigFile(cfgFile)
	}

	viper.SetConfigName(".scram") // name of config file (without extension)
	viper.AddConfigPath("$HOME")  // adding home directory as first search path
	viper.AutomaticEnv()          // read in environment variables that match

	// If a config file is found, read it in.
	if err := viper.ReadInConfig(); err == nil {
		fmt.Println("Using config file:", viper.ConfigFileUsed())
	}
}
