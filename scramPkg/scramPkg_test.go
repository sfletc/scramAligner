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
	"reflect"
	"testing"
)

func TestIndvSeqLoad_single(t *testing.T) {
	seqFiles := []string{"./test_data/test_seq_1.fa"}
	testSeq, colNames := IndvSeqLoad(seqFiles, "cfa", "nil", 18, 32, 2.0, false)
	shouldBe := map[string][]float64{
		"AAAAAAAAAAAAAAAAAAAAAAAA": {500000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGG": {250000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGC": {250000.0},
	}
	if !reflect.DeepEqual(testSeq, shouldBe) {
		t.Error("IndvSeqLoad not working for single cfa file")
	}
	if !reflect.DeepEqual(colNames, []string{"test_seq_1.fa"}) {
		t.Error("IndvSeqLoad column names incorrect")
	}
}

func TestIndvSeqLoad_clean(t *testing.T) {
	seqFiles := []string{"./test_data/test_seq_5.fa"}
	testSeq, _ := IndvSeqLoad(seqFiles, "clean", "nil", 18, 32, 2.0, false)
	shouldBe := map[string][]float64{
		"AAAAAAAAAAAAAAAAAAAAAAAA": {500000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGG": {250000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGC": {250000.0},
	}
	if !reflect.DeepEqual(testSeq, shouldBe) {
		t.Error("IndvSeqLoad not working for single clean file")
	}
}

// Loading more than one clean-format file used to hang forever (only the first
// file's goroutine was ever launched); this guards the fix.
func TestIndvSeqLoad_clean_multi(t *testing.T) {
	seqFiles := []string{"./test_data/test_seq_5.fa", "./test_data/test_seq_5.fa"}
	testSeq, colNames := IndvSeqLoad(seqFiles, "clean", "nil", 18, 32, 2.0, false)
	shouldBe := map[string][]float64{
		"AAAAAAAAAAAAAAAAAAAAAAAA": {500000.0, 500000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGG": {250000.0, 250000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGC": {250000.0, 250000.0},
	}
	if !reflect.DeepEqual(testSeq, shouldBe) {
		t.Error("IndvSeqLoad not working for multiple clean files")
	}
	// identical base names fall back to full paths so CSV columns stay distinguishable
	if !reflect.DeepEqual(colNames, seqFiles) {
		t.Error("IndvSeqLoad column names not de-duplicated for identical base names")
	}
}

func TestIndvSeqLoad_fasta(t *testing.T) {
	seqFiles := []string{"./test_data/test_fasta.fasta"}
	testSeq, _ := IndvSeqLoad(seqFiles, "fa", "nil", 18, 32, 2.0, false)
	shouldBe := map[string][]float64{
		"AAAAAAAAAAAAAAAAAAAAAAAAA": {500000.0},
		"TAAAAAAAAAAAAAAAAAAAAAAAA": {500000.0},
	}
	if !reflect.DeepEqual(testSeq, shouldBe) {
		t.Error("IndvSeqLoad not working for single fasta file")
	}
}

func TestIndvSeqLoad_multi(t *testing.T) {
	seqFiles := []string{"./test_data/test_seq_1.fa", "./test_data/test_seq_2.fa"}
	testSeq, _ := IndvSeqLoad(seqFiles, "cfa", "nil", 18, 32, 1.0, false)
	shouldBe := map[string][]float64{
		"AAAAAAAAAAAAAAAAAAAAAAAA": {500000.0, 250000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGG": {250000.0, 500000.0},
		"GGGGGGGGGGGGGGGGGGGGGGGC": {250000.0, 250000.0},
	}
	if !reflect.DeepEqual(testSeq, shouldBe) {
		t.Error("IndvSeqLoad not working for multiple cfa files")
	}
}

func TestIndvSeqLoad_multi_minCount(t *testing.T) {
	seqFiles := []string{"./test_data/test_seq_1.fa", "./test_data/test_seq_4.fa"}
	testSeq, _ := IndvSeqLoad(seqFiles, "cfa", "nil", 18, 32, 26, false)
	// reads under minCount are dropped per-file; reads absent from any file are then removed
	shouldBe := map[string][]float64{
		"AAAAAAAAAAAAAAAAAAAAAAAA": {1000000.0, 500000.0},
	}
	if !reflect.DeepEqual(testSeq, shouldBe) {
		t.Error("IndvSeqLoad not working with minCount filter")
	}
}

func TestRefLoad(t *testing.T) {
	testRef := RefLoad("./test_data/test_ref.fa")
	shouldBe := []*HeaderRef{
		{"ref_1", "AAAAAAAAAAAAAAAAAAAAAAAAA"},
		{"ref_2", "GGGGGGGGGGGGGGGGGGGGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAG"},
		{"ref_3", ""},
	}
	if len(testRef) != len(shouldBe) {
		t.Fatal("Wrong no of refs in test_ref.fa")
	}
	for i := range testRef {
		if testRef[i].Header != shouldBe[i].Header {
			t.Error("Headers don't match")
		}
		if testRef[i].Seq != shouldBe[i].Seq {
			t.Errorf("Seqs don't match: got %v, want %v", testRef[i].Seq, shouldBe[i].Seq)
		}
	}
}

func TestAlign(t *testing.T) {
	testRef := RefLoad("./test_data/test_ref_align.fa")
	seqFiles := []string{"./test_data/test_seq_1.fa"}
	testSeq, _ := IndvSeqLoad(seqFiles, "cfa", "nil", 18, 32, 1.0, false)
	testAlign := AlignReads(testSeq, testRef, 24)
	shouldBe := map[string]map[string][]int{
		"ref_1": {"AAAAAAAAAAAAAAAAAAAAAAAA": {1, 2}},
		"ref_2": {"GGGGGGGGGGGGGGGGGGGGGGGG": {1}, "AAAAAAAAAAAAAAAAAAAAAAAA": {26}},
		"ref_3": {"AAAAAAAAAAAAAAAAAAAAAAAA": {-2, -1}},
	}
	if !reflect.DeepEqual(testAlign, shouldBe) {
		t.Error("Alignments are not equal")
	}
}

func TestProfileAlignIndv(t *testing.T) {
	testRef := RefLoad("./test_data/test_ref_align.fa")
	seqFiles := []string{"./test_data/test_seq_1.fa"}
	testSeq, _ := IndvSeqLoad(seqFiles, "cfa", "nil", 18, 32, 1.0, false)
	testAlign := AlignReads(testSeq, testRef, 24)

	testAlignSplit := ProfileSplit(testAlign, testSeq)
	shouldBeSplit := map[string]*singleAlignments{
		"ref_1": {
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "+", []float64{100000}},
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "+", []float64{100000}},
		},
		"ref_2": {
			&singleAlignment{"GGGGGGGGGGGGGGGGGGGGGGGG", 1, 1, "+", []float64{250000}},
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 26, "+", []float64{100000}},
		},
		"ref_3": {
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "-", []float64{100000}},
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "-", []float64{100000}},
		},
	}
	if !reflect.DeepEqual(shouldBeSplit, testAlignSplit) {
		t.Error("Profile split alignment is incorrect")
	}

	testAlignNoSplit := ProfileNoSplit(testAlign, testSeq)
	shouldBeNoSplit := map[string]*singleAlignments{
		"ref_1": {
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "+", []float64{500000}},
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "+", []float64{500000}},
		},
		"ref_2": {
			&singleAlignment{"GGGGGGGGGGGGGGGGGGGGGGGG", 1, 1, "+", []float64{250000}},
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 26, "+", []float64{500000}},
		},
		"ref_3": {
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "-", []float64{500000}},
			&singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "-", []float64{500000}},
		},
	}
	if !reflect.DeepEqual(shouldBeNoSplit, testAlignNoSplit) {
		t.Error("Profile no split alignment is incorrect")
	}
}
