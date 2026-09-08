// Regression tests for fixes to reference/read loading and the k-mer alignment core.

package scramPkg

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

func TestLoadFasta_wrappedRecords(t *testing.T) {
	seq := "ACGTACGTACGTACGTACGTA"
	for _, fileType := range []string{"fa", "cfa", "clean"} {
		t.Run(fileType, func(t *testing.T) {
			header := func(id int) string {
				switch fileType {
				case "cfa":
					return fmt.Sprintf(">%d-1\n", id)
				case "clean":
					return fmt.Sprintf(">%d 1\n", id)
				default:
					return fmt.Sprintf(">%d\n", id)
				}
			}
			// The long record must not be counted as a 21-nt read, and the last
			// wrapped record must be flushed even without a final newline.
			content := header(1) + seq[:10] + "\n" + seq[10:] + "\n" +
				header(2) + seq + "\n" + seq + "\n" +
				header(3) + strings.ToLower(seq[:15]) + "\n" + seq[15:]
			path := writeTempFile(t, "wrapped.fa", content)
			got, _ := IndvSeqLoad([]string{path}, fileType, "nil", 18, 32, 1, true)
			want := map[string][]float64{seq: {2}}
			if !reflect.DeepEqual(got, want) {
				t.Fatalf("got %v, want %v", got, want)
			}
		})
	}
}

func TestLoadFastx_trimsCompleteRecords(t *testing.T) {
	seq := "ACGTACGTACGTACGTACGTA"
	adapter := "AGATCGGAAGAG"
	for _, fileType := range []string{"fa", "fq"} {
		t.Run(fileType, func(t *testing.T) {
			content := ">r1\n" + seq + adapter[:5] + "\n" + adapter[5:]
			if fileType == "fq" {
				// Quality length must be checked against the untrimmed sequence.
				content = "@r1\n" + seq + adapter + "\n+r1\n" + strings.Repeat("@", len(seq+adapter))
			}
			path := writeTempFile(t, "reads."+fileType, content)
			got, _ := IndvSeqLoad([]string{path}, fileType, adapter, 18, 32, 1, true)
			want := map[string][]float64{seq: {1}}
			if !reflect.DeepEqual(got, want) {
				t.Fatalf("got %v, want %v", got, want)
			}
		})
	}
}

func TestLoadCfa_minCountAfterAggregation(t *testing.T) {
	seq := "ACGTACGTACGTACGTACGTA"
	for _, fileType := range []string{"cfa", "clean"} {
		for _, noNorm := range []bool{true, false} {
			t.Run(fmt.Sprintf("%s/noNorm=%v", fileType, noNorm), func(t *testing.T) {
				sep := "-"
				if fileType == "clean" {
					sep = " "
				}
				content := ">1" + sep + "1\n" + seq + "\n>2" + sep + "1\n" +
					strings.ToLower(seq) + "\n>3" + sep + "1\n" + strings.Repeat("C", 21) + "\n"
				path := writeTempFile(t, "duplicate.fa", content)
				got, _ := IndvSeqLoad([]string{path}, fileType, "nil", 18, 32, 2, noNorm)
				count := 2.0
				if !noNorm {
					count = 1000000
				}
				want := map[string][]float64{seq: {count}}
				if !reflect.DeepEqual(got, want) {
					t.Fatalf("got %v, want %v", got, want)
				}
			})
		}
	}
}

func TestLoadFastq_rejectsMalformed(t *testing.T) {
	// The public loader exits on malformed input, so exercise it in a child
	// test process and assert that callers receive a nonzero exit status.
	if path := os.Getenv("SCRAM_TEST_FASTQ_FILE"); path != "" {
		IndvSeqLoad([]string{path}, "fq", "nil", 18, 32, 1, true)
		return
	}
	seq := "ACGTACGTACGTACGTACGTA"
	valid := "@r1\n" + seq + "\n+\n" + strings.Repeat("F", len(seq)) + "\n"
	for name, content := range map[string]string{
		"missing_sequence":  "@r1\n",
		"missing_separator": "@r1\n" + seq + "\n",
		"missing_quality":   "@r1\n" + seq + "\n+\n",
		"short_quality":     "@r1\n" + seq + "\n+\nFFF\n",
		"long_quality":      "@r1\n" + seq + "\n+\n" + strings.Repeat("F", len(seq)+1) + "\n",
		"truncated_second":  valid + "@r2\n" + seq + "\n",
	} {
		t.Run(name, func(t *testing.T) {
			path := writeTempFile(t, "bad.fq", content)
			cmd := exec.Command(os.Args[0], "-test.run=^TestLoadFastq_rejectsMalformed$")
			cmd.Env = append(os.Environ(), "SCRAM_TEST_FASTQ_FILE="+path)
			output, err := cmd.CombinedOutput()
			if exitErr, ok := err.(*exec.ExitError); !ok || exitErr.ExitCode() != 1 {
				t.Fatalf("expected exit status 1, got %v: %s", err, output)
			}
			if !strings.Contains(string(output), "Read file format problem") {
				t.Fatalf("missing format error: %s", output)
			}
		})
	}
}

func writeTempFile(t *testing.T, name, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

// A reference line longer than bufio.Scanner's 64KB default token limit used to be
// silently truncated.
func TestRefLoad_longSingleLine(t *testing.T) {
	seq := strings.Repeat("ACGT", 50000) // 200KB on one line
	path := writeTempFile(t, "long.fa", ">chr1\n"+seq+"\n")
	refs := RefLoad(path)
	if len(refs) != 1 {
		t.Fatalf("expected 1 reference, got %d", len(refs))
	}
	if len(refs[0].Seq) != len(seq) {
		t.Errorf("reference truncated: got %d nt, want %d", len(refs[0].Seq), len(seq))
	}
}

// Lowercase and RNA (U) bases in FASTA/FASTQ reads used to silently never align.
func TestLoadFastq_normalizesCase(t *testing.T) {
	fq := "@read1\nacguacguacguacguacguac\n+\nFFFFFFFFFFFFFFFFFFFFFF\n" +
		"@read2\nACGTACGTACGTACGTACGTAC\n+read2\nFFFFFFFFFFFFFFFFFFFFFF\n"
	path := writeTempFile(t, "reads.fq", fq)
	seqMap, _ := IndvSeqLoad([]string{path}, "fq", "nil", 18, 32, 1.0, true)
	want := map[string][]float64{"ACGTACGTACGTACGTACGTAC": {2.0}}
	if !reflect.DeepEqual(seqMap, want) {
		t.Errorf("got %v, want %v", seqMap, want)
	}
}

// A '+' separator line repeating the header, and quality strings starting with '@'
// or '+', used to confuse the first-character FASTQ sniffing.
func TestLoadFastq_trickySeparatorAndQuality(t *testing.T) {
	fq := "@read_00012345_length_22\nACGTACGTACGTACGTACGTAC\n" +
		"+read_00012345_length_22\n@@FFFFFFFFFFFFFFFFFF+@\n" +
		"@read2\nTTTTACGTACGTACGTACGTAC\n+\n+FFFFFFFFFFFFFFFFFFFF@\n"
	path := writeTempFile(t, "tricky.fq", fq)
	seqMap, _ := IndvSeqLoad([]string{path}, "fq", "nil", 18, 32, 1.0, true)
	want := map[string][]float64{
		"ACGTACGTACGTACGTACGTAC": {1.0},
		"TTTTACGTACGTACGTACGTAC": {1.0},
	}
	if !reflect.DeepEqual(seqMap, want) {
		t.Errorf("got %v, want %v", seqMap, want)
	}
}

// Duplicate sequences in a collapsed FASTA used to overwrite instead of sum,
// desynchronizing counts from the RPM denominator.
func TestLoadCfa_duplicateSequencesSum(t *testing.T) {
	cfa := ">1-10\nACGTACGTACGTACGTACGTAC\n>2-5\nacgtacgtacgtacgtacgtac\n"
	path := writeTempFile(t, "dup.fa", cfa)
	seqMap, _ := IndvSeqLoad([]string{path}, "cfa", "nil", 18, 32, 1.0, true)
	want := map[string][]float64{"ACGTACGTACGTACGTACGTAC": {15.0}}
	if !reflect.DeepEqual(seqMap, want) {
		t.Errorf("got %v, want %v", seqMap, want)
	}
}

// A blank line following a length-filtered read used to abort the whole run.
func TestLoadCfa_blankLineAfterFilteredRead(t *testing.T) {
	cfa := ">1-10\nACGT\n\n>2-5\nACGTACGTACGTACGTACGTAC\n"
	path := writeTempFile(t, "blank.fa", cfa)
	seqMap, _ := IndvSeqLoad([]string{path}, "cfa", "nil", 18, 32, 1.0, true)
	want := map[string][]float64{"ACGTACGTACGTACGTACGTAC": {5.0}}
	if !reflect.DeepEqual(seqMap, want) {
		t.Errorf("got %v, want %v", seqMap, want)
	}
}

// Adapter trimming: the first 12 nt of the adapter are the seed; reads without
// the seed are dropped (and reported), trimmed reads are counted.
func TestLoadFasta_adapterTrimming(t *testing.T) {
	adapter := "AGATCGGAAGAGCACACGTCT"
	fa := ">r1\nACGTACGTACGTACGTACGTAC" + adapter + "\n" +
		">r2\nTTTTACGTACGTACGTACGTAC" + adapter[:12] + "GARBAGE\n" +
		">r3\nCCCCACGTACGTACGTACGTAC\n" // no adapter: dropped
	path := writeTempFile(t, "adapter.fa", fa)
	seqMap, _ := IndvSeqLoad([]string{path}, "fa", adapter, 18, 32, 1.0, true)
	want := map[string][]float64{
		"ACGTACGTACGTACGTACGTAC": {1.0},
		"TTTTACGTACGTACGTACGTAC": {1.0},
	}
	if !reflect.DeepEqual(seqMap, want) {
		t.Errorf("got %v, want %v", seqMap, want)
	}
}

// AlignReadsMulti in one pass must equal AlignReads run per length.
func TestAlignReadsMulti_matchesPerLength(t *testing.T) {
	ref := ">refA\nGGCATTCAGGATCGATCGGATCGATTAGGCATCAGTCAGGACTCAGGACTGACTGATCAGTTACGGACTAGC\n" +
		">refB\nTTACGGACTAGCAAGGCATTCAGGATCGATCGGATCG\n"
	refPath := writeTempFile(t, "ref.fa", ref)
	refs := RefLoad(refPath)

	seqMap := map[string][]float64{
		"GGCATTCAGGATCGATCGGAT":    {10.0}, // 21 nt, fwd in refA and refB
		"CTAGTCCGTAACTGATCAGTCAGT": {4.0},  // 24 nt, reverse complement match in refA
		"AAAAAAAAAAAAAAAAAAAAA":    {1.0},  // 21 nt, no match
	}
	lengths := []int{21, 24}
	multi := AlignReadsMulti(IndexReads(seqMap, lengths), refs, lengths)
	for _, nt := range lengths {
		single := AlignReads(seqMap, refs, nt)
		if !reflect.DeepEqual(multi[nt], single) {
			t.Errorf("length %d: multi %v != single %v", nt, multi[nt], single)
		}
	}
	if len(multi[21]["refA"]["GGCATTCAGGATCGATCGGAT"]) == 0 {
		t.Error("expected forward alignment in refA")
	}
	if len(multi[24]["refA"]["CTAGTCCGTAACTGATCAGTCAGT"]) == 0 {
		t.Error("expected reverse alignment in refA")
	}
}

// Reference windows containing N (or any non-ACGT base) must not produce
// alignments, and reads containing N are excluded from the index.
func TestAlign_ambiguityHandling(t *testing.T) {
	ref := ">amb\nACGTACGTACGTACGTACGTACNNNNACGTACGTACGTACGTACGTAC\n"
	refPath := writeTempFile(t, "amb.fa", ref)
	refs := RefLoad(refPath)
	seqMap := map[string][]float64{
		"ACGTACGTACGTACGTACGTAC": {5.0}, // matches cleanly at both ends
		"ACGTACGTACGTACGTACGTNN": {2.0}, // contains N: excluded from alignment
	}
	got := AlignReads(seqMap, refs, 22)
	positions := got["amb"]["ACGTACGTACGTACGTACGTAC"]
	for _, p := range positions {
		abs := p
		if abs < 0 {
			abs = -abs
		}
		// no alignment window may overlap the N run at 1-based 23-26
		if abs <= 26 && abs+22-1 >= 23 {
			t.Errorf("alignment at %d overlaps the N run", p)
		}
	}
	if len(positions) == 0 {
		t.Error("expected alignments outside the N run")
	}
	if _, ok := got["amb"]["ACGTACGTACGTACGTACGTNN"]; ok {
		t.Error("read containing N must not align")
	}
}

// 32 nt reads exercise the packed-uint64 mask edge case.
func TestAlign_maxLengthReads(t *testing.T) {
	seq := "GGCATTCAGGATCGATCGGATCGATTAGGCAT" // 32 nt
	ref := ">m\nTTTT" + seq + "TTTT\n"
	refPath := writeTempFile(t, "max.fa", ref)
	refs := RefLoad(refPath)
	seqMap := map[string][]float64{seq: {1.0}}
	got := AlignReads(seqMap, refs, 32)
	want := map[string]map[string][]int{"m": {seq: {5}}}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}
}

// Gzipped read files identified by the .gz suffix must load identically.
func TestColumnNames_deduplication(t *testing.T) {
	names := columnNames([]string{"a/reads.fq", "b/reads.fq", "other.fq"})
	want := []string{"a/reads.fq", "b/reads.fq", "other.fq"}
	if !reflect.DeepEqual(names, want) {
		t.Errorf("got %v, want %v", names, want)
	}
}
