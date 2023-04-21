import pytest
import csv
import numpy as np
from src.scram2Plot import *

def test_init():
    dna = DNA("AGCT")
    assert dna.sequence == "AGCT"

def test_init_uppercase():
    dna = DNA("agct")
    assert dna.sequence == "AGCT"

def test_len():
    dna = DNA("AGCT")
    assert len(dna) == 4

def test_getitem():
    dna = DNA("AGCT")
    assert dna[0] == "A"
    assert dna[1] == "G"
    assert dna[2] == "C"
    assert dna[3] == "T"

def test_hash():
    dna1 = DNA("AGCT")
    dna2 = DNA("AGCT")
    assert hash(dna1) == hash(dna2)

def test_repr():
    dna = DNA("AGCT")
    assert repr(dna) == "AGCT"

def test_eq():
    dna1 = DNA("AGCT")
    dna2 = DNA("AGCT")
    dna3 = DNA("AGCA")
    assert dna1 == dna2
    assert dna1 != dna3

def test_single_alignment_init():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s.srna == "ACGT"
    assert s.position == 10
    assert s.strand == "+"
    assert s.times_aligned == 5
    assert np.array_equal(s.indv_alignments, [1, 2, 3, 4, 5])

def test_single_alignment_srna_len():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s.srna_len() == 4

def test_single_alignment_mean_alignments():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s.mean_alignments() == 3.0

def test_single_alignment_str():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert str(s) == "ACGT\t10\t+\t5\t[1, 2, 3, 4, 5]"

def test_single_alignment_eq():
    s1 = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    s2 = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    s3 = SingleAlignment("ACGT", 10, "-", 5, [1, 2, 3, 4, 5])
    s4 = SingleAlignment("ACGT", 10, "+", 6, [1, 2, 3, 4, 5, 6])
    s5 = SingleAlignment("TGCA", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s1 == s2
    assert s1 != s3
    assert s1 != s4
    assert s1 != s5

def test_single_ref_profile_init():
    s = SingleRefProfile()
    assert s.ref_len == 0
    assert s.all_alignments == []

def test_single_ref_profile_str():
    s = SingleRefProfile()
    assert str(s) == "0\t[]"

def test_single_ref_profile_eq():
    s1 = SingleRefProfile()
    s2 = SingleRefProfile()
    s3 = SingleRefProfile()
    s3.ref_len = 10
    s4 = SingleRefProfile()
    s4.all_alignments = [1, 2, 3]
    assert s1 == s2
    assert s1 != s3
    assert s1 != s4

def test_ref_profiles_init():
    r = RefProfiles()
    assert r.srna_len == 0
    assert r.replicates == 0
    assert r.single_ref_profiles == {}

def test_ref_profiles_str():
    r = RefProfiles()
    assert str(r) == "0\t0\t{}"

def test_ref_profiles_eq():
    r1 = RefProfiles()
    r2 = RefProfiles()
    r3 = RefProfiles()
    r3.srna_len = 20
    r4 = RefProfiles()
    r4.replicates = 5
    r5 = RefProfiles()
    r5.single_ref_profiles["ref1"] = SingleRefProfile()
    r5.single_ref_profiles["ref2"] = SingleRefProfile()
    assert r1 == r2
    assert r1 != r3
    assert r1 != r4
    assert r1 != r5

def test_ref_profiles_load_single_ref_profiles(tmp_path):
    r = RefProfiles()
    file_contents = "Header,ref_len,srna,position,strand,times_aligned,1,2,3\nr1,20,ACGT,10,+,5,1.0,2.0,3.0\nr1,20,ACGT,20,-,5,4.0,5.0,6.0\nr2,25,CGTA,5,+,5,7.0,8.0,9.0\n"
    file_path = tmp_path / "test_ref_profiles.csv"
    with open(file_path, "w") as f:
        f.write(file_contents)
    r.load_single_ref_profiles(file_path)
    assert r.srna_len == 4
    assert r.replicates == 3
    assert "r1" in r.single_ref_profiles
    assert "r2" in r.single_ref_profiles
    assert r.single_ref_profiles["r1"].ref_len == 20
    assert len(r.single_ref_profiles["r1"].all_alignments) == 2
    assert r.single_ref_profiles["r1"].all_alignments[0] == SingleAlignment(DNA("ACGT"), 10, "+", 5, [1.0, 2.0, 3.0])
    assert r.single_ref_profiles["r1"].all_alignments[1] == SingleAlignment(DNA("ACGT"), 20, "-", 5, [4.0, 5.0, 6.0])
    assert r.single_ref_profiles["r2"].ref_len == 25
    assert len(r.single_ref_profiles["r2"].all_alignments) == 1
    assert r.single_ref_profiles["r2"].all_alignments[0] == SingleAlignment(DNA("CGTA"), 5, "+", 5, [7.0, 8.0, 9.0])

#TODO: tests for DataForPlot