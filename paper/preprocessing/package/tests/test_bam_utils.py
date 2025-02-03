import atactoolbox.atac_repository as repo
import atactoolbox.bam_utils as bam_utils
import os
import pytest
import pysam

def test_add_CB_tag():

    bam_path = os.path.join(os.path.dirname(__file__),"data", "bam", "cropped_146.bam")
    cb_tagged = bam_utils.add_CB_tag(bam_path)
    #check if bamfile has CB tags
    read = next(pysam.AlignmentFile(cb_tagged))

    assert  read.has_tag("CB")

def test_is_bam():

    bam_path = os.path.join(os.path.dirname(__file__),"data", "bam", "cropped_146.bam")

    assert bam_utils.is_bam(bam_path)

def test_call_peak():

    tmp_dir = os.path.join(os.path.dirname(__file__),"data", "repo")
    tmp_dir = repo.make_tmp(tmp_dir)

    bam_path = os.path.join(os.path.dirname(__file__),"data", "bam")

    paths = repo.Repository(tmp_dir)
    paths.setSample("sample")
    paths.setClusterBamDir(bam_path)

    bam_utils.call_peak(paths, ['cropped_146.bam'], print)

    peakcalling = paths.getClusterPeakDir()
    files = os.listdir(peakcalling)

    assert len(files) == 5

    repo.rm_tmp(tmp_dir)
