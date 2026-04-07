"""Test functions of the insertsizes module."""
# Author: Jan Detleffsen (jan.detleffsen@mpi-bn.mpg.de)

import os

import numpy as np
import pandas as pd
import pytest
from multiprocessing import Lock
import peakqc.insertsizes as insertsizes


@pytest.fixture
def tsv_head_tail():
    """Return a tuple of the header section mixed with body and body only."""
    fragments = os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', '200_lines_atac_pbmc_5k_nextgem_fragments.tsv.gz')

    iterator = pd.read_csv(fragments,
                           delimiter='\t',
                           header=None,
                           names=['chr', 'start', 'stop', 'barcode', 'count'],
                           dtype={'barcode': str},
                           iterator=True,
                           chunksize=100)

    chunk_with_header = next(iterator)
    chunk_without_header = next(iterator)

    return chunk_with_header, chunk_without_header


@pytest.fixture
def barcodes():
    """Return a barcodes file."""
    path_barcodes = os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'barcodes.txt')

    with open(path_barcodes, 'r') as file:
        return [line.strip() for line in file]


@pytest.fixture
def fragments_file():
    """Return a fragments file."""

    return os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'fragments_heart_left_ventricle_head_100k.bed')


@pytest.fixture
def int_barcodes_fragments_file():
    """Return a fragments file with integer barcodes."""

    return os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'int_barcodes_fragments.bed')


@pytest.fixture
def bam_file():
    """Return a bam file."""

    return os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'heart_left_ventricle_1mio.bam')


@pytest.fixture
def fragments():
    """Return a fragments file."""

    return pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'fragments_heart_left_ventricle_head_100k.bed'))


@pytest.fixture
def chunk():
    """Return a chunk of a fragments file."""

    return pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'example_chunk.csv'))


def test_init_worker_sets_global_shard_locks():
    """Test that the init_worker function sets the GLOBAL_SHARD_LOCKS variable."""

    test_locks = [Lock() for _ in range(3)]
    insertsizes._init_worker(test_locks)
    assert insertsizes.GLOBAL_SHARD_LOCKS == test_locks


def test_check_in_list():
    """Check the _check_in_list function."""

    # Create a mock list of cell barcodes
    mock_list = ['ATTGCTAACCGGC', 'ACAAGGCTTGGCA', 'ACGTTGCTTGGCA', 'ACGTTGCTTGGCA']

    # Check that the function returns True for a barcode in the list
    assert insertsizes._check_in_list('ATTGCTAACCGGC', mock_list) is True


def test_check_true():
    """Check the _check_true function."""
    mock_list = ['ATTGCTAACCGGC', 'ACAAGGCTTGGCA', 'ACGTTGCTTGGCA', 'ACGTTGCTTGGCA']
    # check that the function returns True for any input
    assert insertsizes._check_true('NOT in the LIST', mock_list) is True


def test_custom_callback(capfd):
    """Test the custom_callback function."""

    # Create an exception to pass to the callback
    test_exception = Exception("Test exception")

    # Call the _custom_callback function with the test exception
    insertsizes._custom_callback(test_exception)

    # Capture the output
    out, err = capfd.readouterr()

    # Assert that the output contains the exception message
    assert "Test exception" in out


def test_count_fragments_worker(chunk):
    """Test the count_fragments_worker function."""

    # Call the function with the mock lock
    insertsizes._init_worker([Lock()])

    # Init a dictionary to store the results
    insertsizes_dicts = [{}]

    # Call the function with the mock lock and the mock dictionary
    insertsizes._count_fragments_worker(chunk, 0, insertsizes_dicts)

    # Check that the dictionary contains the expected results
    assert len(insertsizes_dicts[0]) == 18881  # Number of unique cell barcodes in the chunk
    assert round(insertsizes_dicts[0]['AGGGATAAACCACCGAAGGTCA']['mean_insertsize']) == 140  # Mean insert size for a specific cell barcode
    assert insertsizes_dicts[0]['AGGGATAAACCACCGAAGGTCA']['insertsize_count'] == 38  # Number of fragments for a specific cell barcode
    assert int(insertsizes_dicts[0]['AGGGATAAACCACCGAAGGTCA']['dist'].sum()) == 38  # Number of fragments in the distribution for a specific cell barcode


def test_add_fragment_counts():
    """Test the add_fragment_counts function."""

    # Create a mock dictionary of results
    count_dict = {'AGGGATAAACCACCGAAGGTCA': {'mean_insertsize': 140,
                                             'insertsize_count': 38,
                                             'dist': np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                                               0., 0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 1., 0., 0., 0.,
                                                               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])}}

    # Call the function with the mock dictionary
    count_dict = insertsizes._add_fragment(count_dict=count_dict,
                                           barcode='AGGGATAAACCACCGAAGGTCA',
                                           size=25,
                                           count=3,
                                           max_size=50)

    # Check that the dictionary contains the expected results
    assert count_dict['AGGGATAAACCACCGAAGGTCA']['insertsize_count'] == 41  # Number of fragments for a specific cell barcode
    assert int(count_dict['AGGGATAAACCACCGAAGGTCA']['dist'][25] == 3)  # Number of fragments in the distribution for a specific cell barcode


def test_update_count_dict():
    """Test the update_count_dict function."""

    # Create a mock dictionary with a single cell barcode
    count_dict_1 = {'AGGGATAAACCACCGAAGGTCA': {'mean_insertsize': 5.5,
                                               'insertsize_count': 2,
                                               'dist': np.array([0., 0., 0., 0., 1., 1., 0., 0., 0., 0.])}}

    # Create a mock dictionary with multiple cell barcodes
    count_dict_2 = {'AGGGATAAACCACCGAAGGTCA': {'mean_insertsize': 4.5,
                                               'insertsize_count': 2,
                                               'dist': np.array([0., 0., 0., 1., 1., 0., 0., 0., 0., 0.])},
                    'AGGGATAAACCACCGAAGGTCC': {'mean_insertsize': 10,
                                               'insertsize_count': 1,
                                               'dist': np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1.])},
                    }

    # Call the function with the mock dictionaries
    count_dict = insertsizes._update_count_dict(count_dict_1, count_dict_2)

    # Check that the dictionary contains the expected results
    assert len(count_dict) == 2  # Number of unique cell barcodes in the dictionary
    assert count_dict['AGGGATAAACCACCGAAGGTCA']['insertsize_count'] == 4  # Number of fragments for a specific cell barcode
    assert count_dict['AGGGATAAACCACCGAAGGTCC']['insertsize_count'] == 1  # Number of fragments for a specific cell barcode
    assert count_dict['AGGGATAAACCACCGAAGGTCA']['dist'][3] == 1  # Number of fragments in the distribution for a specific cell barcode
    assert count_dict['AGGGATAAACCACCGAAGGTCA']['dist'][4] == 2


def test_insertsize_from_fragments(fragments_file, barcodes):
    """Test the insertsize_from_fragments function."""

    table = insertsizes.insertsize_from_fragments(fragments=fragments_file,
                                                  barcodes=None,
                                                  n_threads=8)

    assert table.shape[0] == 11219  # Number of unique cell barcodes in the table
    assert round(table.loc['AGGGATAAACCACCGAAGGTCA', 'mean_insertsize']) == 130  # Mean insert size for a specific cell barcode
    assert table.loc['AGGGATAAACCACCGAAGGTCA', 'insertsize_count'] == 6  # Number of fragments for a specific cell barcode
    assert int(table.loc['AGGGATAAACCACCGAAGGTCA', 'dist'].sum()) == 6  # Number of fragments in the distribution for a specific cell barcode

    table = insertsizes.insertsize_from_fragments(fragments=fragments_file,
                                                  barcodes=barcodes,
                                                  n_threads=8)

    assert table.shape[0] == 7798  # Number of unique cell barcodes in the table


def test_insertsize_from_bam(bam_file, barcodes):
    """Test the insertsize_from_bam function."""
    table = insertsizes.insertsize_from_bam(bam_file,
                                            barcodes=None,
                                            barcode_tag='CB',
                                            chunk_size=100000,
                                            regions=None)

    assert table.shape[0] == 16182  # Number of unique cell barcodes in the table
    assert round(table.loc['AGGGATAAACCACCGAAGGTCA', 'mean_insertsize']) == 92  # Mean insert size for a specific cell barcode
    assert table.loc['AGGGATAAACCACCGAAGGTCA', 'insertsize_count'] == 280  # Number of fragments for a specific cell barcode
    assert int(table.loc['AGGGATAAACCACCGAAGGTCA', 'dist'].sum()) == 280  # Number of fragments in the distribution for a specific cell barcode

    table = insertsizes.insertsize_from_bam(bam_file,
                                            barcodes=barcodes,
                                            barcode_tag='CB',
                                            chunk_size=100000,
                                            regions=None)

    assert table.shape[0] == 9079  # Number of unique cell barcodes in the table

    table = insertsizes.insertsize_from_bam(bam_file,
                                            barcodes=barcodes,
                                            barcode_tag='CB',
                                            chunk_size=100000,
                                            regions='chr1:1-100000')

    assert table.shape[0] == 2705  # Number of unique cell barcodes in the table


def test_insertsize_from_fragments_int_barcodes(int_barcodes_fragments_file):
    """Test that insertsize_from_fragments works when barcodes are integers."""

    table = insertsizes.insertsize_from_fragments(fragments=int_barcodes_fragments_file,
                                                  barcodes=None,
                                                  n_threads=1)

    assert '19401152' in table.index
    assert table.loc['19401152', 'insertsize_count'] == 2

    table_filtered = insertsizes.insertsize_from_fragments(fragments=int_barcodes_fragments_file,
                                                           barcodes=['19401152', '56420551'],
                                                           n_threads=1)
    assert table_filtered.shape[0] == 2
    assert '14581843' not in table_filtered.index


@pytest.mark.parametrize('chunk, response', [(0, 49), (1, 100)])
def test_clean_chunk(tsv_head_tail, chunk, response):
    """Test the clean_chunk function."""
    cleaned = insertsizes._clean_chunk(tsv_head_tail[chunk])
    assert len(cleaned) == response
