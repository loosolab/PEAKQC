"""Tools for counting insertsizes from fragments or bam files."""
# Author: Jan Detleffsen (jan.detleffsen@mpi-bn.mpg.de)

import pandas as pd
import numpy as np
import datetime
from multiprocessing import Manager, Lock, Pool
from tqdm import tqdm

import peakqc.general as utils
import os
import re

from beartype import beartype
import numpy.typing as npt
from beartype.typing import Any, Optional


@beartype
def _check_in_list(element: Any, alist: list[Any] | set[Any]) -> bool:
    """
    Check if element is in list.

    Parameters
    ----------
    element : Any
        Element that is checked for.
    alist : list[Any] | set[Any]
        List or set in which the element is searched for.

    Returns
    -------
    bool
        True if element is in list else False
    """

    return element in alist


@beartype
def _check_true(element: Any, alist: Optional[list[Any]] = None) -> bool:  # true regardless of input
    """
    Return True regardless of input.

    Parameters
    ----------
    element : Any
        Element that is checked for.
    alist: Optional[list[Any]]
        List or set in which the element is searched for.

    Returns
    -------
    bool
        True if element is in list else False
    """

    return True


@beartype
def _custom_callback(error: Exception) -> None:
    """
    Error callback function for multiprocessing.

    Parameters
    ----------
    error : Exception
        Error that is raised.

    Returns
    -------
    None
    """
    print(error, flush=True)


GLOBAL_SHARD_LOCKS = None


def init_worker(shard_locks: list[Lock]) -> None:
    """
    Initialize global locks for worker processes.

    Parameters
    ----------
    shard_locks : list[Lock]
        List of locks for each shard.

    Returns
    -------
    None
    """
    global GLOBAL_SHARD_LOCKS
    GLOBAL_SHARD_LOCKS = shard_locks


@beartype
def insertsize_from_fragments(fragments: str,
                              barcodes: Optional[list[str]] = None,
                              chunk_size: int = 5000000,
                              n_threads: int = 8,
                              n_shards: int = 4,
                              verbose: bool = True) -> pd.DataFrame:
    """
    Count the insertsizes of fragments in a fragments file to obtain the insertsize size distribution, beside basic statistics (mean and total count) per barcode.

    Parameters
    ----------
    fragments : str
        Path to fragments file.
    barcodes : list[str], optional
        List of barcodes to count. If None, all barcodes are counted.
    chunk_size : int, default 5000000
        Size of chunks to split the genome into.
    n_threads : int, default 8
        Number of threads to use for multiprocessing.
    n_shards : int, default 4
        Number of shards to use for multiprocessing.
    verbose : bool, default True
        Print process shard affiliation.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the mean insertsizes and total counts per barcode.
    """
    print('Count insertsizes from fragments...')

    # Prepare function for checking against barcodes list
    if barcodes is not None:
        barcodes = set(barcodes)
        check_in = _check_in_list
    else:
        check_in = _check_true

    # Initialize iterator
    iterator = pd.read_csv(fragments,
                           delimiter='\t',
                           header=None,
                           names=['chr', 'start', 'stop', 'barcode', 'count'],
                           iterator=True,
                           chunksize=chunk_size)

    # start timer
    start_time = datetime.datetime.now()

    # Initialize multiprocessing
    # Initialize Manager and sharded dictionaries
    m = Manager()
    managed_dicts = [m.dict() for _ in range(n_shards)]
    # Initialize each managed dictionary as empty
    for i in range(n_shards):
        managed_dicts[i].update({})
    # Create one lock per shard
    shard_locks = [Lock() for _ in range(n_shards)]

    pool = Pool(processes=n_threads, initializer=init_worker, initargs=(shard_locks,))
    tasks = []
    max_active_tasks = n_shards + 1  # Limit active tasks to this number
    shard_counter = 0
    # Distribute chunks to workers, assigning shard index in round-robin fashion
    for i, chunk in enumerate(iterator):
        chunk = clean_chunk(chunk)
        shard_index = shard_counter % n_shards
        if verbose:
            print(f"Processing chunk {i} on shard {shard_index}")

        tasks.append(pool.apply_async(_count_fragments_worker,
                                      args=(chunk, shard_index, managed_dicts, barcodes, check_in),
                                      error_callback=_custom_callback))
        shard_counter += 1

        # Check active tasks in an unordered fashion
        while len(tasks) >= max_active_tasks:
            # Poll all tasks to see if any are ready
            ready_tasks = [t for t in tasks if t.ready()]
            if ready_tasks:
                for t in ready_tasks:
                    t.get()  # Retrieve result and catch exceptions if any
                    tasks.remove(t)
            else:
                # If none are ready, wait briefly on one task (using get with timeout or simply on the first one)
                tasks[0].get()  # This will wait for the first task to complete
                tasks.pop(0)

    # Wait for all worker tasks to complete
    for t in tasks:
        t.wait()

    pool.close()
    pool.join()

    if verbose:
        print('Merge results...')

    # Final merge: combine all shard dictionaries into one final result
    final_dict = {}
    for d in managed_dicts:
        final_dict = _update_count_dict(final_dict, dict(d))

    end_time = datetime.datetime.now()
    elapsed = end_time - start_time
    print("Final merging complete. Elapsed time:", str(elapsed).split(".")[0])

    # Convert dict to pandas dataframe
    print("Converting counts to dataframe...")
    table = pd.DataFrame.from_dict(final_dict, orient="index")
    # round mean_insertsize to 2 decimals
    table["mean_insertsize"] = table["mean_insertsize"].round(2)

    print("Done getting insertsizes from fragments!")

    return table


def clean_chunk(chunk: pd.DataFrame) -> pd.DataFrame:
    """
    Clean a chunk of a fragments file by removing rows with missing or malformed data.

    Removes rows where 'start' or 'stop' columns cannot be converted to numeric values (header or malformed rows).

    Parameters
    ----------
    chunk: pd.DataFrame

    Returns
    -------
    pd.DataFrame
        Cleaned chunk
    """
    # Attempt to convert the 'start' and 'stop' columns to numeric values
    chunk[['start', 'stop']] = chunk[['start', 'stop']].apply(pd.to_numeric, errors='coerce')
    # Drop rows where conversion failed (i.e., header or malformed rows)
    clean = chunk.dropna(subset=['start', 'stop'])
    # Convert columns back to int if necessary
    clean.loc[:, 'start'] = clean['start'].astype(int)
    clean.loc[:, 'stop'] = clean['stop'].astype(int)

    return clean


def _count_fragments_worker(chunk: pd.DataFrame,
                            shard_index: int,
                            managed_dicts: list[dict],
                            barcodes: Optional[list[str]] = None,
                            check_in: Any = _check_true) -> None:
    """
    Worker function to count fragments from a fragments file.

    Parameters
    ----------
    chunk : pd.DataFrame
        Chunk of the fragments file.
    shard_index : int
        Index of the shard to update.
    managed_dicts : list[dict], default None
        Dictionary for multiprocessing.
    barcodes : list[str], optional
        List of barcodes to count. If None, all barcodes are counted.
    check_in : Any, default _check_true
        Function for checking if barcode is in barcodes list.

    Returns
    -------
    None

    """

    # Initialize count_dict
    local_update = {}
    # Iterate over chunk
    for row in chunk.itertuples():
        start = int(row[2])
        end = int(row[3])
        barcode = row[4]
        count = int(row[5])
        size = end - start - 9  # length of insertion (-9 due to to shifted cutting of Tn5)

        # Only add fragment if check is true
        if check_in(barcode, barcodes) is True:
            local_update = _add_fragment(local_update, barcode, size, count)  # add fragment to count_dict

    # Update managed_dict
    with GLOBAL_SHARD_LOCKS[shard_index]:
        current_dict = dict(managed_dicts[shard_index])
        new_dict = _update_count_dict(current_dict, local_update)
        managed_dicts[shard_index].clear()
        managed_dicts[shard_index].update(new_dict)


@beartype
def _add_fragment(count_dict: dict,
                  barcode: str,
                  size: int,
                  count: int = 1,
                  max_size: int = 1000) -> dict:
    """
    Add fragment of size 'size' to count_dict.

    Parameters
    ----------
    count_dict : dict[str, int]
        Dictionary containing the counts per insertsize.
    barcode : str
        Barcode of the read.
    size : int
        Insertsize to add to count_dict.
    count : int, default 1
        Number of reads to add to count_dict.
    max_size : int, default 1000
        Maximum insertsize to consider.

    Returns
    -------
    dict
        Updated count_dict.
    """

    # Add read to dict
    if size > 0 and size <= max_size:  # do not save negative insertsize, and set a cap on the maximum insertsize to limit outlier effects a

        # Initialize if barcode is seen for the first time
        if barcode not in count_dict:
            count_dict[barcode] = {"mean_insertsize": 0, "insertsize_count": 0}

        count_dict[barcode]["insertsize_count"] += count

        # Update mean
        mu = count_dict[barcode]["mean_insertsize"]
        total_count = count_dict[barcode]["insertsize_count"]
        diff = (size - mu) / total_count
        count_dict[barcode]["mean_insertsize"] = mu + diff

        # Save to distribution
        if 'dist' not in count_dict[barcode]:  # initialize distribution
            count_dict[barcode]['dist'] = np.zeros(max_size + 1)
        count_dict[barcode]['dist'][size] += count  # add count to distribution

    return count_dict


@beartype
def _update_count_dict(dict1: dict, dict2: dict) -> dict:
    """
    Update the managed dict with the new counts.

    Parameters
    ----------
    dict1 : dict
        Dictionary containing the counts per insertsize.
    dict2 : dict
        Dictionary containing the counts per insertsize.

    Returns
    -------
    dict
        Updated count_dict.
    """
    merged_dict = {}
    # Gather all barcodes from both dictionaries
    all_barcodes = set(dict1.keys()) | set(dict2.keys())
    for barcode in all_barcodes:
        if barcode in dict1 and barcode in dict2:
            count1 = dict1[barcode]['insertsize_count']
            count2 = dict2[barcode]['insertsize_count']
            total_count = count1 + count2
            mean1 = dict1[barcode]['mean_insertsize']
            mean2 = dict2[barcode]['mean_insertsize']
            # Compute weighted mean
            weighted_mean = (mean1 * count1 + mean2 * count2) / total_count if total_count != 0 else 0

            # Merge distributions using element-wise addition (assumes numpy arrays)
            dist1 = dict1[barcode]['dist']
            dist2 = dict2[barcode]['dist']
            merged_dist = dist1 + dist2

            merged_dict[barcode] = {
                'insertsize_count': total_count,
                'mean_insertsize': weighted_mean,
                'dist': merged_dist
            }
        elif barcode in dict1:
            merged_dict[barcode] = dict1[barcode]
        else:
            merged_dict[barcode] = dict2[barcode]
    return merged_dict


@beartype
def _update_dist(dist_1: npt.ArrayLike, dist_2: npt.ArrayLike) -> npt.ArrayLike:
    """
    Update the Insertsize Distributions.

    Parameters
    ----------
    dist_1 : npt.ArrayLike
        Insertsize distribution 1.
    dist_2 : npt.ArrayLike
        Insertsize distribution 2.

    Returns
    -------
    npt.ArrayLike
        Updated insertsize distribution.
    """
    # check if both distributions are not empty
    if not np.isnan(dist_1).any() and not np.isnan(dist_2).any():
        updated_dist = dist_1 + dist_2  # add distributions
        return updated_dist.astype(int)
    # if one of the distributions is empty, return the other one
    elif np.isnan(dist_1).any():
        return dist_2.astype(int)
    elif np.isnan(dist_2).any():
        return dist_1.astype(int)


def insertsize_from_bam(bamfile: str,
                        barcodes: Optional[list[str]] = None,
                        barcode_tag: Optional[str] = 'CB',
                        chunk_size: int = 100000,
                        regions: Optional[list[str]] = None) -> pd.DataFrame:
    """
    Count the insertsizes of fragments in a bam file to obtain the insertsize size distribution, beside basic statistics (mean and total count) per barcode.

    Parameters
    ----------
    bamfile : str
        Path to bam file.
    barcodes : Optional[list[str]], default None
        List of barcodes to count. If None, all barcodes are counted.
    barcode_tag : Optional[str], default 'CB'
        Tag in bam file that contains the barcode.
    chunk_size : int, default 100000
        Size of chunks to split the genome into.
    regions : Optional[list[str]], default None
        List of regions to count. If None, the whole genome is counted.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the mean insertsizes and total counts per barcode.
    """

    # Check if pysam is installed
    utils.check_module("pysam")
    import pysam

    # Check if regions is a list
    if isinstance(regions, str):
        regions = [regions]

    # Prepare function for checking against barcodes list
    if barcodes is not None:
        barcodes = set(barcodes)  # convert to set for faster lookup
        check_in = _check_in_list
    else:
        check_in = _check_true

    # Open bamfile
    print("Opening bam file...")
    if not os.path.exists(bamfile + ".bai"):
        print("Bamfile has no index - trying to index with pysam...")
        pysam.index(bamfile)

    bam_obj = utils.open_bam(bamfile, "rb", require_index=True)
    # Get chromosome lengths
    chromosome_lengths = dict(zip(bam_obj.references, bam_obj.lengths))

    # Create chunked genome regions:
    print(f"Creating chunks of size {chunk_size}bp...")

    if regions is None:
        regions = [f"{chrom}:0-{length}" for chrom, length in chromosome_lengths.items()]
    elif isinstance(regions, str):
        regions = [regions]

    # Create chunks from larger regions
    regions_split = []
    for region in regions:
        chromosome, start, end = re.split("[:-]", region)
        start = int(start)
        end = int(end)
        for chunk_start in range(start, end, chunk_size):
            chunk_end = chunk_start + chunk_size
            if chunk_end > end:
                chunk_end = end
            regions_split.append(f"{chromosome}:{chunk_start}-{chunk_end}")

    # start timer
    start_time = datetime.datetime.now()

    # Count insertsize per chunk using multiprocessing
    print(f"Counting insertsizes across {len(regions_split)} chunks...")
    count_dict = {}  # initialize count_dict
    read_count = 0  # initialize read count
    # Iterate over chunks
    for region in tqdm(regions_split):
        chrom, start, end = re.split("[:-]", region)
        for read in bam_obj.fetch(chrom, int(start), int(end)):
            read_count += 1
            try:
                barcode = read.get_tag(barcode_tag)
            except Exception:  # tag was not found
                barcode = "NA"

            # Add read to dict
            if check_in(barcode, barcodes) is True:
                size = abs(read.template_length) - 9  # length of insertion
                count_dict = _add_fragment(count_dict, barcode, size)  # add fragment to count_dict

    # convert count_dict type float to int
    for barcode in count_dict:
        if 'dist' in count_dict[barcode]:
            count_dict[barcode]['dist'] = count_dict[barcode]['dist'].astype(int)

    # Close file and print elapsed time
    end_time = datetime.datetime.now()
    bam_obj.close()
    elapsed = end_time - start_time
    print("Done reading file - elapsed time: {0}".format(str(elapsed).split(".")[0]))

    # Convert dict to pandas dataframe
    print("Converting counts to dataframe...")
    table = pd.DataFrame.from_dict(count_dict, orient="index")
    # round mean_insertsize to 2 decimals
    table["mean_insertsize"] = table["mean_insertsize"].round(2)

    table = table[~pd.isna(table['dist'])]  # remove nan rows

    print("Done getting insertsizes from fragments!")

    return table


if __name__ == "__main__":
    print(os.getcwd())
    fragments = '/mnt/workspace2/jdetlef/experimental/16-peakqc/50m.tsv.gz'
    dist = insertsize_from_fragments(fragments, n_threads=8)
