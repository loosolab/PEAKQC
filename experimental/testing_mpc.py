# individual imports
import pandas as pd
import numpy as np
import gzip
import datetime
from multiprocessing import Manager, Lock, Pool
from tqdm import tqdm

from beartype import beartype
from beartype.typing import Any, Optional

@beartype
def _is_gz_file(filepath: str) -> bool:
    """
    Check wheather file is a compressed .gz file.

    Parameters
    ----------
    filepath : str
        Path to file.

    Returns
    -------
    bool
        True if the file is a compressed .gz file.
    """

    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

class MPFragmentCounter():
    """
    """

    def __init__(self):
        """Init class variables."""
        pass

    def init_pool_processes(self, the_lock):
        '''
        Initialize each process with a global variable lock.
        '''
        global lock
        lock = the_lock

    def _check_in_list(self, element: Any, alist: list[Any] | set[Any]) -> bool:
        """
        Check if element is in list.

        TODO Do we need this function?

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

    def _check_true(element: Any, alist: Optional[list[Any]] = None) -> bool:  # true regardless of input
        """
        Return True regardless of input

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

    def insertsize_from_fragments(self, fragments: str,
                                  barcodes: Optional[list[str]] = None,
                                  n_threads: int = 8) -> pd.DataFrame:

        print('Count insertsizes from fragments...')
        # Open fragments file
        if _is_gz_file(fragments):
            f = gzip.open(fragments, "rt")
        else:
            f = open(fragments, "r")

        # Prepare function for checking against barcodes list
        if barcodes is not None:
            barcodes = set(barcodes)
            check_in = self._check_in_list
        else:
            check_in = self._check_true

        iterator = pd.read_csv(fragments,
                               delimiter='\t',
                               header=None,
                               names=['chr', 'start', 'stop', 'barcode', 'count'],
                               iterator=True,
                               chunksize=5000000)

        # start timer
        start_time = datetime.datetime.now()

        # Initialize multiprocessing
        m = Manager()
        lock = Lock()
        managed_dict = m.dict()
        managed_dict['output'] = {}
        pool = Pool(processes=n_threads, initializer=self.init_pool_processes, initargs=(lock,), maxtasksperchild=48)
        jobs = []
        print('Starting counting fragments...')
        # split fragments into chunks
        for chunk in tqdm(iterator, desc="Processing Chunks"):
            # apply async job wit callback function
            job = pool.apply_async(self._count_fragments_worker, args=(chunk, barcodes, check_in, managed_dict))
            jobs.append(job)
        # monitor progress
        # utils.monitor_jobs(jobs, description="Progress")
        # close pool
        pool.close()
        # wait for all jobs to finish
        pool.join()
        # reset settings
        count_dict = managed_dict['output']

        # Close file and print elapsed time
        end_time = datetime.datetime.now()
        f.close()

        elapsed = end_time - start_time
        print("Done reading file - elapsed time: {0}".format(str(elapsed).split(".")[0]))

        # Convert dict to pandas dataframe
        print("Converting counts to dataframe...")
        table = pd.DataFrame.from_dict(count_dict, orient="index")
        #table = table[["insertsize_count", "mean_insertsize"] + sorted(table.columns[2:])]
        table["mean_insertsize"] = table["mean_insertsize"].round(2)

        print("Done getting insertsizes from fragments!")

        return table

    def _count_fragments_worker(self, chunk, barcodes, check_in, managed_dict):
        """
        Worker function for counting fragments.
        Parameters
        ----------
        chunk
        barcodes
        check_in
        managed_dict

        Returns
        -------

        """

        # Initialize count_dict
        count_dict = {}
        for row in chunk.itertuples():
            start = int(row[2])
            end = int(row[3])
            barcode = row[4]
            count = int(row[5])
            size = end - start - 9  # length of insertion (-9 due to to shifted cutting of Tn5)

            # Only add fragment if check is true
            if check_in(barcode, barcodes) is True:
                count_dict = self._add_fragment(count_dict, barcode, size, count)

        with lock:
            latest = managed_dict['output']
            managed_dict['output'] = self._update_count_dict(latest, count_dict)

    def _add_fragment(self, count_dict: dict[str, int],
                      barcode: str,
                      size: int,
                      count: int = 1,
                      max_size=1000):
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
        """

        # Initialize if barcode is seen for the first time
        if barcode not in count_dict:
            count_dict[barcode] = {"mean_insertsize": 0, "insertsize_count": 0}

        # Add read to dict
        if size >= 0 and size <= max_size:  # do not save negative insertsize, and set a cap on the maximum insertsize to limit outlier effects

            count_dict[barcode]["insertsize_count"] += count

            # Update mean
            mu = count_dict[barcode]["mean_insertsize"]
            total_count = count_dict[barcode]["insertsize_count"]
            diff = (size - mu) / total_count
            count_dict[barcode]["mean_insertsize"] = mu + diff

            # Save to distribution
            if size not in count_dict[barcode]:  # first time size is seen
                count_dict[barcode]['dist'] = np.zeros(max_size)
            count_dict[barcode]['dist'][size] += count

        return count_dict

    def _update_count_dict(self, count_dict_1, count_dict_2):
        """
        updates
        """
        # Check if count_dict_1 is empty:
        if len(count_dict_1) == 0:
            return count_dict_2

        # make Dataframes for computation
        df1 = pd.DataFrame(count_dict_1).T
        df2 = pd.DataFrame(count_dict_2).T

        # merge distributions
        combined_dists = df1['dist'].combine(df2['dist'], func=self._update_dist)

        # merge counts
        merged_counts = pd.merge(df1["insertsize_count"], df2["insertsize_count"], left_index=True, right_index=True,
                                 how='outer').fillna(0)
        # sum total counts/barcode
        updated_counts = merged_counts.sum(axis=1)

        # calculate scaling factors
        x_scaling_factor = merged_counts["insertsize_count_x"] / updated_counts
        y_scaling_factor = merged_counts["insertsize_count_y"] / updated_counts

        # merge mean insertsizes
        merged_mean_insertsizes = pd.merge(df1["mean_insertsize"], df2["mean_insertsize"], left_index=True,
                                           right_index=True, how='outer').fillna(0)

        # scale mean insertsizes
        merged_mean_insertsizes["mean_insertsize_x"] = merged_mean_insertsizes["mean_insertsize_x"] * x_scaling_factor
        merged_mean_insertsizes["mean_insertsize_y"] = merged_mean_insertsizes["mean_insertsize_y"] * y_scaling_factor

        # sum the scaled means
        updated_means = merged_mean_insertsizes.sum(axis=1)

        # build the updated dictionary
        updated_dict = pd.DataFrame(
            {'mean_insertsize': updated_means, 'insertsize_count': updated_counts, 'dist': combined_dists}).T.to_dict()

        return updated_dict


    def _update_dist(self, dist_1, dist_2):
        """Updates the Insertsize Distributions"""
        if not np.isnan(dist_1).any() and not np.isnan(dist_2).any():
            updated_dist = dist_1 + dist_2
            return updated_dist
        elif np.isnan(dist_1).any():
            return dist_2
        elif np.isnan(dist_2).any():
            return dist_1


if __name__ == "__main__":
    # Test
    import episcanpy as epi

    #fragments = "/mnt/workspace2/jdetlef/data/public_data/fragments_heart_left_ventricle_194_sorted.bed"
    #h5ad_file = "/mnt/workspace2/jdetlef/data/public_data/heart_lv_SM-JF1NY.h5ad"#

    #fragments = "/home/jan/Workspace/bio_data/small_fragments.bed"
    fragments = "/home/jan/Workspace/bio_data/fragments_heart_left_ventricle_194_sorted.bed"
    h5ad_file = "/home/jan/Workspace/bio_data/heart_lv_SM-JF1NY.h5ad"
    adata = epi.read_h5ad(h5ad_file)
    adata_barcodes = adata.obs.index.tolist()
    # split index for barcodes CBs
    barcodes = []
    for entry in adata_barcodes:
        barcode = entry.split('+')[1]
        barcodes.append(barcode)

    counter = MPFragmentCounter()
    table = counter.insertsize_from_fragments(fragments, barcodes, n_threads=16)
    print(table)