# individual imports
import episcanpy as epi
import pandas as pd
import gzip
import datetime
import multiprocessing as mp
from collections import Counter

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

@beartype
class MPFragmentCounter():
    """
    """

    def __init__(self):
        """Init class variables."""

        self.merged_dict = None

    @beartype
    def _check_in_list(element: Any, alist: list[Any] | set[Any]) -> bool:
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

    @beartype
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
                               chunksize=1000)

        # start timer
        start_time = datetime.datetime.now()

        pool = mp.Pool(n_threads, maxtasksperchild=48)
        jobs = []
        # split fragments into chunks
        for chunk in iterator:
            # apply async job wit callback function
            job = pool.apply_async(self._count_fragments_worker, args=(chunk, barcodes, check_in),
                                   callback=self._log_result)
            jobs.append(job)
        # monitor progress
        # utils.monitor_jobs(jobs, description="Progress")
        # close pool
        pool.close()
        # wait for all jobs to finish
        pool.join()
        # reset settings
        count_dict = self.merged_dict.copy()
        self.merged_dict = None

        # Fill missing sizes with 0
        max_fragment_size = 1001

        for barcode in count_dict:
            for size in range(max_fragment_size):
                if size not in count_dict[barcode]:
                    count_dict[barcode][size] = 0

        # Close file and print elapsed time
        end_time = datetime.datetime.now()
        f.close()

        elapsed = end_time - start_time
        print("Done reading file - elapsed time: {0}".format(str(elapsed).split(".")[0]))

        # Convert dict to pandas dataframe
        print("Converting counts to dataframe...")
        table = pd.DataFrame.from_dict(count_dict, orient="index")
        table = table[["insertsize_count", "mean_insertsize"] + sorted(table.columns[2:])]
        table["mean_insertsize"] = table["mean_insertsize"].round(2)

        print("Done getting insertsizes from fragments!")

        return table

    def _count_fragments_worker(self, chunk, barcodes, check_in):

        for i in range(len(chunk)):
            row = chunk.iloc[i]
            start = int(row['start'])
            end = int(row['stop'])
            barcode = row['barcode']
            count = int(row['count'])
            size = end - start - 9  # length of insertion (-9 due to to shifted cutting of Tn5)

            # Only add fragment if check is true
            if check_in(barcode, barcodes) is True:
                count_dict = self._add_fragment(count_dict, barcode, size, count)

        return count_dict

    @beartype
    def _add_fragment(count_dict: dict[str, int],
                      barcode: str,
                      size: int,
                      count: int = 1) -> dict[str, int]:
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

        Returns
        -------
        dict[str, int]
            Updated count_dict
        """

        # Initialize if barcode is seen for the first time
        if barcode not in count_dict:
            count_dict[barcode] = {"mean_insertsize": 0, "insertsize_count": 0}

        # Add read to dict
        if size >= 0 and size <= 1000:  # do not save negative insertsize, and set a cap on the maximum insertsize to limit outlier effects

            count_dict[barcode]["insertsize_count"] += count

            # Update mean
            mu = count_dict[barcode]["mean_insertsize"]
            total_count = count_dict[barcode]["insertsize_count"]
            diff = (size - mu) / total_count
            count_dict[barcode]["mean_insertsize"] = mu + diff

            # Save to distribution
            if size not in count_dict[barcode]:  # first time size is seen
                count_dict[barcode][size] = 0
            count_dict[barcode][size] += count

        return count_dict

    def _log_result(self, result: Any) -> None:
        """Log results from mp_counter."""

        if self.merged_dict:
            self.merged_dict = dict(Counter(self.merged_dict) + Counter(result))
            # print('merging')
        else:
            self.merged_dict = result

if __name__ == "__main__":
    # Test
    import episcanpy as epi

    fragments = "/mnt/workspace2/jdetlef/data/public_data/fragments_heart_left_ventricle_194_sorted.bed"
    h5ad_file = "/mnt/workspace2/jdetlef/data/public_data/heart_lv_SM-JF1NY.h5ad"

    adata = epi.read_h5ad(h5ad_file)
    adata_barcodes = adata.obs.index.tolist()
    # split index for barcodes CBs
    barcodes = []
    for entry in adata_barcodes:
        barcode = entry.split('+')[1]
        barcodes.append(barcode)

    counter = MPFragmentCounter()
    table = counter.insertsize_from_fragments(fragments, barcodes, n_threads=8)
    print(table)