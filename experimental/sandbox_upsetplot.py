
# imports
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
import upsetplot
import matplotlib.pyplot as plt
# Imports
import sctoolbox
import sctoolbox.tools.qc_filter as qc
import sctoolbox.utils as utils
import sctoolbox.tools as tools
import sctoolbox.plotting as pl

import matplotlib.pyplot as plt
import episcanpy as epi
import pandas as pd
import scrublet as scr
import os
import numpy as np
import scanpy as sc

# type hint imports
from beartype.typing import Tuple, Dict, Optional, Literal, Callable, Iterable, Any  # , Union, List


def _upset_select_cells(adata: sc.AnnData,
                        thresholds: dict[str, dict[str, dict[Literal["min", "max"], int | float]] | dict[
                            Literal["min", "max"], int | float]],
                        groupby: Optional[str] = None,
                        direction: Literal['passed', 'filtered'] = 'passed') -> pd.DataFrame:
    """
    Select cells based on thresholds for UpSet Plot.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix object.
    thresholds : dict[str, dict[str, dict[Literal["min", "max"], int | float]] | dict[Literal["min", "max"], int | float]]
        Dictionary containing thresholds for each column. If groupby is given, thresholds are set per group.
    groupby : Optional[str], default None
        Name of the column in adata.obs to group cells by.

    Returns
    -------
    selection : pd.DataFrame
        DataFrame containing boolean values for each cell based on thresholds.
    """
    selection = {}
    # loop over all columns
    for column_name, values in thresholds.items():
        # loop over all groups
        if groupby:
            # initialize an array of False values
            accumulate_results = np.zeros(adata.obs.shape[0], dtype=bool)
            # loop over all samples
            for sample, cutoffs in values.items():
                # select cells based on the sample
                sample_selection = np.array(adata.obs[groupby] == sample)
                # select cells based on the cutoffs
                if direction == 'filtered':
                    cutoff_selection = np.array(
                        (adata.obs[column_name] < cutoffs['min']) | (adata.obs[column_name] > cutoffs['max']))
                elif direction == 'passed':
                    cutoff_selection = np.array(
                        (adata.obs[column_name] > cutoffs['min']) & (adata.obs[column_name] < cutoffs['max']))
                else:
                    # TODO add proper error handling
                    print('Error')

                # add up the selected cells
                accumulate_results = accumulate_results + np.logical_and(sample_selection, cutoff_selection)
            # add the results to the selection
            accumulate_series = pd.Series(accumulate_results, index=adata.obs.index)
            selection[column_name] = accumulate_series

        else:
            # select cells based on the cutoffs
            if direction == 'filtered':
                selection[column_name] = (adata.obs[column_name] < values['min']) | (
                            adata.obs[column_name] > values['max'])
            elif direction == 'passed':
                selection[column_name] = (adata.obs[column_name] > values['min']) & (
                            adata.obs[column_name] < values['max'])
            else:
                # TODO add proper error handling
                print('Error')

    # convert to DataFrame
    selection = pd.DataFrame(selection)

    return selection


def upset_plot_filter_impacts(adata: sc.AnnData,
                              thresholds: dict[str, dict[str, dict[Literal["min", "max"]]]],
                              limit_combinations: Optional[int] = None,
                              groupby: Optional[int] = None,
                              direction: Literal['passed', 'filtered'] = 'passed',
                              log: Optional[bool] = False) -> Optional[dict]:
    """
    Plot the impact of filtering cells based on thresholds in an UpSet Plot.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix object.
    thresholds : dict[str, dict[str, dict[Literal["min", "max"], int | float]] | dict[Literal["min", "max"], int | float]]
        Dictionary containing thresholds for each column. If groupby is given, thresholds are set per group.
    limit_combinations : Optional[int], default None
        Limit the number of combinations to show in the plot.
    groupby : Optional[str], default None
        Name of the column in adata.obs to group cells by.

    Returns
    -------
    plot_result : Optional[dict]
    """
    if len(thresholds) <= 1:
        logger.info("Skipping UpSet Plot as only one threshold is given.")

        return None

    selection = _upset_select_cells(adata, thresholds, groupby, direction)

    # Number of variables
    n = len(selection.columns)

    # Generate all combinations of True/False
    raw_combinations = np.array(np.meshgrid(*[[False, True]] * n)).T.reshape(-1, n)

    # Sort the combinations first by the number of True values, then lexicographically
    sorted_combinations = np.array(sorted(raw_combinations, key=lambda x: (np.sum(x), list(x))))

    # exclude empty combinations
    mask = np.sum(sorted_combinations, axis=1) != 0
    combinations = sorted_combinations[mask]

    # make a dataframe
    combinations_df = pd.DataFrame(combinations, columns=selection.columns)

    counts = []
    barcodes = []

    # loop over all combinations
    for _, combination in combinations_df.iterrows():
        # get a list of the column names selected
        metric_combination = list(selection.columns[combination])
        # initialize a results array for the combination with len n-cells
        single_result = np.zeros(selection.shape[0], dtype=bool)
        # loop over the selected combination
        for metric in metric_combination:
            # check if initialized
            if sum(single_result) == 0:
                single_result = selection[metric]
            else:
                # add up the selected cells
                if direction == 'filtered':
                    single_result = np.logical_or(single_result, selection[metric])
                elif direction == 'passed':
                    single_result = np.logical_and(single_result, selection[metric])

        # append result of the combination
        counts.append(sum(single_result))
        # append the barcodes of the cells
        barcodes.append(single_result.iloc[single_result.values].index)
    # add to dataframe
    combinations_df['counts'] = counts
    combinations_df['barcodes'] = barcodes

    # limit combinations
    if limit_combinations:
        # select all combinations with a grade less or equal the limit
        limit_mask = np.array(np.sum(combinations_df[selection.columns], axis=1) <= limit_combinations)
        # always include the total counts
        limit_mask[-1] = True
        # index by the mask
        combinations_df = combinations_df[limit_mask]

    # set the combinations as index
    combinations_df.set_index(list(selection.columns), inplace=True)

    if log:
        combinations_df['counts'] = np.log(combinations_df['counts'])

    with warnings.catch_warnings():  # TODO remove when this is merged https://github.com/jnothman/UpSetPlot/pull/278
        warnings.filterwarnings("ignore",
                                message="A value is trying to be set on a copy of a DataFrame or Series through chained assignment using an inplace method.")

        # plot the UpSet Plot
        plot_result = upsetplot.plot(combinations_df['counts'], totals_plot_elements=0)

    if direction == 'filtered':
        plot_result["intersections"].set_ylabel("Cells Filtered")
    elif direction == 'passed':
        plot_result["intersections"].set_ylabel("Cells Passed")
    plt.show()

    return plot_result, combinations_df

if __name__=='__main__':

    adata = utils.adata.load_h5ad('/mnt/workspace2/jdetlef/peakqc_paperprep/benchmarking/assembled/right_atrium_auricular_region_IOBHN.h5ad')

    adata = adata[adata.X.sum(axis=1) > 0]
    adata = adata[:, adata.X.sum(axis=0) > 0]

    adata = tools.qc_filter.calculate_qc_metrics(adata, var_type='features')

    # Decide whether to estimate thresholds individual per condition (False) or globally (True)
    global_threshold = False

    # Set the column in adata.obs containing the biological condition to evaluate
    condition_column = "sample"

    thresholds = {
                  'n_features': {'min': 150, 'max': 5000},
                  #'log1p_n_features': {'min': None, 'max': None},
                  'nucleosome_signal': {'min': 0.1, 'max': 0.7},
                  'fld_score': {'min': 100, 'max': 5000},
                  'frip': {'min': 0.1, 'max': 2}
                  #'tsse_score': {'min': 2, 'max': 1000}
                  # add additional threshold based on the available columns shown above
                  # format: '<obs clolumn>': {'min': <threshold|None>, 'max': <threshold|None>}
                 }

    obs_columns = list(thresholds.keys())

    groupby = condition_column if global_threshold is False else None
    initial_thresholds = tools.qc_filter.get_thresholds(adata,
                                                        thresholds,
                                                        groupby=groupby)
    tools.qc_filter.thresholds_as_table(thresholds)

    _, combinations_df = upset_plot_filter_impacts(adata,
                                                   thresholds=initial_thresholds,
                                                   groupby=groupby,
                                                   limit_combinations=3,
                                                   direction='filtered',
                                                   log=False)

    print(combinations_df.index)

    print(combinations_df)