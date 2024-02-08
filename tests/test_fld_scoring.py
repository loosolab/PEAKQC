
import pytest
import os
import scanpy as sc
import numpy as np
import peakqc.fld_scoring as fld
import peakqc.insertsizes as ins

# ------------------------------- Fixtures and data -------------------------------- #


# Get the paths to the test data; not as fixtures because they are used in parametrized tests
fragments = os.path.join(os.path.dirname(__file__), 'data', 'fld_scoring_related', 'mm10_atac_fragments.bed')
bamfile = os.path.join(os.path.dirname(__file__), 'data', 'fld_scoring_related', 'mm10_atac.bam')

@pytest.fixture
def count_table():
    """Return fragment count table."""
    return ins.insertsize_from_fragments(fragments, barcodes=None)


@pytest.fixture
def disturbed_sine(freq=3.1415 * 2):
    """Return list of disturbed sine wave and sine wave."""
    in_array = np.linspace(0, freq, 1000)
    sine_wave = np.sin(in_array)
    in_array = np.linspace(0, 500, 1000)
    disturbance = np.sin(in_array)
    scaled_disturbance = disturbance / 10
    disturbed_sine = sine_wave + scaled_disturbance

    return disturbed_sine, sine_wave


@pytest.fixture
def stack_sines(disturbed_sine):
    """Return multiple sine waves and disturbed sine waves."""

    sines = []
    disturbed_sine_waves = []
    for i in range(10):
        disturbed_sine_wave, sine_wave = disturbed_sine
        sines.append(sine_wave)
        disturbed_sine_waves.append(disturbed_sine_wave)

    sines = np.array(sines)
    disturbed_sine_waves = np.array(disturbed_sine_waves)

    return sines, disturbed_sine_waves


@pytest.fixture
def modulation():
    """Create a modulation curve."""
    def gaussian(x, mu, sig):  # Gaussian function
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    curves = []
    x_values = np.linspace(-3, 3, 1000)
    for mu, sig in [(-1, 0.5), (0.2, 0.3), (-1, 3)]:  # Gaussian curves with different means and standard deviations
        curves.append(gaussian(x_values, mu, sig))

    curves[1] = curves[1] / 1  # Peak 1
    curves[1] = curves[1] / 5  # Peak 2
    curves[2] = curves[2] / 5  # Bias
    sum_c = np.sum(curves, axis=0)  # Sum of the curves

    return sum_c


@pytest.fixture
def fragment_distributions():
    """Load nucleosomal test data."""
    testdata = np.loadtxt(os.path.join(os.path.dirname(__file__), 'data', 'fld_scoring_related', 'nucleosomal_score.csv'), delimiter=None)

    return testdata


@pytest.fixture
def adata():
    """Fixture for an AnnData object."""
    adata = sc.read_h5ad(os.path.join(os.path.dirname(__file__), 'data', 'fld_scoring_related', 'mm10_atac.h5ad'))
    return adata


# --------------------------------- Tests ---------------------------------------- #

def test_moving_average(disturbed_sine):
    """
    Test that the moving average function works as expected.

    Compares a smoothed disturbed sine wave to the original and inspects the difference.
    """
    disturbed_sine_wave, sine_wave = disturbed_sine
    smoothed_sine = fld.moving_average(disturbed_sine_wave, n=10)

    diff_smooth = np.sum(abs(sine_wave - smoothed_sine))
    diff_disturbed = np.sum(abs(sine_wave - disturbed_sine_wave))

    print("\t")
    print("smoothed difference: " + str(diff_smooth))
    print("disturbed difference: " + str(diff_disturbed))
    assert diff_smooth < 15


def test_multi_ma(stack_sines):
    """Test that the multi_ma function works as expected by comparing a smoothed disturbed sine wave to the original."""
    sine_stack, dist_stack = stack_sines
    smoothed = fld.multi_ma(dist_stack)

    diff_ori = abs(sine_stack - dist_stack)
    diff_smooth = abs(sine_stack - smoothed)

    sum_ori = np.sum(diff_ori, axis=1)
    sum_smooth = np.sum(diff_smooth, axis=1)

    print("\t")
    print("smoothed difference: " + str(sum_smooth))
    print("disturbed difference: " + str(sum_ori))

    assert np.all(sum_smooth < 15)


def test_scale(count_table):
    """Test that the scale function works as expected by checking that the max value is 1 and the min value is 0."""
    table = count_table

    dist = table[[c for c in table.columns if isinstance(c, int)]]
    dists_arr = dist.to_numpy()
    dists_arr = np.nan_to_num(dists_arr)

    scaled = fld.scale(dists_arr)
    scaled_single = fld.scale(dists_arr[0])

    assert np.max(scaled) == 1
    assert np.min(scaled) == 0

    assert np.max(scaled_single) == 1
    assert np.min(scaled_single) == 0


def test_get_dist_df(count_tabele):
