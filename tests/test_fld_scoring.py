
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
    pre_table = ins.insertsize_from_fragments(fragments, barcodes=None)

    return np.array(pre_table['dist'].tolist(), dtype=np.int64)


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
def good_modulation():
    """Create a modulation curve."""

    mus = [45, 200, 360]
    sigs = [45, 55, 100]
    divs = [1, 2, 6]

    return modulation(mus, sigs, divs)


@pytest.fixture
def bad_modulation():
    """Create a modulation curve."""

    mus = [45, 100, 360]
    sigs = [45, 80, 100]
    divs = [1, 8, 10]

    return modulation(mus, sigs, divs)


def modulation(mus, sigs, divs):
    """Build a modulation curve."""
    def gaussian(x, mu, sig):  # Gaussian function
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    curves = []
    x_values = np.linspace(0, 1000, 1000)
    for mu, sig in [(mus[0], sigs[0]), (mus[1], sigs[1]), (mus[2], sigs[2])]:  # Gaussian curves with different means and standard deviations
        curves.append(gaussian(x_values, mu, sig))

    curves[1] = curves[1] / divs[0]  # Peak 1
    curves[1] = curves[1] / divs[1] # Peak 2
    curves[2] = curves[2] / divs[2]  # Bias
    sum_c = np.sum(curves, axis=0)  # Sum of the curves

    return (sum_c * 100)


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

    scaled = fld.scale(count_table)
    scaled_single = fld.scale(count_table[0])

    assert np.max(scaled) == 1
    assert np.min(scaled) == 0

    assert np.max(scaled_single) == 1
    assert np.min(scaled_single) == 0


def test_call_peaks_worker(good_modulation):
    """Test that the call_peaks_worker function works as expected."""
    peaks = fld.call_peaks_worker(good_modulation)

    assert (peaks == np.array([46, 204])).all()


def test_call_peaks(good_modulation):
    """Test that the call_peaks function works as expected."""
    peaks = fld.call_peaks([good_modulation, good_modulation])

    assert (peaks[0] == np.array([46, 204])).all()
    assert len(peaks) == 2


def test_filter_peaks(disturbed_sine):
    """Test that the filter_peaks function works as expected."""
    peaks = np.array([50, 250, 400, 500, 999])
    disturbed_sine_wave, sine_wave = disturbed_sine

    filtered_peaks = fld.filter_peaks(peaks, sine_wave, peaks_thr=0.75, operator="bigger")
    filtered_peaks_smaller = fld.filter_peaks(peaks, sine_wave, peaks_thr=0.75, operator="smaller")

    assert len(filtered_peaks) == 1
    assert filtered_peaks[0] == 250

    assert len(filtered_peaks_smaller) == 4
    assert np.all(filtered_peaks_smaller == np.array([50, 400, 500, 999]))


def test_density_plot(count_table):
    """Tests the density_plot function."""
    figure = fld.density_plot(count_table)

    ax = figure[0]
    ax_type = type(ax).__name__

    assert ax_type.startswith("Axes")


def test_gauss():
    """Test that the gauss function works as expected."""
    x = np.linspace(-4, 4, 1000)
    mu = 0
    sig = 0.5
    gauss = fld.gauss(x, mu, sig)

    assert gauss[350] < 0.05
    assert gauss[499] == np.max(gauss)
    assert gauss[650] < 0.05

    x = np.linspace(-4, 4, 1000)
    mu = 2
    sig = 0.4
    gauss = fld.gauss(x, mu, sig)

    assert gauss[600] < 0.05
    assert gauss[749] == np.max(gauss)
    assert gauss[900] < 0.05

    x = np.linspace(-4, 4, 1000)
    mu = 0
    sig = 1
    gauss = fld.gauss(x, mu, sig)

    assert gauss[200] < 0.05
    assert gauss[499] == np.max(gauss)
    assert gauss[800] < 0.05


def test_build_score_mask():
    """Test that the build_score_mask function works as expected."""

    mask, plot_related = fld.build_score_mask(plot=True,
                                              save=None,
                                              mu_list=[42, 200, 360, 550],
                                              sigma_list=[25, 35, 45, 25])

    assert (np.array([42, 200, 360, 549]) == np.concatenate(fld.call_peaks(mask))).all()

    ax_type = type(plot_related[0]).__name__
    assert ax_type.startswith("Figure")
    ax_type = type(plot_related[1]).__name__
    assert ax_type.startswith("Axes")

    mask = fld.build_score_mask(plot=False,
                                save=None,
                                mu_list=[42, 200, 360, 550],
                                sigma_list=[25, 35, 45, 25])

    assert (np.array([42, 200, 360, 549]) == np.concatenate(fld.call_peaks(mask))).all()


def test_score_mask(good_modulation, bad_modulation):
    """Test that the score_mask function works as expected."""
    pass


def test_custom_conv():
    """Test that the custom_conv function works as expected."""
    pass


def test_cos_wavelet():
    """Test that the cos_wavelet function works as expected."""
    pass

def test_get_wavelets():
    """Test that the get_wavelet function works as expected."""
    pass


