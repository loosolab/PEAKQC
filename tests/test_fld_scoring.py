"""Test functions of the fld_scoring module."""
# Author: Jan Detleffsen (jan.detleffsen@mpi-bn.mpg.de)

import pytest
import os
import scanpy as sc
import numpy as np
import peakqc.fld_scoring as fld
import peakqc.insertsizes as ins

# ------------------------------- Fixtures and data -------------------------------- #


@pytest.fixture
def bamfile():
    """Return path to bam file."""

    return os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'heart_left_ventricle_1mio.bam')


@pytest.fixture
def fragments():
    """Return path to fragment file."""

    return os.path.join(os.path.dirname(__file__), 'data', 'insertsizes_related', 'fragments_heart_left_ventricle_head_100k.bed')


@pytest.fixture
def count_table(fragments):
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
    curves[1] = curves[1] / divs[1]  # Peak 2
    curves[2] = curves[2] / divs[2]  # Bias
    sum_c = np.sum(curves, axis=0)  # Sum of the curves

    return (sum_c * 100)


@pytest.fixture
def cosine_modulation():
    """Create a cosine based modulation curve."""

    def cosine(wavelength, amplitude=1, phase_shift=0):
        # Scale the wavelength and sigma with the scale
        wavl_scale = int(wavelength * 1.5)
        frequency = 1.5 / wavl_scale  # This ensures the frequency is scaled with scale

        x = np.arange(0, 1000)

        return amplitude * np.cos(2 * np.pi * frequency * x + phase_shift)

    cosine_modulation = cosine(25, 0.1) + cosine(150, 0.5) + cosine(75, 0.2) + cosine(100, 0.4)

    return cosine_modulation


@pytest.fixture
def fragment_distributions():
    """Load nucleosomal test data."""
    testdata = np.loadtxt(os.path.join(os.path.dirname(__file__), 'data', 'fld_scoring_related', 'nucleosomal_score.csv'), delimiter=None)

    return testdata


@pytest.fixture
def adata():
    """Fixture for an AnnData object."""
    adata = sc.read_h5ad(os.path.join(os.path.dirname(__file__), 'data', 'fld_scoring_related', 'subsetted_adata.h5ad'))
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

    assert gauss[350] < 0.1
    assert gauss[499] == 1
    assert gauss[650] < 0.1

    x = np.linspace(-4, 4, 1000)
    mu = 2
    sig = 0.4
    gauss = fld.gauss(x, mu, sig)

    assert gauss[600] < 0.1
    assert gauss[749] == 1
    assert gauss[900] < 0.1
    x = np.linspace(-4, 4, 1000)
    mu = 0
    sig = 1
    gauss = fld.gauss(x, mu, sig)

    assert gauss[200] < 0.1
    assert round(gauss[500], 1) == 1
    assert gauss[800] < 0.1


def test_build_score_mask():
    """Test that the build_score_mask function works as expected."""

    mask, figure = fld.build_score_mask(plot=True,
                                        save='density_plot_test.png',
                                        mu_list=[42, 200, 360, 550],
                                        sigma_list=[25, 35, 45, 25])

    assert (np.array([42, 200, 360, 549]) == np.concatenate(fld.call_peaks(mask))).all()

    ax_type = type(figure[0]).__name__
    assert ax_type.startswith("Figure")
    ax_type = type(figure[1]).__name__
    assert ax_type.startswith("Axes")

    assert os.path.isfile('density_plot_test.png')
    os.remove('density_plot_test.png')

    mask = fld.build_score_mask(plot=False,
                                save=None,
                                mu_list=[42, 200, 360, 550],
                                sigma_list=[25, 35, 45, 25])

    assert (np.array([42, 200, 360, 549]) == np.concatenate(fld.call_peaks(mask))).all()


def test_score_mask(good_modulation, bad_modulation):
    """Test that the score_mask function works as expected."""
    good_fit = fld.custom_conv(np.array([good_modulation]))
    good_peaks = fld.call_peaks(good_fit)
    good_peaks = fld.filter_peaks(good_peaks, reference=good_fit, peaks_thr=150, operator="bigger")
    good_score = fld.score_mask(good_peaks, good_fit)
    bad_fit = fld.custom_conv([bad_modulation])
    bad_peaks = fld.call_peaks(bad_fit)
    bad_peaks = fld.filter_peaks(bad_peaks, reference=bad_fit, peaks_thr=150, operator="bigger")
    bad_score = fld.score_mask(bad_peaks, bad_fit)

    assert good_score[0] > bad_score[0]


def test_custom_conv(good_modulation, bad_modulation):
    """Test that the custom_conv function works as expected."""
    good_fit = fld.custom_conv(good_modulation)
    bad_fit = fld.custom_conv(bad_modulation)

    assert np.sum(good_fit) > np.sum(bad_fit)


def test_cos_wavelet():
    """Test that the cos_wavelet function works as expected."""
    # check for the correct wavelength
    wav_1 = fld.cos_wavelet(wavelength=100,
                            amplitude=1.0,
                            phase_shift=0,
                            mu=0.0,
                            sigma=10,
                            plot=False,
                            save=None)

    peaks_1 = fld.call_peaks([wav_1])  # call peaks to get the center

    assert 99 == peaks_1[0][1] - peaks_1[0][0]  # check if the wavelength is correct
    assert 100 == peaks_1[0][2] - peaks_1[0][1]  # check if the wavelength is correct

    # check if the wavelet is centered
    wav_2 = fld.cos_wavelet(wavelength=100,
                            amplitude=1.0,
                            phase_shift=0,
                            mu=0.0,
                            sigma=0.4,
                            plot=False,
                            save=None)

    peaks_2 = fld.call_peaks([wav_2])  # call peaks to get the center
    assert np.where(np.max(wav_2) == wav_2)[0][0] == 149  # check if centered

    # check if sigma scales the cosine
    assert round(wav_1[peaks_1[0][0]] / wav_1[peaks_1[0][1]], 1) == 1
    assert round(wav_2[peaks_2[0][0]] / wav_2[peaks_2[0][1]], 1) == 0.3

    # check if its shifting
    wav_3 = fld.cos_wavelet(wavelength=100,
                            amplitude=1.0,
                            phase_shift=np.pi,
                            mu=0.0,
                            sigma=10,
                            plot=False,
                            save=None)

    peaks_3 = fld.call_peaks([wav_3])  # call peaks to get the center
    peaks_3[0][0] == 100
    peaks_3[0][1] == 199  # half the wavelength shift compared to before

    wav_4 = fld.cos_wavelet(wavelength=100,
                            amplitude=1.0,
                            phase_shift=0,
                            mu=100,
                            sigma=0.4,
                            plot=False,
                            save=None)

    assert np.where(np.max(wav_4) == wav_4)[0][0] == 249  # one wavelength shift for the gauss curve compared to before

    # test if plotting works
    _, figure = fld.cos_wavelet(wavelength=100,
                                amplitude=1.0,
                                phase_shift=0,
                                mu=0.0,
                                sigma=10,
                                plot=True,
                                save='cos_wavelet_test.png')

    ax_type = type(figure[0]).__name__
    assert ax_type.startswith("Figure")
    ax_type = type(figure[1]).__name__
    assert ax_type.startswith("Axes")

    # test if saving works
    assert os.path.isfile('cos_wavelet_test.png')
    os.remove('cos_wavelet_test.png')


def test_get_wavelets():
    """Test that the get_wavelet function works as expected."""
    wavelengths = [100, 200, 300]
    wavelets = fld.get_wavelets(wavelengths, sigma=0.4)

    assert len(wavelets) == 3


def test_wavelet_tranformation_fld(cosine_modulation):
    """Test that the wavelet_tranformation function works as expected."""
    wavelengths = [25, 50, 75, 100, 125, 150, 175, 200]
    wav_t = fld.wavelet_transform_fld(dists_arr=[cosine_modulation, cosine_modulation], wavelengths=wavelengths)

    assert round(np.mean(np.diff(fld.call_peaks([wav_t[0][0]], distance=1, width=0)))) == 25
    assert round(np.mean(np.diff(fld.call_peaks([wav_t[0][2]], distance=1, width=0)))) == 75
    assert round(np.mean(np.diff(fld.call_peaks([wav_t[0][4]], distance=1, width=0)))) == 100
    assert round(np.mean(np.diff(fld.call_peaks([wav_t[0][6]], distance=1, width=0)))) == 151


def test_score_by_conv(good_modulation, bad_modulation):
    """Test that the score_by_conv function works as expected."""
    good_score = fld.score_by_conv([good_modulation])
    bad_score = fld.score_by_conv([bad_modulation])

    assert good_score > bad_score


def test_plot_wavelet_transformation(cosine_modulation):
    """Test that the plot_wavelet_transformation function works as expected."""
    wavelengths = [25, 50, 75, 100, 125, 150, 175, 200]
    wav_t = fld.wavelet_transform_fld(dists_arr=[cosine_modulation], wavelengths=wavelengths)

    axes = fld.plot_wavelet_transformation(wav_t[0], wavelengths, cosine_modulation, save='wavelet_transformation_test.png')

    ax_type = type(axes[0]).__name__
    assert ax_type.startswith("Axes")

    assert os.path.isfile('wavelet_transformation_test.png')
    os.remove('wavelet_transformation_test.png')


def test_plot_custom_conv(good_modulation):
    """Test that the plot_custom_conv function works as expected."""
    convolution = fld.custom_conv([good_modulation])
    peaks = fld.call_peaks(convolution)
    scores = fld.score_mask(peaks, convolution)
    axes = fld.plot_custom_conv(convolution, [good_modulation], peaks, scores, save='custom_conv_test.png')

    ax_type = type(axes[0]).__name__
    assert ax_type.startswith("Axes")

    assert os.path.isfile('custom_conv_test.png')
    os.remove('custom_conv_test.png')


def test_add_fld_metrices(adata, fragments, bamfile):
    """Test that the add_fld_score function works as expected."""
    adata_f = adata.copy()
    adata_b = adata.copy()

    fld.add_fld_metrics(adata=adata_b,
                        fragments=bamfile,
                        barcode_col=None,
                        barcode_tag="CB",
                        chunk_size_bam=1000000,
                        regions=None,
                        peaks_thr=0.5,
                        wavelength=150,
                        sigma=0.4,
                        plot=False,
                        save_density=None,
                        save_overview=None,
                        sample=0)

    assert 'fld_score' in adata_b.obs.columns
    assert 'mean_fragment_size' in adata_b.obs.columns
    assert 'n_fragments' in adata_b.obs.columns

    fld.add_fld_metrics(adata=adata_f,
                        fragments=fragments,
                        barcode_col=None,
                        chunk_size_fragments=5000000,
                        peaks_thr=0.5,
                        wavelength=150,
                        sigma=0.4,
                        plot=False,
                        save_density=None,
                        save_overview=None,
                        sample=0,
                        n_threads=8)

    assert 'fld_score' in adata_f.obs.columns
    assert 'mean_fragment_size' in adata_f.obs.columns
    assert 'n_fragments' in adata_f.obs.columns
