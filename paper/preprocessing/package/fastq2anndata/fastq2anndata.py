import os
import subprocess
import sys
import pandas as pd
from pathlib import Path
from snaptools.alignment import index_ref
from snaptools.alignment import run_align_pe
from snaptools.snap_pre import snap_pre
from snaptools.add_bmat import snap_bmat
from snaptools.add_pmat import snap_pmat
import fastq2anndata.converter as cv
import pysam
from sinto import fragments

def preprocessing(index_settings, align_settings, peakcalling_settings, snap_settings, converter_settings, general_settings):
    '''

    :param genome:
    :param genome_prefix:
    :param aligner:
    :param path_to_aligner:
    :param num_threads:
    :return:
    '''

    verbose = general_settings['verbose']
    tmp_folder = general_settings['tmp_folder']
    num_threads = general_settings['n_threads']

    general_settings = auto_path(general_settings)

    verbose_1 = general_settings['index_genome']
    print(f'index genome: {verbose_1}')

    if general_settings['index_genome']:
        index_ref(input_fasta=index_settings['genome'],
                  output_prefix=index_settings['genome_prefix'],
                  aligner=index_settings['aligner'],
                  path_to_aligner=index_settings['path_to_aligner'],
                  num_threads=num_threads)

    else:
        print('Index-Genome: SKIPPED')

    verbose_1 = general_settings['align']
    verbose_2 = align_settings['multi']
    print(f'align fastq: {verbose_1}, multi: {verbose_2}')
    if general_settings['align']:
        # align paired end
        if align_settings['multi']:

            samples_dict = align_settings['samples_dict']
            for sample, value in samples_dict.items():
                fastq1 = os.path.join(align_settings['fastq_dir'], value[0])
                fastq2 = os.path.join(align_settings['fastq_dir'], value[1])

                output_bam = general_settings['output_bam'] + '/' + sample + '.bam'

                print("Writing bam: " + output_bam)
                run_align_pe(input_reference=index_settings['genome'],
                             input_fastq1=fastq1,
                             input_fastq2=fastq2,
                             output_bam=output_bam,
                             aligner=index_settings['aligner'],
                             path_to_aligner=index_settings['path_to_aligner'],
                             num_threads=num_threads,
                             if_sort=align_settings['if_sort'],
                             min_cov=align_settings['min_cov'],
                             aligner_options=align_settings['args'],
                             read_fastq_command=align_settings['read_fastq_command'],
                             tmp_folder=tmp_folder,
                             overwrite=align_settings['overwrite'],
                             verbose=verbose)

        else:

            output_bam = general_settings['output_bam'] + '/' + align_settings['sample'] + '.bam'

            print("Writing bam: " + output_bam)
            run_align_pe(input_reference=index_settings['genome'],
                         input_fastq1=align_settings['input_fastq1'],
                         input_fastq2=align_settings['input_fastq2'],
                         output_bam=output_bam,
                         aligner=index_settings['aligner'],
                         path_to_aligner=index_settings['path_to_aligner'],
                         num_threads=num_threads,
                         if_sort=align_settings['if_sort'],
                         min_cov=align_settings['min_cov'],
                         aligner_options=align_settings['args'],
                         read_fastq_command=align_settings['read_fastq_command'],
                         tmp_folder=tmp_folder,
                         overwrite=align_settings['overwrite'],
                         verbose=verbose)

    else:
        print('Alignment: SKIPPED')

    verbose_1 = general_settings['process_alignment']
    print(f'process alignment: {verbose_1}')
    if general_settings['process_alignment']:
        process_alignment(general_settings, align_settings)

    else:
        print('bam preprocessing: SKIPPED')

    if general_settings['callpeak']:
        print('performing peakcalling')
        if peakcalling_settings['optional_input']:
            print('using optional input')
            call_peak(general_settings['output_peaks'], align_settings['sample'], peakcalling_settings['optional_bam'], print)
        else:
            if align_settings['multi']:
                print('mode: multi')
                samples_dict = align_settings['samples_dict']
                for sample in samples_dict.keys():

                    input_bam = general_settings['output_bam'] + '/' + sample + '.bam'
                    call_peak(general_settings['output_peaks'], sample, input_bam, print)

            else:
                print('mode: single')
                input_bam = general_settings['output_bam'] + '/' + align_settings['sample'] + '.bam'
                call_peak(general_settings['output_peaks'], align_settings['sample'], input_bam, print)

        # remove all columns other than chr, start, end from narrowPeak file and save as bed file
        # save bedfiles in a dict with sample name as key
        bedfiles = {}
        for file in os.listdir(general_settings['output_peaks']):
            if file.endswith('.narrowPeak'):
                output = general_settings['output_peaks'] + '/' + file.split('.narrowPeak')[0] + '.bed'
                #check if file already exists
                if os.path.isfile(output):
                    print('File already exists: ' + output)
                else:
                    print('crop file: ' + file)
                    df = pd.read_csv(general_settings['output_peaks'] + '/' + file, sep='\t', header=None)
                    # change first column to b'*' to make it compatible with snaptools
                    df.iloc[:, 0] = df.iloc[:, 0].apply(lambda x: b'' + x.encode())

                    df = df.iloc[:, 0:3]
                    df.to_csv(output, sep='\t', header=False, index=False) # save bed file
                    bedfiles[file.split('.')[0]] = output
    else:
        print('Peakcalling: SKIPPED')

    if general_settings['build_snap']:

        if align_settings['multi']:
            samples_dict = align_settings['samples_dict']
            for sample in samples_dict.keys():

                input_bam = general_settings['output_bam'] + '/' + sample + '.bam'
                output_snap = general_settings['output_snap'] + '/' + sample + '.snap'
                snap_pre(input_file=input_bam,
                         output_snap=output_snap,
                         genome_name=index_settings['genome_prefix'],
                         genome_size=snap_settings['genome_size'],
                         min_mapq=snap_settings['snap_min_mapq'],
                         min_flen=snap_settings['snap_min_flen'],
                         max_flen=snap_settings['snap_max_flen'],
                         min_cov=snap_settings['snap_min_cov'],
                         max_num=snap_settings['snap_max_num'],
                         barcode_file=snap_settings['snap_barcode_file'],
                         keep_chrm=snap_settings['snap_keep_chrm'],
                         keep_single=snap_settings['snap_keep_single'],
                         keep_secondary=snap_settings['snap_keep_secondary'],
                         keep_discordant=snap_settings['snap_keep_discordant'],
                         tmp_folder=tmp_folder,
                         overwrite=snap_settings['snap_overwrite'],
                         qc_file=snap_settings['snap_qc_file'],
                         verbose=verbose)

        else:
            input_bam = general_settings['output_bam'] + '/' + align_settings['sample'] + '.bam'
            output_snap = general_settings['output_snap'] + '/' + align_settings['sample'] + '.snap'

            snap_pre(input_file=input_bam,
                     output_snap=output_snap,
                     genome_name=index_settings['genome_prefix'],
                     genome_size=snap_settings['genome_size'],
                     min_mapq=snap_settings['snap_min_mapq'],
                     min_flen=snap_settings['snap_min_flen'],
                     max_flen=snap_settings['snap_max_flen'],
                     min_cov=snap_settings['snap_min_cov'],
                     max_num=snap_settings['snap_max_num'],
                     barcode_file=snap_settings['snap_barcode_file'],
                     keep_chrm=snap_settings['snap_keep_chrm'],
                     keep_single=snap_settings['snap_keep_single'],
                     keep_secondary=snap_settings['snap_keep_secondary'],
                     keep_discordant=snap_settings['snap_keep_discordant'],
                     tmp_folder=tmp_folder,
                     overwrite=snap_settings['snap_overwrite'],
                     qc_file=snap_settings['snap_qc_file'],
                     verbose=verbose)
    else:
        print('Build Snap: SKIPPED')

    if general_settings['add_peak_matrix']:
        print('Adding peak-matrix')
        if align_settings['multi']:
            samples_dict = align_settings['samples_dict']
            for sample in samples_dict.keys():

                peak_file = general_settings['output_peaks'] + '/' + sample + '_peaks.bed'
                snap_file = general_settings['output_snap'] + '/' + sample + '.snap'
                try:
                    snap_pmat(snap_file=snap_file,
                              peak_file=peak_file,
                              buffer_size=snap_settings['pmat_buffer_size'],
                              tmp_folder=tmp_folder,
                              verbose=verbose)
                except Exception as e:
                    print(e)
                    print('Error: Could not add peak matrix to snap file')
                    print('Error adding peak matrix: ' + sample + '.snap')

        else:
            peak_file = general_settings['output_peaks'] + '/' + align_settings['sample'] + '_peaks.bed'
            snap_file = general_settings['output_snap'] + '/' + align_settings['sample'] + '.snap'
            try:
                snap_pmat(snap_file=snap_file,
                          peak_file=peak_file,
                          buffer_size=snap_settings['pmat_buffer_size'],
                          tmp_folder=tmp_folder,
                          verbose=verbose)
            except Exception as e:
                print('Error adding peak matrix: ' + align_settings['sample'] + '.snap')
                print(e)

    if general_settings['add_bin_matrix']:
        print('Adding bin-matrix')

        if align_settings['multi']:
            samples_dict = align_settings['samples_dict']
            for sample in samples_dict.keys():

                snap_file = general_settings['output_snap'] + '/' + sample + '.snap'
                try:
                    snap_bmat(snap_file=snap_file,
                              bin_size=snap_settings['bin_size_list'],
                              tmp_folder=tmp_folder,
                              verbose=verbose)
                except Exception as e:
                    print('Error adding bin matrix: ' + sample + '.snap')
                    print(e)

        else:
            snap_file = general_settings['output_snap'] + '/' + align_settings['sample'] + '.snap'
            try:
                snap_bmat(snap_file=snap_file,
                          bin_size=snap_settings['bin_size_list'],
                          tmp_folder=tmp_folder,
                          verbose=verbose)
            except Exception as e:
                print('Error adding bin matrix: ' + align_settings['sample'] + '.snap')
                print(e)

    if general_settings['convert_to_anndata']:
        print('Converting Snap to Anndata')
        if align_settings['multi']:
            samples_dict = align_settings['samples_dict']
            for sample in samples_dict.keys():
                snap_file = general_settings['output_snap'] + '/' + sample + '.snap'
                basename = cv.make_rds(snap_path=snap_file)

                output_anndata = general_settings['output_anndata'] + '/' + sample + '.h5ad'
                cv.snap2anndata(rds_file_name=(basename+'.rds'),
                                path=general_settings['output_rds'],
                                mtx_name=converter_settings['mtx_name'],
                                feature_names=converter_settings['feature_name'],
                                jaccard=converter_settings['jaccard'],
                                jaccard_key=converter_settings['jaccard_key'],
                                save_all_rds=converter_settings['save_all_rds'],
                                verbose=verbose,
                                copy=converter_settings['copy'],
                                save=output_anndata
                                )
        else:
            snap_file = general_settings['output_snap'] + '/' + align_settings['sample'] + '.snap'
            rds_file = general_settings['output_rds'] + '/' + align_settings['sample'] + '.rds'
            basename = cv.make_rds(snap_path=snap_file, sample=align_settings['sample'])
            output_anndata = general_settings['output_anndata'] + '/' + align_settings['sample'] + '.h5ad'

            cv.snap2anndata(rds_file_name=(basename+'.rds'),
                            path=general_settings['output_rds'],
                            mtx_name=converter_settings['mtx_name'],
                            feature_names=converter_settings['feature_name'],
                            jaccard=converter_settings['jaccard'],
                            jaccard_key=converter_settings['jaccard_key'],
                            save_all_rds=converter_settings['save_all_rds'],
                            verbose=verbose,
                            copy=converter_settings['copy'],
                            save=output_anndata
                            )
    else:
        print('Convert to AnnData: SKIPPED')

def has_CB_tag(bam, tmp_files):
    """Check if the CB tag is present in the bam file"""
    read = next(pysam.AlignmentFile(bam))
    if read.has_tag("CB"):
        print("Bamfile has tag: 'CB'")
    else:
        bam = add_CB_tag(bam)
    return bam, tmp_files

def add_CB_tag(bam,
               output_suffix='_CB_tagged',
               replace=True):
    """
    Add a CB tag to a bam file aligned by bwa.
    The barcode is therefore extracted from the first column of the reads and the CB tag is added to the read afterwards.

    :param bam:
    :param output_suffix:
    :param replace:
    :return:
    """

    print('Adding CB Tag')

    path_split = os.path.split(bam)
    output_präfix = path_split[1].split('.bam')[0]
    output = path_split[0] + '/' + output_präfix + output_suffix + '.bam'

    bam_obj = open_bam(bam, "rb", verbosity=0)
    out_bam = pysam.AlignmentFile(output, "wb", template=bam_obj)

    for read in bam_obj:
        if read.has_tag("CB"):
            os.remove(output)
            raise Exception("bamfile already has CB tags!")
        barcode = read.qname.split(":")[0].upper();
        read.tags += [('CB', barcode)]
        out_bam.write(read)

    bam_obj.close()
    out_bam.close()

    return output

def sort_bam(bam_file, output_file, verbose, by_name=False):
    """
    Sort bam file
    :param bam_file: bam file
    :param output_path: output path
    :param overwrite: overwrite
    :param verbose: verbose
    :return: None
    """
    if verbose:
        print('Sorting bam file: ' + bam_file)
    if by_name:
        pysam.sort('-n', '-o', output_file , bam_file)
    else:
        pysam.sort('-o', output_file , bam_file)


def index_bam(bam_file, verbose):
    """
    Index bam file
    :param bam_file: bam file
    :param output_path: output path
    :param overwrite: overwrite
    :param verbose: verbose
    :return: None
    """
    if verbose:
        print('Indexing bam file: ' + bam_file)

    pysam.index(bam_file)


def collate_bam(bam_file, output_bam, verbose):

    if verbose:
        print(f'Collate bam file: {bam_file}')

    pysam.collate(bam_file, output_bam)


def fixmate(bam_file, output_bam, verbose):

    if verbose:
        print(f'Fixmate bam file: {bam_file}')

    pysam.fixmate('-m', bam_file, output_bam)


def markdup(bam_file, output_bam, threads, verbose):

    if verbose:
        print(f'Markdup bam file: {bam_file}')

    pysam.markdup('-r', '--barcode-tag', 'CB', '--threads', str(threads), bam_file, output_bam)

def process_alignment(general_settings, align_settings):
    """
    Generate fragment file and markdup
    :return: None
    """

    print('Generating fragment file')

    if align_settings['multi']:
        samples_dict = align_settings['samples_dict']
        for sample in samples_dict.keys():
            bam_processing(general_settings, sample)


    else:
        sample = align_settings['sample']
        bam_processing(general_settings, sample)


def bam_processing(general_settings, sample):
    """
    Process Bamfile to generate Fragments and filter PCR duplicates
    """
    tmp_files = []

    bam_file = general_settings['output_bam'] + '/' + sample + '.bam'
    print(f'starting bam processing. final output: {bam_file}')
    cb_tagged_file, tmp_files = has_CB_tag(bam_file, tmp_files)  # Check if CB tag is present if not add it
    tmp_files.append(cb_tagged_file)
    # output file collate
    collated_bam = general_settings['output_bam'] + '/' + 'collated_' + sample
    tmp_files.append(collated_bam + '.bam')
    # output file fixmate
    fixmate_bam = general_settings['output_bam'] + '/' + 'fixmate_' + sample + '.bam'
    tmp_files.append(fixmate_bam)
    # output file sort
    sorted_bam = general_settings['output_bam'] + '/' + 'sorted_' + sample + '.bam'
    tmp_files.append(sorted_bam)
    # output file markdup
    markdup_bam = general_settings['output_bam'] + '/' + 'markdup_' + sample + '.bam'

    # output region sorted
    region_sort_bam_file = general_settings['output_bam'] + '/' + 'region_sorted' + sample + '.bam'

    # output fragments
    fragment_file = general_settings['output_bam'] + '/' + 'fragments_'+ sample + '.bed'

    # collate
    print(f'collating: {cb_tagged_file}')
    collate_bam(bam_file=cb_tagged_file,
                output_bam=collated_bam,
                verbose=general_settings['verbose'])

    # fixmate
    input_fixmate = collated_bam + '.bam'
    print(f'fixmate on: {input_fixmate}')
    fixmate(bam_file=input_fixmate,
            output_bam=fixmate_bam,
            verbose=general_settings['verbose'])

    # sort
    print(f'sorting: {fixmate_bam}')
    sort_bam(bam_file=fixmate_bam,
             output_file=sorted_bam,
             verbose=general_settings['verbose'])

    # markdup
    print(f'markdup: {sorted_bam}')
    markdup(bam_file=sorted_bam,
            output_bam=markdup_bam,
            threads=general_settings['n_threads'],
            verbose=general_settings['verbose'])

    print('exchange old bam with markdup...')
    os.remove(bam_file)
    os.rename(markdup_bam, region_sort_bam_file)

    print(f'indexing: {region_sort_bam_file}')
    index_bam(bam_file=region_sort_bam_file,
              verbose=general_settings['verbose'])

    print(f'sinto making fragments (input: {region_sort_bam_file}, output: {fragment_file})')
    fragments.fragments(bam=region_sort_bam_file,
                        fragment_path=fragment_file)

    print('removing temporary files')
    for file in tmp_files:
        print(f'try removing {file}')
        try:
            os.remove(file)
        except Exception as e:
            print(f'Error while cleaning up temporaray files {e}')

    # sort by read name for peakcalling
    print(f'sorting by name: {region_sort_bam_file}')
    sort_bam(bam_file=region_sort_bam_file,
             output_file=bam_file,
             verbose=general_settings['verbose'],
             by_name=True)



def auto_path(general_settings):

    dirs = []

    if not general_settings['custom_output']:

        output_directory = general_settings['output_directory']
        general_settings['output_bam'] = os.path.join(output_directory, 'bamfiles')
        dirs.append(general_settings['output_bam'])
        general_settings['output_peaks'] = os.path.join(output_directory, 'peaks')
        dirs.append(general_settings['output_peaks'])
        general_settings['output_snap'] = os.path.join(output_directory, 'snap')
        dirs.append(general_settings['output_snap'])
        general_settings['output_rds'] = os.path.join(output_directory, 'rds')
        dirs.append(general_settings['output_rds'])
        general_settings['output_anndata'] = os.path.join(output_directory, 'anndata')
        dirs.append(general_settings['output_anndata'])

    else:
        dirs.append(general_settings['bamfiles'])
        dirs.append(general_settings['output_peaks'])
        dirs.append(general_settings['output_snap'])
        dirs.append(general_settings['output_rds'])
        dirs.append(general_settings['output_anndata'])

    makeDir(dirs)

    return general_settings

def makeDir(to_build):
    '''
    Method to make directories given by a list of paths
    :param to_build: list of str
    :return: None
    '''

    for path in to_build:
        if path != None and not os.path.isdir(path):
            try:
                Path(path).mkdir(parents=True)
                print(path + ': NEWLY SETUP')

            except Exception as e:
                print(e)

def call_peak(output_path,
              output_prefix,
              bam,
              print,
              tmp_folder=None,
              path_to_macs=None,
              macs_options=None,
              ):
    """
    Call peaks on selected bamfiles

    :param output_path:
    :param bams: list of str
        List of bamfiles to call peaks on
    :param tmp_folder:
        temporary directory to store intermediate files
    :param path_to_macs: str
        Path to MACS2 peakcalling software
    :param macs_options: list of str
        a list of strings indicating options you'd like passed to aligner.
        (default: "--nomodel --qval 1e-2 -B --SPMR --call-summits --keep-dup all");
    :return:
    """

    # if the path_to_macs path given, need to check the existance of MACS1
    if path_to_macs != None:
        path_to_macs += "/"
        if not os.path.isdir(path_to_macs):
            print(('Error: ' + path_to_macs + ' is not a folder'));
            sys.exit(1);
        if not os.path.exists(path_to_macs + "macs2"):
            print('Error: macs2 does not exist')
            sys.exit(1);
    else:
        try:
            # pipe output to /dev/null for silence
            null = open("/dev/null", "w")
            subprocess.Popen("macs2", stdout=null, stderr=null)
            null.close()
        except OSError as e:
            print('Error: macs2 does not exist!');
            sys.exit(1);
        path_to_macs = ""

    # check temp folder
    if(tmp_folder!=None):
        if not os.path.isdir(tmp_folder):
            print("Error: 'tmp_folder' is not a folder or does not exist")
            sys.exit(1);

    # default aligner option
    if macs_options is None:
        #macs_options = ["--nomodel", "--qval 1e-2", "-B", "--SPMR", "--call-summits", "--keep-dup all"];
        macs_options = []

    print(f'calling MACS2 on: {bam}')
    # call peaks using macs
    args = [path_to_macs + "macs2"];
    args.append("callpeak");
    args.append("-f BAM");
    args.append("-t " + bam);
    args.append("-n " + output_prefix);
    # args.append("-g " + gsize);
    args.append("--outdir " + output_path)
    args.extend(macs_options);
    #ftmp = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder)

    try:
        subprocess.check_call(" ".join(args), stdout=None, shell=True, executable='/bin/bash');
    except subprocess.CalledProcessError as e:
        sys.exit('error: fail to run macs2!');

def open_bam(file, mode, verbosity=3, **kwargs):
    """
    Open bam file with pysam.AlignmentFile. On a specific verbosity level.

    Parameters
    ----------
    file : str
        Path to bam file.
    mode : str
        Mode to open the file in. See pysam.AlignmentFile
    verbosity : int, default 3
        Set verbosity level. Verbosity level 0 for no messages.
    **kwargs :
        Forwarded to pysam.AlignmentFile

    Returns
    -------
    pysam.AlignmentFile :
        Object to work on SAM/BAM files.
    """

    # save verbosity, then set temporary one
    former_verbosity = pysam.get_verbosity()
    pysam.set_verbosity(verbosity)

    # open file
    handle = pysam.AlignmentFile(file, mode, **kwargs)

    # return to former verbosity
    pysam.set_verbosity(former_verbosity)

    return handle

if __name__ == '__main__':

    where_am_i = os.path.dirname(os.path.realpath(__file__))
    where_am_i = os.path.split(where_am_i)[0]
    general_settings = {}
    general_settings['output_directory'] = where_am_i + '/data'
    general_settings['custom_output'] = False
    general_settings['output_bam'] = where_am_i + '/data/bamfiles'
    general_settings['output_peaks'] = where_am_i + '/data/peaks'
    general_settings['output_snap'] = where_am_i + '/data/snap'
    general_settings['output_rds'] = where_am_i + '/data/rds'
    general_settings['output_anndata'] = where_am_i + '/data/anndata'
    general_settings['verbose'] = True
    general_settings['tmp_folder'] = './'
    general_settings['n_threads'] = 8
    general_settings['index_genome'] = False
    general_settings['align'] = False
    general_settings['callpeak'] = False
    general_settings['build_snap'] = False
    general_settings['add_peak_matrix'] = False
    general_settings['add_bin_matrix'] = False
    general_settings['convert_to_anndata'] = False

    general_settings['create_fragments_bed'] = True

    index_settings = {}
    index_settings['genome'] = where_am_i + '/genomes/hg38/hg38.fa'
    index_settings['genome_prefix'] = 'hg38'
    index_settings['aligner'] = 'bwa'
    index_settings['path_to_aligner'] = '/app/bwa'

    align_settings = {}
    align_settings['fastq_dir'] = where_am_i + '/data/fastq'
    align_settings['multi'] = False
    align_settings['sample'] = 'Esophagus'
    align_settings['input_fastq1'] = where_am_i + '/data/fastq/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R1.fastq.gz'
    align_settings['input_fastq2'] = where_am_i + '/data/fastq/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R2.fastq.gz'
    align_settings['samples_dict'] = samples_dict = {}
    align_settings['read_fastq_command'] = 'zcat'
    align_settings['min_cov'] = 0
    align_settings['if_sort'] = True
    align_settings['overwrite'] = True

    peakcalling_settings = {}
    peakcalling_settings['path_to_macs'] = None
    peakcalling_settings['macs_options'] = None
    peakcalling_settings['optional_input'] = None
    peakcalling_settings['optional_bam'] = where_am_i + '/data/bam/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.bam'

    snap_genome_size = where_am_i + '/genomes/hg38/hg38.gs'
    snap_min_mapq = 10
    snap_min_flen = 0
    snap_max_flen = 1000
    snap_min_cov = 0
    snap_max_num = 1000000
    snap_barcode_file = None
    snap_keep_chrm = True
    snap_keep_single = True
    snap_keep_secondary = True
    snap_keep_discordant = False
    snap_overwrite = True
    snap_qc_file = None
    pmat_buffer_size = 1000
    bin_size_list = [5000]

    snap_settings = {}
    snap_settings['genome_size'] = snap_genome_size
    snap_settings['snap_min_mapq'] = snap_min_mapq
    snap_settings['snap_min_flen'] = snap_min_flen
    snap_settings['snap_max_flen'] = snap_max_flen
    snap_settings['snap_min_cov'] = snap_min_cov
    snap_settings['snap_max_num'] = snap_max_num
    snap_settings['snap_barcode_file'] = snap_barcode_file
    snap_settings['snap_keep_chrm'] = snap_keep_chrm
    snap_settings['snap_keep_single'] = snap_keep_single
    snap_settings['snap_keep_secondary'] = snap_keep_secondary
    snap_settings['snap_keep_discordant'] = snap_keep_discordant
    snap_settings['snap_overwrite'] = snap_overwrite
    snap_settings['snap_qc_file'] = snap_qc_file
    snap_settings['pmat_buffer_size'] = pmat_buffer_size
    snap_settings['bin_size_list'] = bin_size_list

    converter_settings = {}
    converter_settings['mtx_name'] = 'pmat'
    converter_settings['feature_name'] = 'peak'
    converter_settings['metadata_name'] = 'metaData'
    converter_settings['jaccard'] = False
    converter_settings['jaccard_key'] = 'jmat'
    converter_settings['save_all_rds'] = False
    converter_settings['copy'] = True


    preprocessing(index_settings, align_settings, peakcalling_settings, snap_settings, converter_settings, general_settings)

    '''
    snaptools align-paired-end --input-reference=/home/jan/python-workspace/sc-atac/data/genome/hg38/hg38 --input-fastq1=/home/jan/python-workspace/sc-atac/data/fastq/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R1.fastq.gz --input-fastq2=/home/jan/python-workspace/sc-atac/data/fastq/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R2.fastq.gz --output-bam=/home/jan/python-workspace/sc-atac/data/bam/output.bam --aligner=bwa --path-to-aligner=/home/jan/python-workspace/sc-atac/biotools/bwa/bin/ --read-fastq-command=zcat --min-cov=0 --num-threads=8 --if-sort=True --tmp-folder=./ --overwrite=True
    ./bwa mem -t 8 /home/jan/python-workspace/sc-atac/data/genome/hg38/hg38.fa /home/jan/python-workspace/sc-atac/data/fastq/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R1.fastq.gz /home/jan/python-workspace/sc-atac/data/fastq/ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R2.fastq.gz 

    '''