#!/bin/bash
#Syntax fastqToSnap dir/R1.fastq dir/R2fastq dir/anotherR1.fastq dir/anotherR2fastq .....

#06.07.2022
#Input(names of the fastq files):
path_array=(\
#"ENC-1JKYN-012-SM-JF1O6_snATAC_body_of_pancreas_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1JKYN-012-SM-JF1O6_snATAC_body_of_pancreas_Rep1.demultiplexed.R2.fastq.gz" \
"ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta_Rep1.demultiplexed.R1.fastq.gz" \
"ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta_Rep1.demultiplexed.R2.fastq.gz" \
"ENC-1JKYN-066-SM-CSSD4_snATAC_omental_fat_pad_Rep1.demultiplexed.R1.fastq.gz" \
"ENC-1JKYN-066-SM-CSSD4_snATAC_omental_fat_pad_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1JKYN-118-SM-IQYCP_snATAC_lower_leg_skin_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1JKYN-118-SM-IQYCP_snATAC_lower_leg_skin_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1JKYN-146-SM-A8CPH_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1LVAN-288-SM-JF1O3_snATAC_stomach_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1LVAN-288-SM-JF1O3_snATAC_stomach_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1LGRB-229-SM-IOBHV_snATAC_stomach_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1LGRB-229-SM-IOBHV_snATAC_stomach_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1JKYN-179-SM-ACCPU_snATAC_upper_lobe_of_left_lung_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1JKYN-179-SM-ACCPU_snATAC_upper_lobe_of_left_lung_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1JKYN-191-SM-JF1O8_snATAC_colon_sigmoid_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1JKYN-191-SM-JF1O8_snATAC_colon_sigmoid_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1JKYN-206-SM-ACCQ1_snATAC_colon_transverse_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1JKYN-206-SM-ACCQ1_snATAC_colon_transverse_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1K2DA-037-SM-JF1NS_snATAC_body_of_pancreas_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1K2DA-037-SM-JF1NS_snATAC_body_of_pancreas_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1K2DA-044-SM-A62E9_snATAC_upper_lobe_of_left_lung_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1K2DA-044-SM-A62E9_snATAC_upper_lobe_of_left_lung_Rep1.demultiplexed.R2.fastq.gz" \
#"ENC-1K2DA-063-SM-CTD24_snATAC_gastroesophageal_sphincter_Rep1.demultiplexed.R1.fastq.gz" \
#"ENC-1K2DA-063-SM-CTD24_snATAC_gastroesophageal_sphincter_Rep1.demultiplexed.R2.fastq.gz" \
"ENC-1K2DA-070-SM-AZPYJ_snATAC_esophagus_squamous_epithelium_Rep1.demultiplexed.R1.fastq.gz" \
"ENC-1K2DA-070-SM-AZPYJ_snATAC_esophagus_squamous_epithelium_Rep1.demultiplexed.R2.fastq.gz" \
"ENC-1K2DA-084-SM-IQYD1_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R1.fastq.gz" \
"ENC-1K2DA-084-SM-IQYD1_snATAC_esophagus_muscularis_mucosa_Rep1.demultiplexed.R2.fastq.gz" \
);
echo $path_array

#Parameter
#To use a different directory for the fastq files simply uncomment the fastq_dir variable  
#fastq_dir=/some/other/directory
#data directory
data_dir=/home/rstudio/data
#reference genome:
genome_input=/home/rstudio/genomes/hg38/hg38.fa
genome_output_prefix=hg38
#Aligner:
aligner=bwa
path_to_aligner=/home/rstudio/biotools/bwa
#input by command line:
commandLineInput=false
#index Genome:
index_genome=true
#perform alignment:
perform_alignment=true
#build snapfile
build_snapfile=true
#perform peakcalling
run_macs=true
#add PMAT to snap
snapAddPMAT=true
#add BMAT to snap
snapAddBMAT=true
#convert snap to anndata
snapToAnndata=true

#----------------------------------------------------------------------------------------
#Check if the data directory exits and ask if it should be created
if [ ! -d $data_dir ]
then 
    echo "Directory '"data_dir"' does not exist"
    while true; do
        read -p "Should the directory be created? [y/n] " yn
        case $yn in
            [Yy]* ) mkdir -pv $data_dir/{peaks/{peaks_bed,output_macs},anndata,snap,bam,fastq,analysis,}; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi
#----------------------------------------------------------------------------------------
#Setting Arguments by command line
if $commandLineInput
then
i=0;
path_array=()
#set -A path_array;
for file in "$@"
do
   echo "FastQ - $i: $file";
   path_array[$i]=$file;
   i=$((i + 1));
done
echo ""
#Argument Check
if !(((i+1) % 2))
then
	echo 'Only pairs as argument valid'
	exit $E_BADARGS
fi
fi
#-----------------------------------------------------------------------------------------
#Index reference genome
if $index_genome
then 
 snaptools index-genome  \
 	--input-fasta=$genome_input  \
 	--output-prefix=$genome_output_prefix  \
 	--aligner=$aligner  \
 	--path-to-aligner=$path_to_aligner  \
 	--num-threads=8
fi
#-----------------------------------------------------------------------------------------
#build paths and loop through files
input_reference=$genome_input

output_bam_dir=$data_dir"/bam"

output_dir_pre_snap=$data_dir"/snap"

#TODO anndata path

c=0;
count_alignments=0;

for ((c>0 ; c<${#path_array[@]} ; c=c+2 ));
do
if [ -z "$fastq_dir" ]
then
path_1=$data_dir"/fastq/"${path_array[c]};
path_2=$data_dir"/fastq/"${path_array[(c+1)]};
else
path_1=$fastq_dir"/"${path_array[c]};
path_2=$fastq_dir"/"${path_array[(c+1)]};
fi
filename_1=$(basename $path_1 .demultiplexed.R1.fastq.gz)
filename_2=$(basename $path_2 .demultiplexed.R2.fastq.gz)


echo "Path to R1 FastQ: "$path_1;
echo "Filename R1: "$filename_1;
echo ""
echo "Path to R2 FastQ: "$path_2;
echo "Filename R2: "$filename_2;
echo ""
output_path_bam=$output_bam_dir/$filename_1.bam

echo "Output Path Alignment: "$output_path_bam
echo ""
#Align command
if $perform_alignment
then
 snaptools align-paired-end  \
 	--input-reference=$input_reference  \
 	--input-fastq1=$path_1  \
 	--input-fastq2=$path_2  \
 	--output-bam=$output_path_bam  \
 	--aligner=$aligner  \
 	--path-to-aligner=$path_to_aligner  \
 	--read-fastq-command=zcat  \
 	--min-cov=0  \
 	--num-threads=8  \
 	--if-sort=True  \
 	--tmp-folder=./  \
 	--overwrite=TRUE    
fi
#---------------------------------------------------------------------------------------
if $run_macs
then
echo "Test run_macs: "
path_peaks_bed=$(bash run_macs.sh $data_dir"/peaks" $output_path_bam)
echo $path_peaks_bed
fi
path_peaks_bed=$data_dir"/peaks/peaks_bed/"$filename_1\_peaks.bed
echo $path_peaks_bed
#---------------------------------------------------------------------------------------
#Build snap file
echo ""
echo "snap directory: "
output_path_pre_snap=$output_dir_pre_snap/$filename_1.snap
echo "Path to the snapfile: "$output_path_pre_snap
if $build_snapfile
then
echo "building snap"
#Full command

#snaptools snap-pre  \
#	--input-file=$output_path_bam  \
#	--output-snap=$output_path_pre_snap  \
#	--genome-name=$genome_output_prefix  \
#	--genome-size=$genome_output_prefix.chrom.sizes  \
#	--min-mapq=30  \
#	--min-flen=0  \
#	--max-flen=1000  \
#	--keep-chrm=TRUE  \
#	--keep-single=FALSE  \
#	--keep-secondary=FALSE  \
#	--overwrite=True  \
#	--min-cov=100  \
#	--verbose=True

#minimum Arguments Variant:

 snaptools snap-pre \
 	--input-file=$output_path_bam \
 	--output-snap=$output_path_pre_snap \
 	--genome-name=$genome_output_prefix \
 	--genome-size=$data_dir"/genome/"$genome_output_prefix".gs"
fi
#--------------------------------------------------------------------------------------
if $snapAddPMAT
then
echo "add PMAT to snap"
echo $output_path_pre_snap
echo $path_peaks_bed
echo ""
bash snapAddPMAT.sh $output_path_pre_snap $path_peaks_bed
fi
#--------------------------------------------------------------------------------------
if $snapAddBMAT
then 
echo "add BMAT to snap"
snaptools snap-add-bmat --snap-file $output_path_pre_snap
echo ""
fi

if $snapToAnndata
then
echo "executing converter.py"
python3 converter.py $output_path_pre_snap
echo "done"
echo ""
fi

count_alignments=$((count_alignments + 1));
echo "Number of runs performed: "  $count_alignments;
done
