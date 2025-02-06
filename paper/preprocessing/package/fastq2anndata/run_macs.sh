#!/bin/bash
#Syntax run_macs.sh destination/ dir/*bam dir/*bam dir/*bam dir/*bam .....
#06.07.2022
#Input:

#path_array=(\
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LVAN-288-SM-JF1O3_snATAC_stomach_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-229-SM-IOBHV_snATAC_stomach_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-161-SM-JF1NP_snATAC_stomach_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-066-SM-CSSD4_snATAC_omental_fat_pad_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-133-SM-IOBHO_snATAC_heart_left_ventricle_Rep1.demultiplexed.bam" \
#);

#Parameters:
output_dir="/some/dir/macs_out"
#command line input true/false
commandLineInput=true
#verbose default: false (false when called by fastqToSnap.sh)
verbose=false

#FDR
# weak:
q_value="0.01"
# q_value="0.001"
# medium:
# q_value="0.0001"
# q_value="0.00001"
# q_value="0.000001"
# q_value="0.0000001"
# strong:
# q_value="0.00000001"
# q_value="0.000000001"


#----------------------------------------------------------------------------------------
#Setting Arguments by command line
if $commandLineInput
then
i=0;
init_dir=true
path_array=()
for file in "$@"
do
   if $init_dir
   then 
   output_dir=$file
   fi
   if [ $init_dir == false ]
   then

      if $verbose
      then
      echo "BAM - $i: $file";
      fi

   path_array[$i]=$file;
   i=$((i + 1));
   fi
   init_dir=false
done
fi
#----------------------------------------------------------------------------------------
#Call peaks by MACS2
for path in "${path_array[@]}"
do
filename=$(basename $path .bam)
if $verbose
then 
echo "Call peaks on: "
echo $filename
echo ""
fi
macs2 callpeak -t $path -f BAM --outdir $output_dir/macs_out -n $filename -q $q_value
narrowPeak_file=$output_dir/macs_out/$filename\_peaks.narrowPeak
if $verbose
then
echo "File to be cropped to peaks.bed:"
echo $narrowPeak_file
echo ""
fi
peaks_bed=$output_dir/peaks_bed/$filename\_peaks.bed
epi_peaks_bed=$output_dir/peaks_bed/$filename\_epi_peaks.bed
#----------------------------------------------------------------------------------------
#Crop narrow_peak to peaks.bed
awk 'BEGIN{FS=OFS="\t"}{gsub($1, k=("b'\''"$1"'\''")); print$1,$2,$3};' $narrowPeak_file > $peaks_bed
awk 'BEGIN{FS=OFS="\t"}{print$1"_"$2"_"$3};' $narrowPeak_file > $epi_peaks_bed
echo $peaks_bed
done
