#!/bin/bash

#Maintainer Jan Detleffsen
path_array_snapfiles=(\
"/mnt/workspace/jdetlef/snapATAC_Analysis/snapFiles/ENC-1LVAN-288-SM-JF1O3_snATAC_stomach_Rep1.demultiplexed.snap" \
"/mnt/workspace/jdetlef/snapATAC_Analysis/snapFiles/ENC-1LGRB-229-SM-IOBHV_snATAC_stomach_Rep1.demultiplexed.snap" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/snapFiles/ENC-1JKYN-161-SM-JF1NP_snATAC_stomach_Rep1.demultiplexed.snap" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-066-SM-CSSD4_snATAC_omental_fat_pad_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-133-SM-IOBHO_snATAC_heart_left_ventricle_Rep1.demultiplexed.bam" \
);

path_array_peakfiles=(\
"/mnt/workspace/jdetlef/snapATAC_Analysis/run_macs2/macs_out/bedfiles/ENC-1LVAN-288-SM-JF1O3_snATAC_stomach_Rep1.demultiplexed_peaks.bed" \
"/mnt/workspace/jdetlef/snapATAC_Analysis/run_macs2/macs_out/bedfiles/ENC-1LGRB-229-SM-IOBHV_snATAC_stomach_Rep1.demultiplexed_peaks.bed" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/run_macs2/macs_out/bedfiles/ENC-1JKYN-161-SM-JF1NP_snATAC_stomach_Rep1.demultiplexed_peaks.bed" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1JKYN-066-SM-CSSD4_snATAC_omental_fat_pad_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver_Rep1.demultiplexed.bam" \
#"/mnt/workspace/jdetlef/snapATAC_Analysis/bamfiles/ENC-1LGRB-133-SM-IOBHO_snATAC_heart_left_ventricle_Rep1.demultiplexed.bam" \
);

#Parameters

#Use command line input:
commandLineInput=true

#----------------------------------------------------------------------------------------
#Setting Arguments by command line
if $commandLineInput;
then
    path_array_snapfiles=()
    path_array_peakfiles=()
    i=0;
    set -A input_array;
    for file in "$@"
        do
        echo "Path - $i: $file";
        input_array[$i]=$file;
        i=$((i + 1));
        done
    #Argument Check
    if !(((i+1) % 2));
    then
        echo 'Only pairs as argument valid'
        exit $E_BADARGS
    fi

    #sort paths of the peakfiles and snapfiles
    #Alignment
    j=0;
    k=0;
    #for k in {0..5..2}
    for ((j>0 ; j<${#input_array[@]} ; j=j+2 ));
        do
        path_array_snapfiles[k]+="${input_array[j]}"
        path_array_peakfiles[k]+="${input_array[(j+1)]}"
        echo ${path_array_snapfiles[k]}
        echo ${path_array_peakfiles[k]}
        k=$((k + 1));
        done
fi

#----------------------------------------------------------------------------------------
#Add PMAT to snap-file
echo "Add PMAT to files:"
n_snapfiles=${#path_array_snapfiles[@]}
n_peakfiles=${#path_array_peakfiles[@]}
echo "number of snapfiles:  "$n_snapfiles
echo "number of peakfiles:  "$n_peakfiles

if [ $n_snapfiles == $n_peakfiles ]
then

c=0;
for ((c>0 ; c<$n_snapfiles ; c=c+1 ));
do
path_snap=${path_array_snapfiles[$c]};
path_peaks=${path_array_peakfiles[$c]};
echo "executing snaptools snap-add-pmat with:"
echo "Snapfile: "$path_snap;
echo "Peakfile: "$path_peaks;
snaptools snap-add-pmat --snap-file $path_snap --peak-file $path_peaks
done
else
echo "unequal number of snap- and peakfiles"
fi
