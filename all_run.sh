#!/bin/bash

project="pathologist_prostate_any"
sections=([0]="1.2" [1]="2.4" [2]="3.3")
#headers=([0]="12"  [1]="24" [2]="33")
barcode="/home/shafighi/calling_from_BAM/barcodes/1000L2_barcodes_header.txt"
#selected="3.3"
#selected_col="33"
st_bam_dir="/home/shafighi/alireza/stpipeline/"
wes_bam_dir="/home/shafighi/alireza/data/bam/"
vcf_dir="/home/shafighi/calling_from_BAM/output/vardict/"
true_spots="/home/shafighi/TUMOROSCOPE/input/CellCountSummary_selected_any.txt"
st_dir="results/${project}/st_reads/"
wes_dir="results/${project}/wes_reads/"
filtered_dir="results/${project}/filtered_vcfs/"
vcf_file_selected="${filtered_dir}1_all.vcf"
tumor_column_vcf="w12"


counting_st () {
for i in "${!sections[@]}";do
tumor_column_vcf="w${sections[i]//[.]/}"
vcf_file="${vcf_dir}w${sections[i]//[.]/}_blood.vcf"
st_bam_file="${st_bam_dir}st_${sections[i]}/annotated.bam"
st_bai_file="${st_bam_dir}st_${sections[i]}/annotated.bam.bai"
output_name="vardict2_st_calls_${sections[i]}"
python st_read_counter.py $st_bam_file  $st_bai_file  $barcode  $vcf_file  $output_name  $tumor_column_vcf $st_dir
done
}

counting_wes () {
for i in "${!sections[@]}"
do
wes_bam_file="${wes_bam_dir}wes${sections[i]}_bowtie.bam"
wes_bai_file="${wes_bam_dir}wes${sections[i]}_bowtie.bam.bai"
output_name="vardict2_wes_st_p-value_${sections[i]}"
echo $i
echo $st_dir
echo $output_name
echo $wes_bam_file 
echo $wes_bai_file 
echo $barcode $vcf_file 
echo $output_name 
echo $tumor_column_vcf 
python wes_st_read_counter.py $wes_bam_file $wes_bai_file $barcode $vcf_file_selected $output_name $tumor_column_vcf $wes_dir
done
}



#bash run.sh ${sections[@]} ${headers[@]} $barcode $vcf_dir $st_bam_dir $st_dir
#python filter_vcf.py $barcode $true_spots $vcf_dir $st_dir $filtered_dir ${sections[@]}
#bash wes_run.sh $barcode $vcf_file_selected $tumor_column_vcf $wes_bam_dir $filtered_dir ${sections[@]}
counting_wes
