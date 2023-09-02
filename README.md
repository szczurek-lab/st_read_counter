# read_counter

Since standard mutation callers cannot be directly applied to Spatial Transcriptomics BAM files, we developed a custom tool called read_counter (https://github.com/szczurek-lab/st_read_counter) to identify mutation alleles for input to Tumoroscope.  


# Usage 

```python st_read_counter.py st_bam_file  st_bai_file  barcode  vcf_file  output_name  tumor_column_vcf st_dir```

The `st_read_counter.py` script is designed to streamline the process of analyzing Spatial Transcriptomics (ST) data to get the allelic read counts. It requires several input files and parameters to operate effectively:

1. `st_bam_file` and `st_bai_file`: These are the BAM and BAI files containing the ST data.

2. `barcode`: This parameter represents the barcode or identifier associated with the st sample, which is crucial for identifying and distinguishing the position of the data points.

3. `vcf_file`: The VCF (Variant Call Format) file contains information about selected mutations from Whole Exome Sequencing (WES). 

4. `output_name`: The name of the output file.

5. `tumor_column_vcf`: This specifies the column in the VCF file that corresponds to the sample.

6. `st_dir`: This parameter defines the directory where the program will generate the necessary output files.

The script's primary function is to process the provided ST BAM files, integrate them with the relevant VCF data, and generate allele counts for each potential variant site. These allele counts serve as essential inputs for downstream analysis tools like Tumoroscope, enabling the inference of subclonal populations across the spatial spots.

To gain comprehensive usage details and further insights into the functionality of `st_read_counter.py`, please refer to the `read_counter` repository and associated documentation.
