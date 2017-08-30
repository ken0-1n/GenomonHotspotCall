# GenomonHotspotCall
identify hotspot mutations


## Dependency
Python (>= 2.7), pysam, 
samtools

## Install

```
wget https://github.com/ken0-1n/GenomonHotspotCall/archive/v0.1.0.tar.gz
tar xzvf v0.1.0.tar.gz
cd GenomonHotspotCall-0.1.0/
python setup.py install
```

## Run

```
$ hotspotCall -h
usage: hotspotCall [-h] [--version] [-S samtools_params] [-t min_tumor_misrate] [-c max_control_misrate] [-R T/N ratio_control] [-m min_lod_score] [-r rna.bam] tumor.bam control.bam output_file hotspot_mutations.bed

positional arguments:
  tumor.bam              The path to the tumor bam file
  control.bam            The path to the control bam file
  output_file            The path to the output file
  hotspot_mutations.bed  The bed format file that lists mutations
                        
optional arguments:
  -h, --help             Show this help message and exit
  --version              Show program's version number and exit
  -S samtools_params     Samtools params
  -t min_tumor_misrate   The minimum amount of tumor allele frequency (default 0.1)
  -c max_control_misrate The maximum amount of control allele frequency (default 0.1)
  -R T/N ratio_control   The maximum value of the ratio between normal and tumor (default 0.1)
  -m min_lod_score       The minimum lod score (default 8.0)
  -r rna.bam             The path to the RNA bam file
```

We use this tool to create the hotspot_mutations.bed.
https://github.com/ken0-1n/GenomonHotspotDatabase

