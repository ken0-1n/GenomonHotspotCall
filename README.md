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
$ hotspotCall -h
usage: hotspotCall [-h] [--version] [-S samtools_params] [-t min_tumor_misrate] [-c max_control_misrate] [-R T/N ratio_control] [-m min_lod_score] [-r rna.bam] tumor.bam control.bam output_file hotspot_mutations.bed

positional arguments:
  tumor.bam             the path to the tumor bam file
  control.bam           the path to the control bam file
  output_file           the path to the output file
  hotspot_mutations.bed
                        the bed format file that lists mutations

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -S samtools_params
  -t min_tumor_misrate
  -c max_control_misrate
  -R T/N ratio_control
  -m min_lod_score
  -r rna.bam            the path to the RNA bam file
```

