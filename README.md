# salmonella_serotyping
Salmonella enterica whole genome assembly and serotyping for Oxford Nanopore reads
Uses dragonflye for whole genome assembly (https://github.com/rpetit3/dragonflye) and  SISTR (https://github.com/phac-nml/sistr_cmd) for serotype prediction

## Usage
```
nextflow run main.nf --input path/to/input --out_dir path/to/output
```
## Required arguments
```
--input Path to directory containing subdirectories with fastq files
--out_dir Path to output directory
```

## Dependencies
* nextflow
* docker
* wsl2

## Docker images and versions used in this pipeline
* porechop - quay.io/biocontainers/porechop:0.2.4--py39h1f90b4d_6
* dragonflye - quay.io/biocontainers/dragonflye:1.1.1--hdfd78af_0
* sistr - staphb/sistr:1.1.1
* busco - nanozoo/busco:5.5.0--427a3e7
* make_report - dhp565/rmd_ampliseq:21022024
