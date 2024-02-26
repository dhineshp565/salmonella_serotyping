# salmonella_serotyping
Salmonella enterica whole genome assembly and serotyping for Oxford Nanopore reads
Uses dragonflye for whole genome assembly (https://github.com/rpetit3/dragonflye) and  SISTR (https://github.com/phac-nml/sistr_cmd) for serotype prediction

## Usage
```
nextflow run main.nf --input path/to/input --out_dir path/to/output
```
## Parameters
```
--input Path to directory containing subdirectories with fastq files
--out_dir Path to output directory
```

## Depdencies
### - Nextflow
### - Docker
