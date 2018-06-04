#### VariantCalling Pipe 
```
Bash Script and Snakefile to call variations using GATK
```
#### Dependency

- picard
- gatk
- samtools
- snakemake

##### Usage
```
1. Run cfg.py to generate config.json for Snakefile
2. Run snakemake --cores N to do variant calling (N is required for running on lsf) 
```

###### cfg.py
```
python cfg.py wd DIR SUFFIX KEY  # add files through wildcard
python cfg.py ad KEY VALUE # add files through their paths 
```
###### snakefile
```
snakemake --cores N 
```
