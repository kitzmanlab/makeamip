### Setup python environment

Install conda
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
# answer questions 
```

Set channel preference
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Set up python environment
```
conda create --name makeamip python=2 future pysam pybedtools cython biopython openpyxl scipy scikit-learn pytables primer3-py bzip2=1.0.6=1 matplotlib

pip install pyvcf

# now, activate this environment
source actviate makeamip
```

### Setup: references

Set env variable to point to location where reference & annotations are staged
```
export REFS=/nfs/kitzman2/jacob/proj/mip_design_pipe/annots_and_files/

```
 
- *${REFS}/indices/* - alignment and jellyfish index files

- *${REFS}/annots/* - annotation files 
	- hg19 knownGene table
	- hg19 refGene table
	- gencode V24 lifted over to hg19 table
	- 1000 Genome SNP variant file

- *${REFS}/refs/* - reference genome sequences

- *${REFS}/parameters/* - parameters for MIP probe scoring model

### Prereq software

Should be in $PATH:
- mrFAST 2.6.1.0
- Jellyfish >= 2.1.4
- bwa >= 0.7.12
- samtools

Also, picard should be installed, with an environmental variable PICARD_DIR pointing to the path where the picard.jar file lives.

### Usage 
See pipeline_script/example_usage.sh