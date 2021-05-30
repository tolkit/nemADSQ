# nemADSQ
A repository containing pipelines for nematode genome assembly, decontamination, scaffolding, and QC.

## Requirements
Most requirements will be satisfied automatically by nextflow using the conda environment specification, but requires Red (REpeat Detector) which you can get [http://toolsmith.ens.utulsa.edu](http://toolsmith.ens.utulsa.edu). Install one of the pre-built binaries for your system or build from source. The `Red` executable needs to be in your $PATH. 


## Example run

```
# clone the repo
git clone https://github.com/tolkit/nemADSQ.git

# To install nextflow
wget -qO- https://get.nextflow.io | bash
# if no java found by nextflow:
# sudo apt install openjdk-8-jre-headless

# Get Red
wget http://toolsmith.ens.utulsa.edu/red/data/DataSet2Unix64.tar.gz
tar zxf DataSet2Unix64.tar.gz 

# add new binaries to a local bin directory
mkdir bin
export PATH=$PWD/bin:$PATH
mv nextflow redUnix64/Red bin/

# Download data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/067/145/GCA_904067145.1_BOKI2/GCA_904067145.1_BOKI2_genomic.fna.gz -O B_okinawaensis.fasta.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR243/DRR243689/DRR243689_1.fastq.gz -O B_okinawaensis.npr.fasta.gz

# Run pipeline
nextflow -C nemADSQ/qc_assem.conf run nemADSQ/qc_assem.nf --assemblies "B_okinawaensis.fasta.gz" --reads "B_okinawaensis.npr.fasta.gz" --busco2nigons nemADSQ/example/gene2Nigon_busco20200927.tsv.gz
```

Example result of *Bursaphelenchus okinawaensis* (GCA_904067145.1) assembly.
<img src="examples/B_okinawaensis.png">