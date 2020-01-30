# peterk87/nf-villumina
**Generic viral Illumina sequence analysis pipeline**

[![Build Status](https://travis-ci.org/peterk87/nf-villumina.svg?branch=master)](https://travis-ci.org/peterk87/nf-villumina)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/peterk87/nf-villumina.svg)](https://hub.docker.com/r/peterk87/nf-villumina)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2925)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with a [Singularity][] container making installation trivial and results highly reproducible.

`nf-villumina` will 

- remove low quality reads ([fastp])
- filter for reads from a taxonomic group of interest (by default superkingdom [Viruses](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10239&lvl=3&lin=f&keep=1&srchmode=1&unlock) (taxid=10239)) using [Kraken2][] and [Centrifuge][] classification results
- perform *de novo* assembly with [Unicycler] and [Shovill] on the taxonomic classification filtered reads
- search all contig sequences using NCBI nucleotide [BLAST] against a database of your choice (we recommend the version 5 NCBI nt DB)


**NOTE:** You will need to create/download databases for [Kraken2], [Centrifuge] and [BLAST] in order to get the most out of this workflow!


### Pre-requisites

#### Taxonomic Classification for [Kraken2] and [Centrifuge]

For taxonomic classification with [Kraken2] and [Centrifuge], you will need to download (or build) databases for these programs so that you may use them within the `nf-villumina` workflow. 

You can point to the Kraken2 and Centrifuge database with `export KRAKEN2_DB=/path/to/kraken2/database` and `export CENTRIFUGE_DB=/path/to/centrifuge/database/prefix` in your `~/.bashrc` so you don't need to specify it each time you run the workflow with `--kraken2_db /path/to/kraken2/standard2 --centrifuge_db /path/to/centrifuge/nt-2018-03-03/nt`

#### Kraken2 DBs

- [MiniKraken2_v2_8GB](ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz): (5.5GB) 8GB Kraken 2 Database built from the Refseq bacteria, archaea, and viral libraries and the GRCh38 human genome
- [GTDB_r89_54k Kraken2 DBs](https://monash.figshare.com/articles/GTDB_r89_54k/8956970): There are multiple Kraken2 DBs of various sizes available for download. For more info, see https://github.com/rrwick/Metagenomics-Index-Correction and the manuscript: Méric, Wick et al. (2019) Correcting index databases improves metagenomic studies. doi: https://doi.org/10.1101/712166 


#### Centrifuge DBs

- [NCBI nucleotide non-redundant sequences (2018-03-03) (64 GB)](ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/nt_2018_3_3.tar.gz)
- [GTDB_r89_54k Centrifuge DB (108 GB tar file)](https://monash.figshare.com/ndownloader/files/16378439): For more info, see https://github.com/rrwick/Metagenomics-Index-Correction and the manuscript: Méric, Wick et al. (2019) Correcting index databases improves metagenomic studies. doi: https://doi.org/10.1101/712166 


### [BLAST] DBs

For nf-villumina, you must have a version 5 BLAST DB with embedded taxonomic information installed, e.g. version 5 `nt` DB (see https://ftp.ncbi.nlm.nih.gov/blast/db/v5/)

You can download pre-built [BLAST] DBs like `nt` and `nr` from [the NCBI FTP site](https://ftp.ncbi.nlm.nih.gov/blast/db/) using the `update_blastdb.pl` script included with your install of BLAST+ to download and/or update your local BLAST databases.

Show all available databases:

```bash
$ update_blastdb.pl --showall
```

Download the BLASTDB version 5 "nt" database to your current directory decompressing files and deleting original compressed archives:

```bash
update_blastdb.pl --blastdb_version 5 nt --decompress
``` 

**NOTE:** For ease of use, all databases should be downloaded to the same directory (e.g. `/opt/DB/blast` set in `$BLASTDB` environment variable in your `~/.bashrc`)


Check that your database has been downloaded properly and has taxids associated with the sequences contained within it:

```bash
$ blastdbcheck -db nt -must_have_taxids -verbosity 3
```


### Documentation
The peterk87/nf-villumina pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->



### Credits
peterk87/nf-villumina was originally written by Peter Kruczkiewicz.

Bootstrapped with [nf-core/tools](https://github.com/nf-core/tools) `nf-core create`. 

Thank you to the [nf-core/tools](https://github.com/nf-core/tools) team for a great tool for bootstrapping creation of a production ready Nextflow workflows.


[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[Singularity]: https://sylabs.io/
[Kraken2]: https://ccb.jhu.edu/software/kraken2/
[Centrifuge]: https://ccb.jhu.edu/software/centrifuge/manual.shtml
[fastp]: https://github.com/OpenGene/fastp
