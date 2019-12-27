# peterk87/nf-villumina
**Generic viral Illumina sequence analysis pipeline**

[![Build Status](https://travis-ci.org/peterk87/nf-villumina.svg?branch=master)](https://travis-ci.org/peterk87/nf-villumina)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/peterk87/nf-villumina.svg)](https://hub.docker.com/r/peterk87/nf-villumina)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2925)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with a [Singularity][] containers making installation trivial and results highly reproducible.

`nf-villumina` will remove low quality reads ([fastp]) and remove non-viral reads using [Kraken2][] and [Centrifuge][] classification results. [Unicycler] and [Shovill] *de novo* assemblies are produced on the taxonomic classification filtered reads. NCBI nucleotide [BLAST] search is performed against a database of your choice (we recommend )


**NOTE:** If you want to perform taxonomic classification with Kraken2 and Centrifuge, you will need to download or build the Kraken2 and Centrifuge classification databases. You can point to the Kraken2 and Centrifuge database with `export KRAKEN2_DB=/path/to/kraken2/database` and `export CENTRIFUGE_DB=/path/to/centrifuge/database/prefix` in your `~/.bashrc` so you don't need to specify it each time you run the workflow with `--kraken2_db /path/to/kraken2/standard2 --centrifuge_db /path/to/centrifuge/nt-2018-03-03/nt`


### Pre-requisites

#### Taxonomic Classification for [Kraken2] and [Centrifuge]



MiniKraken2_v2_8GB: (5.5GB) 8GB Kraken 2 Database built from the Refseq bacteria, archaea, and viral libraries and the GRCh38 human genome (ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz)


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



[Singularity]: https://sylabs.io/
[Kraken2]: https://ccb.jhu.edu/software/kraken2/
[Centrifuge]: https://ccb.jhu.edu/software/centrifuge/manual.shtml
[fastp]: https://github.com/OpenGene/fastp
