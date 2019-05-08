Bootstrap:docker
From:nfcore/base

%labels
    MAINTAINER Peter Kruczkiewicz
    DESCRIPTION Singularity image containing all requirements for the peterk87/nf-villumina pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-villumina-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
