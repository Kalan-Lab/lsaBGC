# This code is adapted from BiG-SCAPE's Dockerfile - if you are using zol - you should definitely check out
# BiG-SCAPE/CORASON incase it suites your needs better (https://www.nature.com/articles/s41589-019-0400-9) - eg.
# you are interested in clustering diverse gene-clusters based on protein domain simlarity or are interested
# in investigating the variability of contexts for a single reference gene!

FROM continuumio/miniconda3
LABEL maintainer="Rauf Salamzade - Kalan Lab, UW-Madison"

WORKDIR /usr/src
SHELL ["/bin/bash", "-c"]
RUN apt-get update && apt-get install -y git wget

# Clone lsaBGC github repo and create conda environment
RUN git clone --branch develop https://github.com/Kalan-Lab/lsaBGC && rm -rf lsaBGC/test_case.tar.gz

# Create and activate conda environment
RUN conda install -c conda-forge mamba && \
    mamba env create -f /usr/src/lsaBGC/lsaBGC_env.yml -p /usr/src/lsaBGC_conda_env/ && \
    mamba create -p /usr/src/antismash_conda_env/ -c bioconda -c conda-forge -c defaults antismash -y && \
    source activate /usr/src/antismash_conda_env/ && download-antismash-databases && conda deactivate && \
    conda clean --all -y && conda remove mamba && \
    echo "source activate /usr/src/lsaBGC_conda_env/" > ~/.bashrc && source ~/.bashrc && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Install lsaBGC
WORKDIR /usr/src/lsaBGC/
ENV PATH /usr/src/lsaBGC_conda_env/bin:$PATH
RUN python setup.py install && pip install -e . && setup_annotation_dbs.py -nk -dsh && setup_bigscape.py && \
    chmod -R 555 /usr/src/lsaBGC/ && chmod 777 /home

USER 1000:1000
RUN mkdir /home/input /home/output
WORKDIR /home
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENTRYPOINT ["lsaBGC-Easy.py"]
CMD ["--help"]