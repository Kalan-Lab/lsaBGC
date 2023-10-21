# This code is adapted from BiG-SCAPE's Dockerfile - if you are using zol - you should definitely check out
# BiG-SCAPE/CORASON incase it suites your needs better (https://www.nature.com/articles/s41589-019-0400-9) - eg.
# you are interested in clustering diverse gene-clusters based on protein domain simlarity or are interested
# in investigating the variability of contexts for a single reference gene!

FROM continuumio/miniconda3
LABEL maintainer="Rauf Salamzade - Kalan Lab, UW-Madison"

WORKDIR /usr/src
SHELL ["/bin/bash", "-c"]


# Clone lsaBGC github repo and create conda environment, then create and activate conda environment,
# and install lsaBGC
RUN apt-get update && apt-get install -y git wget libtiff5 gcc && \
    git clone -b develop-v1.40 https://github.com/Kalan-Lab/lsaBGC && rm -rf lsaBGC/test_case.tar.gz && \
    conda install -n base conda-libmamba-solver && conda config --set solver libmamba && \
    conda env create -f /usr/src/lsaBGC/lsaBGC_env.yml -p /usr/src/lsaBGC_conda_env/ && \
    conda create -p /usr/src/antismash_conda_env/ -c bioconda -c conda-forge -c defaults python=3.11 hmmer2 hmmer diamond fasttree prodigal blast muscle glimmerhmm cxx-compiler meme && \
    source activate /usr/src/antismash_conda_env/ && \
    wget https://dl.secondarymetabolites.org/releases/7.0.1/antismash-7.0.1.tar.gz && tar -zxf antismash-7.0.1.tar.gz && \
    pip install ./antismash-7.0.1 && download-antismash-databases && \
    download-antismash-databases && \
    conda remove --force meme && conda deactivate && \
    conda clean --all -y && \
    echo "source activate /usr/src/lsaBGC_conda_env/" > ~/.bashrc && source ~/.bashrc && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Install lsaBGC
WORKDIR /usr/src/lsaBGC/
ENV PATH /usr/src/lsaBGC_conda_env/bin:$PATH
RUN python setup.py install && pip install -e . && setup_annotation_dbs.py -nk -dsh && setup_bigscape.py && \
    chmod -R 777 /usr/src/lsaBGC/ && chmod 777 /home 
USER 1000:1000
RUN mkdir /home/input /home/output
WORKDIR /home
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENTRYPOINT ["LSABGC"]
CMD ["--help"]


