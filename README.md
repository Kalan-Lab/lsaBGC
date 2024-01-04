# *lsa*BGC
### Lineage Specific Analysis (*lsa*) of Biosynthetic Gene Clusters (BGC)

[![Manuscript](https://img.shields.io/badge/Manuscript-MGen-darkblue?style=flat-square&maxAge=2678400)](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988)
[![Documentation](https://img.shields.io/badge/Documentation-Wiki-darkgreen?style=flat-square&maxAge=2678400)](https://github.com/Kalan-Lab/lsaBGC/wiki)
[![Docker](https://img.shields.io/badge/Docker-DockerHub-darkred?style=flat-square&maxAge=2678400)](https://hub.docker.com/r/raufs/lsabgc)

*lsa*BGC offers modular programs, as well as workflows, designed for investigating and mining for biosynthetic gene cluster diversity across a focal lineage/taxa of interest. 

Compatible with BGC predictions from antiSMASH, GECCO, and DeepBGC. 

![image](https://github.com/Kalan-Lab/lsaBGC/assets/4260723/aa0703e5-299b-4f1c-820a-f4ebae2b264d)

## Documentation and How to Get Started:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background-&-Considerations)
2. [An Overview of Final Results from lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-Overview-of-lsaBGC-AutoAnalyze's-Final-Results)
3. [Quick Start - Using the simple lsaBGC-Easy.py (bacterial) and lsaBGC-Euk-Easy.py (fungal) workflows](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial)
      - *Use on any platform in just 3 steps [via Docker](https://github.com/Kalan-Lab/lsaBGC#using-docker-for-major-workflows-only)* 
4. [Modular Usage - Exploring BGCs in Cutibacterium](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium)
5. [GSeeF - quick and simple visualization of GCFs/BGCs across a species phylogeny](https://github.com/Kalan-Lab/lsaBGC/wiki/17.-GSeeF---Visualizing-GCF-Cluster-Presence-and-Annotation-Along-a-Species-Phylogeny)
6. [visualize_BGC-ome - quick and simple visualization of a sample's BGC-ome](https://github.com/Kalan-Lab/lsaBGC/wiki/19.-Plot-Sample-BGC-ome)
7. [***new***: Investigate a single cluster of related BGCs using the sibling suite *zol*](https://github.com/Kalan-Lab/zol)

## IMPORTANT: PLEASE USE v1.52+:

Please make sure to use v1.52+ of the pipeline - if you are using `lsaBGC-Easy.py` with antiSMASH, the default settings for antiSMASH based BGC prediction from v1.38 to v1.51 included the argument `--taxon fungi` by mistake. It should only be the default for the analagous `lsaBGC-Euk-Easy.py` program. 

## Installation:

### Using Conda (for full usage of suite)

Installation can be performed via conda (see below for Docker) and should take ~5 minutes with mamba or ~10-20 minutes with conda and has been tested on both unix (specifically Ubuntu) and macOS. We are happy to attempt to address issues with installation if any arise, please open a Git Issues case:

```bash
# 1. clone Git repo and cd into it!
git clone https://github.com/Kalan-Lab/lsaBGC
cd lsaBGC/

# 2. create conda environment using yaml file and activate it!
# For a much faster installation replace "conda" in the following
# commands with "mamba" (after installing mamba in your base conda
# environment)
mamba env create -f lsaBGC_env.yml -p /path/to/lsaBGC_conda_env/
conda activate /path/to/lsaBGC_conda_env/

# 3. complete python installation with the following commands:
# since version 1.50, setup.py sonicparanoid no longer included
# via pip based installation - to use SonicParanoid, please use
# Docker - might incorporate as a conda installation in the future.
python setup.py install
pip install -e .
```

Optional, but recommended, command to download KOfams + PGAP HMMs + MIBiG protein FASTA for annotation:

```bash
# Warning: can take >10 minutes! 
# Can skip to run tests first to make sure things are working properly.
# within lsaBGC Git repo with conda environment activated:
setup_annotation_dbs.py
```

If clustering of BGCs into GCFs using BiG-SCAPE is preferred to lsaBGC-Cluster.py, setup BiG-SCAPE using the following:

```
setup_bigscape.py
```
   
A small test case is provided here and can be run after installation by simply issuing (takes around ~7 minutes using 4 cpus/threads):

```bash
# Warning: uses 4 cpus/threads! 
bash run_tests.sh
```

There are also additional [test cases](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases) to demonstrate usage of individual programs along with expected outputs from commands. We also have a [walk-through tutorial Wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium) to showcase the use of the suite and relations between core programs.

The major outputs of the final `lsaBGC-AutoAnalyze.py` run are in the resulting folder `test_case/lsaBGC_AutoAnalyze_Results/Final_Results/` and described on [this wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-Overview-of-lsaBGC-AutoAnalyze's-Final-Results). Examples for the final AutoAnalyze results from an `lsaBGC-Easy.py` run on Cutibacterium avidum can be found [here on Google Drive](https://drive.google.com/drive/u/1/folders/1jHFFOUTd4SbIO-xiGG8MWTZaP1U4RF1j). 
   
### Using Docker (for major workflows only)

A docker image is provided for the `lsaBGC-Easy.py` and `lsaBGC-Euk-Easy.py` workflows together with a wrapper script. The image is pretty large (~25Gb without SonicParanoid, ~33Gb with SonicParanoid) but includes all the databases and dependencies needed for lsaBGC, BiG-SCAPE, antiSMASH, and GECCO analysis. For lsaBGC, to save space, the KOfam database is not included. For antiSMASH, MEME is not incldued, thus RODEO and CASSIS analyses are not available.
   
To use the latest Docker image, please: (1) install Docker and (2) download the wrapper script:

```bash
# 1. download wrapper script for running image with SonicParanoid
wget https://raw.githubusercontent.com/Kalan-Lab/lsaBGC/main/docker/withSonicParanoid/run_LSABGC.sh

# or

# choose version with script for running image without SonicParanoid
wget https://raw.githubusercontent.com/Kalan-Lab/lsaBGC/main/docker/withoutSonicParanoid/run_LSABGC.sh

# 2. run it
bash run_LSABGC.sh
```
   
## Quick Start - using `lsaBGC-Easy.py` and `lsaBGC-Euk-Easy.py` 

Check out how to use `lsaBGC-Easy.py` and `lsaBGC-Euk-Easy.py` on [their wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial)!

![image](https://user-images.githubusercontent.com/4260723/181613839-df183cdc-1103-403f-b5d1-889484f52be9.png)
***
![image](https://github.com/Kalan-Lab/lsaBGC/assets/4260723/d35875e5-4299-4c02-a874-50bbf8f1cf30)
   
## Acknowledgements:

We would like to thank members of the Kalan lab, Currie lab, Kwan lab, Anantharaman, and Pepperell labs at UW Madison for feedback on the development of lsaBGC.

## Feedback:

Issues or suggestions for new features / changes to approaches? Please create an issue/ticket on GitHub issues and let us know!

## License:

```
BSD 3-Clause License

Copyright (c) 2021, Kalan-Lab
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
