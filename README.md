# *lsa*BGC
### Lineage Specific Analysis (*lsa*) of Biosynthetic Gene Clusters (BGC)

[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400)](https://www.biorxiv.org/content/10.1101/2022.04.20.488953v2)
[![Documentation](https://img.shields.io/badge/documentation-wiki-darkgreen?style=flat-square&maxAge=2678400)](https://github.com/Kalan-Lab/lsaBGC/wiki)


*lsa*BGC offers modular programs, as well as workflows, designed for investigating and mining for biosynthetic gene cluster diversity across a focal lineage/taxa of interest. 

Compatible with antiSMASH, GECCO, and DeepBGC. 

<image src="https://github.com/Kalan-Lab/lsaBGC/blob/main/docs/images/lsaBGC1.1_Simplified.png">

## Documentation and How to Get Started:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background-&-Considerations)
2. [An Overview of Final Results from lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-Overview-of-lsaBGC-AutoAnalyze's-Final-Results)
3. [Quick Start - Using the simple lsaBGC-Easy.py (bacterial) and lsaBGC-Euk-Easy.py (fungal) workflows](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial)
4. [Modular Usage - Exploring BGCs in Cutibacterium](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium)
5. [GSeeF - quick and simple visualization of GCFs/BGCs across a species phylogeny](https://github.com/Kalan-Lab/lsaBGC/wiki/17.-GSeeF---Visualizing-GCF-Cluster-Presence-and-Annotation-Along-a-Species-Phylogeny)
6. [visualize_BGC-ome - quick and simple visualization of a sample's BGC-ome](https://github.com/Kalan-Lab/lsaBGC/wiki/19.-Plot-Sample-BGC-ome)
7. [***new***: Investigate a single cluster of related BGCs using the sibling suite *zol*](https://github.com/Kalan-Lab/zol)

## Installation:

Installation can be performed via conda and should take ~5-10 minutes and has been tested on both unix (specifically Ubuntu) and macOS. We are happy to attempt to address issues with installation if any arise, please open a Git Issues case:

```
# 1. clone Git repo and cd into it!
git clone https://github.com/Kalan-Lab/lsaBGC
cd lsaBGC/

# 2. create conda environment using yaml file and activate it!
conda env create -f lsaBGC_env.yml -p /path/to/lsaBGC_conda_env/
conda activate /path/to/lsaBGC_conda_env/

# 3. complete python installation with the following commands:
python setup.py install
pip install -e .
```

Optional, but recommended, command to download KOfams + PGAP HMMs + MIBiG protein FASTA for annotation:

```
# Warning: can take >10 minutes! 
# Can skip to run tests first to make sure things are working properly.
# within lsaBGC Git repo with conda environment activated:
setup_annotation_dbs.py
```

If clustering of BGCs into GCFs using BiG-SCAPE is preferred to lsaBGC-Cluster.py, setup BiG-SCAPE using the following:

```
setup_bigscape.py
```

Additional, information pertaining to installation can be found at: [Installation Guide](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation)

A small test case is provided here and can be run after installation by simply issuing (takes around ~7 minutes using 4 cpus/threads):

```
# Warning: uses 4 cpus/threads! 
bash run_tests.sh
```

There are also additional [test cases](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases) to demonstrate usage of individual programs along with expected outputs from commands. We also have a [walk-through tutorial Wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium) to showcase the use of the suite and relations between core programs.

The major outputs of the final `lsaBGC-AutoAnalyze.py` run are in the resulting folder `test_case/lsaBGC_AutoAnalyze_Results/Final_Results/` and described on [this wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-Overview-of-lsaBGC-AutoAnalyze's-Final-Results). Examples for the final AutoAnalyze results from an `lsaBGC-Easy.py` run on Cutibacterium avidum can be found [here on Google Drive](https://drive.google.com/drive/u/1/folders/1jHFFOUTd4SbIO-xiGG8MWTZaP1U4RF1j). 

## Quick Start - using `lsaBGC-Easy.py` and `lsaBGC-Euk-Easy.py` 

Check out how to use `lsaBGC-Easy.py` and `lsaBGC-Euk-Easy.py` on [their wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial)!

![image](https://user-images.githubusercontent.com/4260723/181613839-df183cdc-1103-403f-b5d1-889484f52be9.png)

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
