# *lsa*BGC
### Lineage Specific Analysis of Biosynthetic Gene Clusters

*lsa*BGC is a modular, software suite designed to provide a comprehensive set of functions for investigating and mining for 
biosynthetic gene cluster diversity across a focal lineage/taxa of interest using AntiSMASH based annotation. It consists of 
8 independent programs: `lsaBGC-Ready.py`, `lsaBGC-Cluster.py`, `lsaBGC-Expansion`, `lsaBGC-Refine.py`, `lsaBGC-See.py`, 
`lsaBGC-PopGene.py`, `lsaBGC-Divergence.py`, and `lsaBGC-DiscoVary.py`.

![](https://github.com/Kalan-Lab/lsaBGC/blob/main/docs/images/lsaBGC1.1_Simplified.png)

## Major Updates 

* Jun 09, 2022 - Fixed issues with `lsaBGC-Ready.py` & New Tutorial check it out [here](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Tutorial:-Exploring-BGCs-in-Cutibacterium)!
* Jun 06, 2022 - Major updates to `lsaBGC-Ready.py` - the new recommended program for setting-up to run the lsaBGC suite.
* May 24, 2022 - `lsaBGC-Ready.py` is now available and can take pre-computed antiSMASH BGC predictions, along with optional BiG-SCAPE clustering results, to produce the required inputs for major lsaBGC analytical programs (`lsaBGC-See.py`, `lsaBGC-Refine.py`, `lsaBGC-PopGene.py`, `lsaBGC-DiscoVary.py`). 

## Documentation:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background)
2. [Extended Installation Guide (see below for more concise version of installation)](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation)
3. [The Object Oriented Core of lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/02.-The-Object-Oriented-Core-of-lsaBGC)
4. [Tutorial: Exploring BGCs in Cutibacterium](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Tutorial:-Exploring-BGCs-in-Cutibacterium)
5. [Generating Required Inputs for lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/04.-Generating-Required-Inputs-for-lsaBGC)
6. [Clustering BGCs into GCFs](https://github.com/Kalan-Lab/lsaBGC/wiki/05.-Clustering-BGCs-into-GCFs)
7. [Refinement of BGCs Belonging to GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/06.-Refinement-of-BGCs-Belonging--to-GCF)
8. [Visualizing GCFs Across Phylogenies](https://github.com/Kalan-Lab/lsaBGC/wiki/07.-Visualizing-GCFs-Across-Phylogenies)
9. [High throughput Detection of New GCF Instances Across Draft Genome Assemblies](https://github.com/Kalan-Lab/lsaBGC/wiki/08.-High-throughput-Detection-of-New-GCF-Instances-Across-Draft-Genome-Assemblies)
10. [Assessing Evolutionary Linkage of BGCs with their Genome wide Contexts](https://github.com/Kalan-Lab/lsaBGC/wiki/09.-Assessing-Evolutionary-Linkage-of-BGCs-with-their-Genome-wide-Contexts)
11. [Population Genetics Analysis of Genes Found in a GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/10.-Population-Genetics-Analysis-of-Genes-Found-in-a-GCF)
12. [Discovering Novel Variations in GCF Genes from Raw Sequencing Reads]()
13. [Benchmarking Gene Detection through Expansion vs. DiscoVary](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-Benchmarking-Gene-Detection-through-Expansion-vs.-DiscoVary)
14. [The lsaBGC AutoAnalyze Workflow](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-The-lsaBGC-AutoAnalyze-Workflow)
15. [Running test datasets for core lsaBGC programs](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases)

*Documentation moving to "Read the Docs" soon!*

## Installation:

Should take ~10 minutes.

#### Note, through these steps, please make sure to change the dummy paths `/path/to/conda_env/` and `/path/to/lsaBGC` with the desired location of the conda environment on your computing system and the location of where you clone the lsaBGC git repository on your system, respectively!!!

To install, please take the following steps:

1. Clone this git repository:

```git clone git@github.com:Kalan-Lab/lsaBGC.git```

2. Setup the conda environment using the yml file.

```
conda env create -f lsaBGC_environment.yml -p /path/to/conda_env/
```

3. Activate the environment and perform setup and pip installation in the git repository:
```
# activate the conda environment for lsaBGC just created
conda activate /path/to/conda_env/

# change directories to where the Git repo for lsaBGC was downloaded
cd /path/to/lsaBGC/

# perform python install within conda environment
python setup.py install
pip install -e .
```

4. lsaBGC uses OrthoFinder (v2.5.4) to cluster proteins in BGCs into homolog groups, 
this process can require lots of files to be written and these limits are often
controlled by certain settings on servers. In my experience, the soft limit is often 
the problem. For more insight into such constraints see: 
https://github.com/davidemms/OrthoFinder/issues/384

While other ortholog grouping software are available, OrthoFinder2 offers several
benefits to ensure the most high quality ortholog grouping.

To automatically set your soft limit to be 1 million files everytime you 
load the conda environment for lsaBGC, please run the following commands (again
make sure to replace the dummy paths!):
```
mkdir -p /path/to/conda_env/etc/conda/activate.d
touch /path/to/conda_env/etc/conda/activate.d/env_vars.sh
echo $'#!/bin/sh\n\ulimit -n 1000000\n' > /path/to/conda_env/etc/conda/activate.d/env_vars.sh
```

We also recommend setting the environment variable $TMPDIR to a larger space than the typical default `/tmp/` which is usually short on space. This becomes needed for large `sort` operations, which can pop up when using CompareM.

```
echo $'export $TMPDIR=/path/to/larger_tmp_dir/' >> /path/to/conda_env/etc/conda/activate.d/env_vars.sh
```

5. Setup database(s) for annotation used by `lsaBGC-Ready.py`. This is currently just the,
KOfam profile HMMs (~5GB). To setup databases, simply run the script:

```
setup_annotation_dbs.py 
```

## Acknowledgements:

We would like to thank members of the Kalan lab, Currie lab, Kwan lab, and Anantharaman lab at UW Madison for feedback on the development of lsaBGC.

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
