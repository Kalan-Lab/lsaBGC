# *lsa*BGC
### Lineage Specific Analysis of Biosynthetic Gene Clusters
#### Rauf Salamzade - Kalan Lab, MMI UW-Madison

*lsa*BGC is a software suite designed to provide a comprehensive set of functions for investigating and mining for 
biosynthetic gene cluster diversity across a focal lineage/taxa of interest using AntiSMASH based annotation.

**_Note - Jun 01, 2022 While everything has passed testing cases (using the example Corynebacterium datasets linked in the repo), we are actively updating the early steps of the framework to be easier to use, require fewer dependency and integrate better with prior clusteirng results from programs such as BiG-SCAPE and BiG-SLICE. We hope to have updated documentation detailing recommendations for setting up lsaBGC analyses and an improved `lsaBGC-Ready.py` version in the next week._** 
## Major Updates:

* May 24, 2022 - `lsaBGC-Ready.py` is now available and can take pre-computed antiSMASH BGC predictions, along with optional BiG-SCAPE clustering results, to produce the required inputs for major lsaBGC analytical programs (`lsaBGC-See.py`, `lsaBGC-Refine.py`, `lsaBGC-PopGene.py`, `lsaBGC-DiscoVary.py`). 

## Documentation:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background)
2. [Extended Installation Guide (see below for more concise version of installation)](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation)
3. [The Object Oriented Core of lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/02.-The-Object-Oriented-Core-of-lsaBGC)
4. [Detailed Walkthrough](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Detailed-Walkthrough)
5. [Generating Required Inputs for lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/04.-Generating-Required-Inputs-for-lsaBGC)
6. [Clustering BGCs into GCFs](https://github.com/Kalan-Lab/lsaBGC/wiki/05.-Clustering-BGCs-into-GCFs)
7. [Refinement of BGCs Belonging to GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/06.-Refinement-of-BGCs-Belonging--to-GCF)
8. [Visualizing GCFs Across Phylogenies](https://github.com/Kalan-Lab/lsaBGC/wiki/07.-Visualizing-GCFs-Across-Phylogenies)
9. [High throughput Detection of New GCF Instances Across Draft Genome Assemblies](https://github.com/Kalan-Lab/lsaBGC/wiki/08.-High-throughput-Detection-of-New-GCF-Instances-Across-Draft-Genome-Assemblies)
10. [Assessing Evolutionary Linkage of BGCs with their Genome wide Contexts](https://github.com/Kalan-Lab/lsaBGC/wiki/09.-Assessing-Evolutionary-Linkage-of-BGCs-with-their-Genome-wide-Contexts)
11. [The lsaBGC AutoAnalyze Workflow](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-The-lsaBGC-AutoAnalyze-Workflow)
12. [Benchmarking Gene Detection through Expansion vs. DiscoVary](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-Benchmarking-Gene-Detection-through-Expansion-vs.-DiscoVary)
13. [Running test datasets for core lsaBGC programs](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases)

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

```
# To automatically set your soft limit to be 1 million files everytime you 
# load the conda environment for , please run the following commands (again
# make sure to replace the dummy paths!):
mkdir -p /path/to/conda_env/etc/conda/activate.d
touch /path/to/conda_env/etc/conda/activate.d/env_vars.sh
echo $'#!/bin/sh\n\ulimit -n 1000000' > /path/to/conda_env/etc/conda/activate.d/env_vars.sh
```

## Dependencies:
As described in the Installation section above, dependencies can be set up easily through the use of a Conda environment and the provided yaml file.

The set of dependencies for the core lsaBGC programs and auxiliary scripts, along with versions used for testing are listed on the wiki in the Installation page.

lsaBGC was developed and tested on UNIX systems; however, there are no apprent reasons users would have difficulty running on OS X or Windows.

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
