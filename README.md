# *lsa*BGC
### Lineage Specific Analysis of Biosynthetic Gene Clusters
#### Rauf Salamzade - Kalan Lab, MMI UW-Madison

*lsa*BGC is a software suite designed to provide a comprehensive set of functions for investigating and mining for 
biosynthetic gene cluster diversity across a focal lineage/taxa of interest using AntiSMASH based annotation.

## Documentation

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

## Installation

Should take < 10 minutes.

To install, please take the following steps:

1. Clone this git repository:

```git clone git@github.com:Kalan-Lab/lsaBGC.git```

2. Setup the conda environment using the yml file.

```conda env create -f lsaBGC_environment.yml -p /path/to/conda_environment/```

3. Activate the environment and perform setup and pip installation in the git repository:
```bash
# activate the conda environment for lsaBGC just created
conda activate /path/to/conda_environment/

# change directories to where the Git repo for lsaBGC was downloaded
cd /path/to/lsaBGC/

# perform python install within conda environment
python setup.py install
pip install .
```

## Dependencies
As described in the Installation section above, dependencies can be set up easily through the use of a Conda environment and the provided yaml file.

The set of dependencies for the core lsaBGC programs and auxiliary scripts, along with versions used for testing are listed on the wiki in the Installation page.

lsaBGC was developed and tested on UNIX systems; however, there are no apprent reasons users would have difficulty running on OS X or Windows.

## Acknowledgements

We would like to thank members of the Kalan lab, Currie lab, Kwan lab, and Anantharaman lab at UW Madison for feedback on the development of lsaBGC.

## License

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
