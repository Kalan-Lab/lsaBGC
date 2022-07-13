# *lsa*BGC
### Lineage Specific Analysis of Biosynthetic Gene Clusters

*lsa*BGC is a modular software suite designed to provide a comprehensive set of functions for investigating and mining for 
biosynthetic gene cluster diversity across a focal lineage/taxa of interest (currently limited to bacterial) using AntiSMASH 
based annotation. It consists of 8 independent programs: `lsaBGC-Ready.py`, `lsaBGC-Cluster.py`, `lsaBGC-Expansion`, 
`lsaBGC-Refine.py`, `lsaBGC-See.py`, `lsaBGC-PopGene.py`, `lsaBGC-Divergence.py`, and `lsaBGC-DiscoVary.py`.

![](https://github.com/Kalan-Lab/lsaBGC/blob/main/docs/images/lsaBGC1.1_Simplified.png)

## Major Updates 

* Jul 12, 2022 - Huge thanks to Martin Larralde for recommendation to update GECCO processing and better identify "protocore"-esque genes in BGCs using CRF weights instead of e-values! Also, have now added a simplified sub-section down below with [notes on the PopGene report table](#notes-on-the-popgene-table-report).
* Jul 10, 2022 - Several updates made. Fixed small issues with smooth running of new framework, `lsaBGC-Ready.py`. Removed some dependencies and have added GToTree for creating species phylogeny + estimated sample to sample amino acid expected divergences. New small test dataset now included in this repo for immediate testing + much simplified installation guide. Most major change is that lsaBGC now works with DeepBGC and GECCO predictions! lsaBGC's backend relies on 'proto-core homolog groups' / 'rule-based key domains' determined by AntiSMASH, to get around the absence of such marker genes/domains in DeepBGC and GECCO predictions, domains in the highest 10% of deebgc_scores or lowest 10% of e-values are treated as "proto-core" and used in `lsaBGC-AutoExpansion.py`/`lsaBGC-DiscoVary.py` as well as highlighted/treated as the "core" in `lsaBGC-PopGene.py` reports.
* Jun 26, 2022 - Added "loose" mode to `lsaBGC-Expansion.py` and option for users to manually define "protocore" homolog groups. Also, "protocore" homolog groups for a GCF now must have "rule-based" marker to exclude MGEs like transposons which insert within protocore regions of BGCs.
* Jun 19, 2022 - Have set MAGUS as the default protein alignment method (highly scalable wrapper of mafft) + updated notes on scalability.
* Jun 18, 2022 - Updated [`lsaBGC-AutoAnalyze.py`](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-The-lsaBGC-AutoAnalyze-Workflow) (automated lsaBGC analysis for each GCF) for better integration into new framework based around `lsaBGC-Ready.py`. 
* Jun 14, 2022 - Added [note on scalability](#user-content-notes-on-scalability), below on this page, and future plans to address them.
* Jun 09, 2022 - Fixed issues with `lsaBGC-Ready.py` & New Tutorial check it out [here](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium)!
* Jun 06, 2022 - Major updates to `lsaBGC-Ready.py` - the new recommended program for setting-up to run the lsaBGC suite.
* May 24, 2022 - `lsaBGC-Ready.py` is now available and can take pre-computed antiSMASH BGC predictions, along with optional BiG-SCAPE clustering results, to produce the required inputs for major lsaBGC analytical programs (`lsaBGC-See.py`, `lsaBGC-Refine.py`, `lsaBGC-PopGene.py`, `lsaBGC-DiscoVary.py`). 

## Documentation:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background)
2. [The Object Oriented Core of lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/02.-The-Object-Oriented-Core-of-lsaBGC)
3. [Quick Start & Tutorial: Exploring BGCs in Cutibacterium](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium)
4. [Generating Required Inputs for lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/04.-Generating-Required-Inputs-for-lsaBGC)
5. [Clustering BGCs into GCFs](https://github.com/Kalan-Lab/lsaBGC/wiki/05.-Clustering-BGCs-into-GCFs)
6. [Refinement of BGCs Belonging to GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/06.-Refinement-of-BGCs-Belonging--to-GCF)
7. [Visualizing GCFs Across Phylogenies](https://github.com/Kalan-Lab/lsaBGC/wiki/07.-Visualizing-GCFs-Across-Phylogenies)
8. [High throughput Detection of New GCF Instances Across Draft Genome Assemblies](https://github.com/Kalan-Lab/lsaBGC/wiki/08.-High-throughput-Detection-of-New-GCF-Instances-Across-Draft-Genome-Assemblies)
9. [Assessing Evolutionary Linkage of BGCs with their Genome wide Contexts](https://github.com/Kalan-Lab/lsaBGC/wiki/09.-Assessing-Evolutionary-Linkage-of-BGCs-with-their-Genome-wide-Contexts)
10. [Population Genetics Analysis of Genes Found in a GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/10.-Population-Genetics-Analysis-of-Genes-Found-in-a-GCF)
11. [Discovering Novel Variations in GCF Genes from Raw Sequencing Reads]()
12. [Benchmarking Gene Detection through Expansion vs. DiscoVary](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-Benchmarking-Gene-Detection-through-Expansion-vs.-DiscoVary)
13. [The lsaBGC AutoAnalyze Workflow](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-The-lsaBGC-AutoAnalyze-Workflow)
14. [Running test datasets for core lsaBGC programs](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases)

*Documentation moving to "Read the Docs" soon!*

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

Optional, but recommended, command to download KOfams (+ other databases in the near future):

```
# Warning: can take ~5-10 minutes!
# within lsaBGC Git repo with conda environment activated:
setup_annotation_dbs.py
```

Additional, information pertaining to installation can be found at: [Installation Guide](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation)

A small test case is provided here and can be run after installation by simply issuing (takes around ~10 minutes using 8 cores):

```
bash run_tests.sh
```

Additionally we suggest checking out additional [test cases](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases) to demonstrate usage of individual programs along with expected outputs from commands. We also have a [quick start + walk-through tutorial Wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium) to showcase the use of the suite and relations between core programs.

The major outputs of the final `lsaBGC-AutoAnalyze.py` run in the resulting folder `test_case/lsaBGC_AutoAnalyze_Results/` include:

* **GCF_Homolog_Group_Information.txt** - a table with each row corresponding to a homolog group for each GCF - gives a quick high-level overview of the pan-secondary metabolome of your microbe of interest that you can filter and scroll through! See below section "Notes on the "lsaBGC-PopGene.py" Table Report" for more information. 
* **Consensus_Sequence_Similarity_of_Homolog_Groups.pdf** - a PDF showing the phylogenetic distribution of GCFs found in >25% of samples. Coloring represents difference to consensus gene sequences across all samples.
* **GCF_Conservation_and_PopStats_Views.pdf** - 2 pg PDF report showing gene schemas. Degree of transparency = conservation of gene across BGCs. Coloring corresponds to either Tajima's D or (will probably be getting updated soon to be)
* **GCF_Divergence.pdf** - a boxplot of the Beta-RD values (the divergence of BGCs between pairs of genomes to the expected divergence, i.e. genome-wide ANI/AAI or divergence of single-copy core genes) for each GCF.
* **See/GCF_X/BGC_Visualization.species_phylogeny.pdf** - for each GCF_X in the `See/` sub-directory you can see a phylogenetic distribution of BGCs belonging to the GCF accounting for BGC splitting within a sample's genome due to assembly fragmentation. These are BGCs detected in the "expansion/additional genomes" set by `lsaBGC-Expansion.py`. 

## Notes on the PopGene Report Table

If you run `lsaBGC-PopGene.py` for an individual GCF or `lsaBGC-AutoAnalyze.py` for all the GCFs one of the major outputs from the suite will be a table with each row corresponding to a homolog group from each GCF.

This report is automatically sorted by the consensus order of homolog groups allowing for easy intuitive scrolling. For full information on each column, check out the [data dictionary table](https://github.com/Kalan-Lab/lsaBGC/wiki/10.-Population-Genetics-Analysis-of-Genes-Found-in-a-GCF#data-dictionary) on the `lsaBGC-PopGene.py` wiki.

For more advanced evolutionary statistics, I recommend primarily consulting the Tajima's D and Median Beta-RD statistics, ***dN/dS & Fst are experimental implementations!!!*** and often times not very informative. 

**Median Beta-RD** is the median value of pairwise calculations of gene divergence to expected divergence for a pair of genomes based on genome-wide values such as ANI/AAI or divergence of single-copy core genes (e.g. ribosomal proteins).

**Tajima's D** can be most intiutively thought of as the proportion of high-frequency to low-frequency/rare minor allele positions along the multiple sequence alignment of a gene. It is thus essential to account for representation bias in your dataset and as such we perform de-replication (check out [drep](https://github.com/MrOlm/drep)) of genomes prior to `lsaBGC-PopGene.py` analysis to properly interpret the statistics.

In the near future, I plan to incorporate a phylogeny-aware calculation of dN/dS, as well as update other evolutionary / pop-gen stats. Will announce these in major updates as they are incorporated.

A highlight of some key columns in the table after loading into Google sheets:

![image](https://user-images.githubusercontent.com/4260723/178626400-9563427f-d98d-48d6-8b4b-0e9673c6d0cd.png)

## Notes on Scalability:

Updated 06/19/2022

lsaBGC strives for large-scalability and is designed to be high-throughput. We previously ran it on >15,000 Staphylococci. However, it is recommended that the number of "Primary" genomes not exceed 300 (or 500 at most!) because de novo ortholog grouping is performed for this set (an all vs. all procedure). It is therefore recommended that users first select a distributed set of representative, "primary" genomes for their taxa of interest (which can be done via genome dereplication, check out [drep](https://github.com/MrOlm/drep)!). Then treat the rest of the genomes as "additional" (this can be 1000s of genomes!).

For evolutionary statistics calculations, codon alignments are now built using MAGUS (a divide and conquer wrapper of MAFFT) by default! 

## Future Updates Planned and of High Priority:

* Get `lsaBGC-Ready.py` working for fungi + plants.
* Incorporate additional/update evolutionary statistics in `lsaBGC-PopGene.py`! 

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
