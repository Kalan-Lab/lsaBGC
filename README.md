# *lsa*BGC
### Lineage Specific Analysis of Biosynthetic Gene Clusters

[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400)](https://www.biorxiv.org/content/10.1101/2022.04.20.488953v2)


*lsa*BGC is a modular software suite designed to provide a comprehensive set of functions for investigating and mining for biosynthetic gene cluster diversity across a focal lineage/taxa of interest using AntiSMASH based annotation. It consists of 8 independent programs: `lsaBGC-Ready.py`, `lsaBGC-Cluster.py`, `lsaBGC-Expansion`, `lsaBGC-Refine.py`, `lsaBGC-See.py`, `lsaBGC-PopGene.py`, `lsaBGC-Divergence.py`, and `lsaBGC-DiscoVary.py`.

![](https://github.com/Kalan-Lab/lsaBGC/blob/main/docs/images/lsaBGC1.1_Simplified.png)

## Documentation:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background-&-Considerations)
2. [An Overview of Final Results from lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-Overview-of-lsaBGC-AutoAnalyze's-Final-Results)
3. [Really Quick Start - Using lsaBGC-Easy.py](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial)
4. [Quick Start & Tutorial: Exploring BGCs in Cutibacterium](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium)
5. [Generating Required Inputs for lsaBGC](https://github.com/Kalan-Lab/lsaBGC/wiki/04.-Generating-Required-Inputs-for-lsaBGC)
6. [Clustering BGCs into GCFs](https://github.com/Kalan-Lab/lsaBGC/wiki/05.-Clustering-BGCs-into-GCFs)
7. [Refinement of BGCs Belonging to GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/06.-Refinement-of-BGCs-Belonging--to-GCF)
8. [Visualizing GCFs Across Phylogenies](https://github.com/Kalan-Lab/lsaBGC/wiki/07.-Visualizing-GCFs-Across-Phylogenies)
9. [High throughput Detection of New GCF Instances Across Draft Genome Assemblies](https://github.com/Kalan-Lab/lsaBGC/wiki/08.-High-throughput-Detection-of-New-GCF-Instances-Across-Draft-Genome-Assemblies)
10. [Assessing Evolutionary Linkage of BGCs with their Genome wide Contexts](https://github.com/Kalan-Lab/lsaBGC/wiki/09.-Assessing-Evolutionary-Linkage-of-BGCs-with-their-Genome-wide-Contexts)
11. [Population Genetics Analysis of Genes Found in a GCF](https://github.com/Kalan-Lab/lsaBGC/wiki/10.-Population-Genetics-Analysis-of-Genes-Found-in-a-GCF)
12. [Discovering Novel Variations in GCF Genes from Raw Sequencing Reads](https://github.com/Kalan-Lab/lsaBGC/wiki/11.-Discovering-Novel-Variations-in-GCF-Genes-from-Raw-Sequencing-Reads)
13. [The lsaBGC AutoAnalyze Workflow](https://github.com/Kalan-Lab/lsaBGC/wiki/12.-The-lsaBGC-AutoAnalyze-Workflow)
14. [Benchmarking Gene Detection through Expansion vs. DiscoVary](https://github.com/Kalan-Lab/lsaBGC/wiki/15.-Benchmarking-Gene-Detection-through-Expansion-vs.-DiscoVary)
15. [Running additional test datasets for core lsaBGC programs](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases)
16. [GSeeF - quick and simple visualization of BGCs across a species phylogeny](https://github.com/Kalan-Lab/lsaBGC/wiki/17.-GSeeF---Visualizing-GCF-Cluster-Presence-and-Annotation-Along-a-Species-Phylogeny)


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

Optional, but recommended, command to download KOfams + PGAP HMMs for annotation:

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

Additionally, we suggest checking out additional [test cases](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases) to demonstrate usage of individual programs along with expected outputs from commands. We also have a [quick start + walk-through tutorial Wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium) to showcase the use of the suite and relations between core programs.

The major outputs of the final `lsaBGC-AutoAnalyze.py` run are in the resulting folder `test_case/lsaBGC_AutoAnalyze_Results/Final_Results/` and described on [this wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/16.-Overview-of-lsaBGC-AutoAnalyze's-Final-Results). Examples for the final AutoAnalyze results from an `lsaBGC-Easy.py` run on Cutibacterium avidum can be found [here on Google Drive](https://drive.google.com/drive/u/1/folders/1jHFFOUTd4SbIO-xiGG8MWTZaP1U4RF1j). 

## Quick Start - using `lsaBGC-Easy.py`

Check out how to use `lsaBGC-Easy.py` on [it's wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial:-Combining-lsaBGC-with-ncbi-genome-download)!

![image](https://user-images.githubusercontent.com/4260723/181613839-df183cdc-1103-403f-b5d1-889484f52be9.png)

## Notes on the PopGene Report Table

If you run `lsaBGC-PopGene.py` for an individual GCF or `lsaBGC-AutoAnalyze.py` for all the GCFs one of the major outputs from the suite will be a table with each row corresponding to a homolog group from each GCF.

This report is automatically sorted by the consensus order of homolog groups allowing for easy intuitive scrolling. For full information on each column, check out the [data dictionary table](https://github.com/Kalan-Lab/lsaBGC/wiki/10.-Population-Genetics-Analysis-of-Genes-Found-in-a-GCF#data-dictionary) on the `lsaBGC-PopGene.py` wiki.

For more advanced evolutionary statistics, I recommend primarily consulting the Tajima's D and Median Beta-RD statistics, ***dN/dS & Fst are experimental implementations!!!*** and often times not very informative. 

**Median Beta-RD** is the median value of pairwise calculations of gene divergence to expected divergence for a pair of genomes based on genome-wide values such as ANI/AAI or divergence of single-copy core genes (e.g. ribosomal proteins).

**Tajima's D** can be most intiutively thought of as the proportion of high-frequency to low-frequency/rare minor allele positions along the multiple sequence alignment of a gene. It is thus essential to account for representation bias in your dataset and as such we perform de-replication (check out [drep](https://github.com/MrOlm/drep)) of genomes prior to `lsaBGC-PopGene.py`/`lsaBGC-AutoAnalyze.py` analysis to properly compute the statistic.

In the near future, I plan to incorporate a phylogeny-aware calculation of dN/dS, as well as update other evolutionary / pop-gen stats. Will announce these in major updates as they are incorporated.

A highlight of some key columns in the table after loading into Google sheets:

![image](https://user-images.githubusercontent.com/4260723/178626400-9563427f-d98d-48d6-8b4b-0e9673c6d0cd.png)

## Notes on Scalability:

Updated 06/19/2022

lsaBGC strives for large-scalability and is designed to be high-throughput. We previously ran it on >15,000 Staphylococci. However, it is recommended that the number of "Primary" genomes not exceed 300 (or 500 at most!) because de novo ortholog grouping is performed for this set (an all vs. all procedure). It is therefore recommended that users first select a distributed set of representative, "primary" genomes for their taxa of interest (which can be done via genome dereplication, check out [drep](https://github.com/MrOlm/drep)!). Then treat the rest of the genomes as "additional" (this can be 1000s of genomes!).

For evolutionary statistics calculations, codon alignments are now built using MAGUS (a divide and conquer wrapper of MAFFT) by default! 

## Future Updates Planned and of High Priority:

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

## Major Updates
* Feb 29, 2023 - Added GSeeF analysis to the end of the lsaBGC-Easy workflow and added support for parsing annotations from DeepBGC and GECCO into GSeeF.
* Feb 28, 2023 - Handle similarly named BGCs from different genomic assemblies (due to similar scaffold names) in lsaBGC-Ready.py/lsaBGC-Easy.py + automated ulimit setting updates in lsaBGC-Ready.py similar to what is employed in lsaBGC-Easy.py 
* Feb 5, 2023 - Corrected flag to make prodigal the default and pyrodigal optional unless requested by the user - we are considering formally making pyrodigal the default gene-caller however because there are a number of improvements/corrections pyprodigal offers and it is being actively maintained.
* Jan 16, 2023 - Added more adjustable parameters to lsaBGC-DiscoVary. Expansion uses bitscores to select best homolog group for genes from new genomes while cutoffs still based on E-values to provide more resolution and sensitivity in distinguishing NRPS/PKS homolog groups with recently diverged protocore genes (in alignment with recent v1.31 switch to use more resolute homolog groups).  
* Nov 25, 2022 - Much improved logging and standard output for lsaBGC-Easy (less yelling to the terminal). Recent switch to use hierarchical ortholog groupings found to be non-deterministic (still default - but can use coarser orthogroups via `-mc` argument for `lsaBGC-Easy.py` and `lsaBGC-Ready.py`). `GSeeF.py` made more robust and default is to only print tracks for most common 50 GCFs. 
* Nov 9, 2022 - Switch from using coarse othogrouping to hierarchical ortholog groupings in OrthoFinder analysis. Add option to provide user custom species tree. 
* Nov 6, 2022 - Introduce [`GSeeF.py`](https://github.com/Kalan-Lab/lsaBGC/wiki/17.-GSeeF---Visualizing-GCF-Cluster-Presence-and-Annotation-Along-a-Species-Phylogeny), improve interface of `lsaBGC-Easy.py`, and fix code for synchronizing gene-calling potential differences between BGCs and prodigal on full genomes (introduced in v1.2). 
* Aug 21, 2022 - Parsing of BiG-SCAPE results updated to account for GCFs spanning multiple class categories. All BGC instances in each class the GCF occurs will be conglomerated. Default usage of lsaBGC-Cluster.py is unaffected.
* Aug 16, 2022 - We recently (in version 1.1) had switched from scipy's `median_absolute_deviation` function to `median_abs_deviation` function - to mimic behavior of the original function and how we benchmarked lsaBGC in our manuscript, have set the scale of the function to `normal`. Though we are considering changing the scale to the default of the new function `1.0` in latter versions after further testing.
* Aug 08, 2022 - Added a couple more columns to PopGene report: `single-copy_in_GCF_context`, `max_beta_rd`. Also renamed column `hg median copy count` to a more appropriate label of `hg proportion multi-copy genome-wide` (calculation is the same!). Added `use_core_only` option to codon-alignment comparison methods to more appropriately calculate Beta-RD from recent switch to not regard sites where one sample/protein/gene has an allele but the comparing sample has a gap as a mismatch, but rather such sites are now ignored by default.  
* Aug 04, 2022 - Added some minor code to allow running `lsaBGC-Easy.py` with only an "Additional_Genomes" directory. Command would look something like: `lsaBGC-Easy.py -n "None" -g Genomes_Directory/ -o Results/`
* Aug 01, 2022 - Polished `lsaBGC-Easy.py` and created more streamlined version: see details on its [Wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial:-Combining-lsaBGC-with-ncbi-genome-download).
* Jul 28, 2022 - Introduced `lsaBGC-Easy.py` which simplifies usage of the suite. The general suite (besides `lsaBGC-DiscoVary.py`) can also now handle fungal genomes if BGCs are predicted using antiSMASH, but investigating fungal genomes is not currently an option in `lsaBGC-Easy.py`. Additional user interface upgrades + dependencies introduced, as such please reconfigure the conda environment and reinstall if you tried lsaBGC previously! Long-overdue but have a notice on citing the many dependencies now in the main folder of the repo titled: `CITATION_NOTICE`. 
* Jul 12, 2022 - Huge thanks to Martin Larralde for recommendation and advice on how to update GECCO processing to better identify "protocore"-esque genes in BGCs! Also, have now added a simplified sub-section down below with [notes on the PopGene report table](#notes-on-the-popgene-table-report).
* Jul 10, 2022 - Several updates made. Fixed small issues with smooth running of new framework, `lsaBGC-Ready.py`. Removed some dependencies and have added GToTree for creating species phylogeny + estimated sample to sample amino acid expected divergences. New small test dataset now included in this repo for immediate testing + much simplified installation guide. Most major change is that lsaBGC now works with DeepBGC and GECCO predictions! lsaBGC's backend relies on 'proto-core homolog groups' / 'rule-based key domains' determined by AntiSMASH, to get around the absence of such marker genes/domains in DeepBGC and GECCO predictions, domains in the highest 10% of deebgc_scpus or lowest 10% of e-values are treated as "proto-core" and used in `lsaBGC-AutoExpansion.py`/`lsaBGC-DiscoVary.py` as well as highlighted/treated as the "core" in `lsaBGC-PopGene.py` reports.
* Jun 26, 2022 - Added "loose" mode to `lsaBGC-Expansion.py` and option for users to manually define "protocore" homolog groups. Also, "protocore" homolog groups for a GCF now must have "rule-based" marker to exclude MGEs like transposons which insert within protocore regions of BGCs.
* Jun 19, 2022 - Have set MAGUS as the default protein alignment method (highly scalable wrapper of mafft) + updated notes on scalability.
* Jun 18, 2022 - Updated [`lsaBGC-AutoAnalyze.py`](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-The-lsaBGC-AutoAnalyze-Workflow) (automated lsaBGC analysis for each GCF) for better integration into new framework based around `lsaBGC-Ready.py`. 
* Jun 14, 2022 - Added [note on scalability](#user-content-notes-on-scalability), below on this page, and future plans to address them.
* Jun 09, 2022 - Fixed issues with `lsaBGC-Ready.py` & New Tutorial check it out [here](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-&-In-Depth-Tutorial:-Exploring-BGCs-in-Cutibacterium)!
* Jun 06, 2022 - Major updates to `lsaBGC-Ready.py` - the new recommended program for setting-up to run the lsaBGC suite.
* May 24, 2022 - `lsaBGC-Ready.py` is now available and can take pre-computed antiSMASH BGC predictions, along with optional BiG-SCAPE clustering results, to produce the required inputs for major lsaBGC analytical programs (`lsaBGC-See.py`, `lsaBGC-Refine.py`, `lsaBGC-PopGene.py`, `lsaBGC-DiscoVary.py`). 
