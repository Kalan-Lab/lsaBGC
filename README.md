# *lsa*BGC
### Lineage Specific Analysis of Biosynthetic Gene Clusters

*lsa*BGC is a modular software suite designed to provide a comprehensive set of functions for investigating and mining for 
biosynthetic gene cluster diversity across a focal lineage/taxa of interest (currently limited to bacterial) using AntiSMASH 
based annotation. It consists of 8 independent programs: `lsaBGC-Ready.py`, `lsaBGC-Cluster.py`, `lsaBGC-Expansion`, 
`lsaBGC-Refine.py`, `lsaBGC-See.py`, `lsaBGC-PopGene.py`, `lsaBGC-Divergence.py`, and `lsaBGC-DiscoVary.py`.

![](https://github.com/Kalan-Lab/lsaBGC/blob/main/docs/images/lsaBGC1.1_Simplified.png)

## Major Updates 
* Jun 18, 2022 - Updated `lsaBGC-AutoAnalyze.py` for better integration into new framework based around `lsaBGC-Ready.py`.
* Jun 14, 2022 - Added [note on scalability](#user-content-notes-on-scalability), below on this page, and future plans to address them.
* Jun 09, 2022 - Fixed issues with `lsaBGC-Ready.py` & New Tutorial check it out [here](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Tutorial:-Exploring-BGCs-in-Cutibacterium)!
* Jun 06, 2022 - Major updates to `lsaBGC-Ready.py` - the new recommended program for setting-up to run the lsaBGC suite.
* May 24, 2022 - `lsaBGC-Ready.py` is now available and can take pre-computed antiSMASH BGC predictions, along with optional BiG-SCAPE clustering results, to produce the required inputs for major lsaBGC analytical programs (`lsaBGC-See.py`, `lsaBGC-Refine.py`, `lsaBGC-PopGene.py`, `lsaBGC-DiscoVary.py`). 

## Documentation:

Documentation can currently be found on this Github repo's wiki: https://github.com/Kalan-Lab/lsaBGC/wiki

1. [Background on lsaBGC - what it does and does not do](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background)
2. [Installation Guide](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation)
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

Installation is performed via Conda and should take ~5-10 minutes. We are happy to attempt to address issues with installation if any arise, please open a Git Issues case.

On the [Installation wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation), you can find a step-by-step guide and bash scripts for automated installation are provided. [Test cases](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases) to demonstrate individual programs are available along with a more [comprehensive tutorial](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Tutorial:-Exploring-BGCs-in-Cutibacterium) to showcase the use of the suite and relations between core programs.

## Notes on Scalability:

Currently, the primary/core set of samples/genomes used to define BGCs should not exceed 300 samples. This is because OrthoFinder2 performs all-vs-all alignments between genomes and will produce a lot of files. To resolve this limitation in the future, I am looking into developing stratedgies to make OrthoFinder2 more scalable (such as performing all-vs-all diamond alignments via a reflexive alignment to a single giant file with all proteins from all genomes - diamond gains a speed boost + avoids writing a lot of files) or using alternate software, though OrthoFinder2 has some immense benefits for accurate clustering that would be great to retain. For most genera, the 300 upper limit should work well after dereplication of genomes to remove redundancy. lsaBGC-AutoExpansion can then be used at high-scale to find homologous instances to GCFs defined from primary genomes in addtiional/draft genomes. 

Additionally, currently `lsaBGC-PopGene.py` can also take a long time to run when more than 300 genomes are provided. Dereplication of input genomes should definitely be performed prior to using this program however to more appropriately caclulate evolutionary and population genetic statistics. We plan to resolve this limitation regardless using [MAGUS](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008950) in the near future, which allows for fast alignment of up to a million sequences through a divide-and-conquer type approach. 

## Future Updates Planned and of High Priority:

* Update algorithms to have options to work with DeepBGC + GECCO BGC predictions.
* Incorporate eukaryotic gene calling in `lsaBGC-Ready.py` to allow application to fungi + plants.
* Incorporate MAGUS for protein sequence alignments. Make codon sequence alignments a seperate program. 

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
