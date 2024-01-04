import os
import sys
from setuptools import setup

setup(name='lsaBGC',
      version='1.52',
      description='Suite for comparative genomic, population genetics and evolutionary analysis, as well as metagenomic mining of micro-evolutionary novelty in BGCs all in the context of a single lineage of interest.',
      url='http://github.com/Kalan-Lab/lsaBGC/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['lsaBGC'],
      scripts=['docker/LSABGC',
               'scripts/setup_annotation_dbs.py',
               'scripts/setup_bigscape.py',
               'scripts/GSeeF.py',
               'scripts/compareBGCtoGenomeCodonUsage.py',
               'scripts/listAllGenomesInDirectory.py',
               'scripts/listAllBGCGenbanksInDirectory.py',
               'scripts/runProdigalAndMakeProperGenbank.py',
               'scripts/visualize_BGC-Ome.py',
               'scripts/processAndReformatNCBIGenbanks.py',
               'scripts/genbankToProkkaGFF.py',
               'scripts/simpleProteinExtraction.py',
               'scripts/createPickleOfSampleAnnotationListingFile.py',
               'scripts/popSizeAndSampleSelector.py',
               'scripts/readifyAdditionalGenomes.py',
               'workflows/lsaBGC-Easy.py',
               'workflows/lsaBGC-Euk-Easy.py',
               'workflows/lsaBGC-AutoAnalyze.py',
               'workflows/lsaBGC-AutoExpansion.py',
               'bin/lsaBGC-Ready.py',
               'bin/lsaBGC-Cluster.py',
               'bin/lsaBGC-See.py',
               'bin/lsaBGC-ComprehenSeeIve.py',
               'bin/lsaBGC-PopGene.py', 
               'bin/lsaBGC-Refiner.py', 
               'bin/lsaBGC-Expansion.py', 
               'bin/lsaBGC-Divergence.py', 
               'bin/lsaBGC-DiscoVary.py',
               'bin/lsaBGC-MIBiGMapper.py'],
      zip_safe=False)

try:
      os.system('pip install cython==3.0.0a10')
      os.system('pip install sonicparanoid==2.0.2')
      os.system('sonicparanoid-get-test-data -o .')
      os.system('sonicparanoid -i sonicparanoid_test/test_input -o sonicparanoid_test/test_output -p my_first_run')
      os.system('rm -rf sonicparanoid_test/')
      os.system('pip install -q --force-reinstall -v "setuptools==58.2.0" 2> /dev/null')
except:
      sys.stderr.write('Warning: unable to install sonicparanoid!\n')
      sys.exit(1)
