from setuptools import setup

setup(name='lsaBGC',
      version='1.0',
      description='Suite for comparative genomic, population genetics and evolutionary analysis, as well as metagenomic mining of micro-evolutionary novelty in BGCs all in the context of a single lineage of interest.',
      url='http://github.com/Kalan-Lab/lsaBGC/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['lsaBGC'],
      scripts=['scripts/createPickleOfSampleAnnotationListingFile.py', 'scripts/popSizeAndSampleSelector.py', 'workflows/lsaBGC-AutoProcess.py', 'workflows/lsaBGC-AutoAnalyze.py', 'workflows/lsaBGC-AutoExpansion.py', 'bin/lsaBGC-Cluster.py', 'bin/lsaBGC-See.py', 'bin/lsaBGC-PopGene.py', 'bin/lsaBGC-Refiner.py', 'bin/lsaBGC-Expansion.py', 'bin/lsaBGC-Divergence.py', 'bin/lsaBGC-DiscoVary.py'],
      zip_safe=False)
