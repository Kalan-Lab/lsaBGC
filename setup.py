from setuptools import setup

setup(name='lsaBGC',
      version='1.0',
      description='Suite for comparative genomic, population genetic, and metagenomic mining of BGCs in the context of single lineages.',
      url='http://github.com/Kalan-Lab/lsaBGC/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['lsaBGC'],
      scripts=['bin/createInputsForLsaBGC.py', 'bin/summarizeLsaBGCReporterResults.py', 'bin/lsaBG-See.py',
			   'bin/lsaBGC-Clust.py', 'bin/lsaBGC-HomologReports.py', 'bin/lsaBGC-MetaNovelty.py',
			   'bin/lsaBGC-Select.py', 'bin/lsaBGC-HMMExpander.py'],
      zip_safe=False)