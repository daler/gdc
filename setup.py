from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='GDC',
      version=version,
      description="Genomic dataset constructor",
      long_description="""\
Make BED, GFF, FASTA, SAM, and BAM files for testing your packages""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Ryan Dale',
      author_email='dalerr@niddk.nih.gov',
      url='',
      license='',
      packages=['gdc','gdc.test', 'gdc.scripts', 'gdc.test.data'],
      package_data={'gdc':['test/data/*']},
      package_dir={'gdc': 'gdc'},
      scripts=['gdc/scripts/gdc_parse_model.py'],
      zip_safe=False,
      install_requires=['argparse', 'GFFutils'],
      )
