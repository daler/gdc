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
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
