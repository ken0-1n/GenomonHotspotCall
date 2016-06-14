#!/usr/bin/env python

from distutils.core import setup

setup(name='hotspotCall',
      version='0.1.0',
      description='Python tools for filtering somatic mutations using beta-binomial sequencing error model.',
      author='Kenichi Chiba',
      author_email='kchiba@hgc.jp',
      url='',
      package_dir = {'': 'lib'},
      packages=['hotspotCall'],
      scripts=['hotspotCall'],
      license='GPL-3'
     )

