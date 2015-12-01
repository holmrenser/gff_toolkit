#!/usr/bin/python
__author__ = 'rensholmer'
__created__ = '17/08/15'

import sys
from setuptools import setup

def main():
	setup(name='gff_toolkit',
	      packages=['gff_toolkit'],
	      author='rens holmer',
	      author_email='rens.holmer@wur.nl',
	      version='0.1.5',
	      py_modules=['gff','gffsubpart','parser','test','__init__'],
	      url='https://github.com/holmrenser/gff_toolkit',
	      install_requires=['intervaltree']
	      )

if __name__ == '__main__':
	main()
