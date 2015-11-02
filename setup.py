#!/usr/bin/python
__author__ = 'rensholmer'
__created__ = '17/08/15'

import sys
from distutils.core import setup


def main():
	setup(name='gff_toolkit',
	      author='rens holmer',
	      version='0.1',
	      py_modules=['gff','gffsubpart','parser','test','__init__']
	      )

if __name__ == '__main__':
	main()