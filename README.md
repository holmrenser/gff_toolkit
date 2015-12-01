[![Build Status](https://travis-ci.org/holmrenser/gff_toolkit.svg)](https://travis-ci.org/holmrenser/gff_toolkit)
[![Coverage Status](https://coveralls.io/repos/holmrenser/gff_toolkit/badge.svg?branch=master&service=github)](https://coveralls.io/github/holmrenser/gff_toolkit?branch=master)
[![PyPI version](https://badge.fury.io/py/gff_toolkit.svg)](https://badge.fury.io/py/gff_toolkit)
# gff_toolkit
tools for handling gff files

To install gff_toolkit simply run:
```pip install gff_toolkit```

To upgrade after a new version run:
```pip install gff_toolkit --upgrade```

```python
#!/usr/bin/python
"""
Example on how to reformat manual annotations from Geneious to proper Gff
"""

__author__ = 'rensholmer'

import sys
import gff_toolkit as gt

def main(manual_file,author_name):
	manual_gff = gt.parser(manual_file,filetype='manual',limit=dict(featuretype='CDS',source='Geneious'),author=author_name)
	for gene in manual_gff.getitems(featuretype='gene'):
		for subs in manual_gff.get_children(gene):  
			print subs.stringify().strip()

if __name__ == '__main__':
	if len(sys.argv) == 3:
		main(*sys.argv[1:])
	else:
		print 'Usage: python {0} <geneious manual annotation gff> <author name>'.format(sys.argv[0])
```
