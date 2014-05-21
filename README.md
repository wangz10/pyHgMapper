# pyHgMapper

A simple python command line tool for mapping chromosome coordinates to genetic elements on Human (or other species) genome.

# How To Install

Via git:

```
	$ git clone https://github.com/wangz10/pyHgMapper
```

## Dependent Python modules:

[Pygal](http://pygal.org/ "Pygal") for drawing interactive SVG images:
```
	$ pip install pygal
```

Poster for streaming HTTP uploads:
```
	$ pip install poster
```

# How To Use

```
	$ mapping.py [-h] [-a ANNOTATION] indel_filename
```

To load the genomic information of human genome, go to [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables? "UCSC Table Browser"),
set track to "RefSeq Genes", table to "refGene", and download the table as 
either plain text or gz format and save it into the same direcotory with the mapping.py script.

# More Documentations...
