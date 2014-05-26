# pyHgMapper

A simple python command line tool for mapping chromosome coordinates to genetic elements on Human (or other species) genome. 

Results include common summary statistics in various of plots; details of each indels in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1 "BED format") format; INDEL affected genes in GMT format, enrichment analysis is performed via [Enrichr](http://amp.pharm.mssm.edu/Enrichr/ "Enrichr") and links to the results also provided. 

## How To Install

Via git:

```
	$ git clone https://github.com/wangz10/pyHgMapper
```

## Dependent Python modules:

This tool is dependent on common scientific computing modules: Numpy, Scipy, Matplotlib. 

As well as other dependent modules:

[Pygal](http://pygal.org/ "Pygal") for drawing interactive SVG images:
```
	$ pip install pygal
```

Poster for streaming HTTP uploads:
```
	$ pip install poster
```

## How To Use

To load the genomic information of human genome, go to [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables? "UCSC Table Browser"),
set track to "RefSeq Genes", table to "refGene", and download the table as 
either plain text or gz format and save it into the same directory with the mapping.py script.

```
	$ python mapping.py [-h] [-a ANNOTATION] indel_filename
```
