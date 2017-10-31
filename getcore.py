#!/usr/bin/python3
""" a tool for extracting differences in pfam files from
	output of sort_refactored:
	
	three important remarks:

	##### 1 check that the common benchmarkstrain is hardcoded ####
	
	##### 2 check input #####
	 to get proper inpuput a composite file must be created of all .csv
	files with appended sequences with the command head:
	
	in file with output pfam.csv:
	$ head -x *.csv > file.txt

	where x is larger than the total amount of strains in the .csv files.

	usage python3 getcore.py reference.txt query_reference.txt [r/q]
	
	#### 3 to test the output: ####
	
	make the first comparison file.:
	$ python3 getcore.py allst.txt allst_ref.txt q > corequery.txt
	
	sort out unique hits:
	$ cat corequery.txt | grep None | cut -f2 >uniqcorequery.txt

	check that the sum is 0 when comparing number of unique + common core - single core.
	$ a=$(cat allst_ref.txt | grep = | wc -l); b=$(cat uniqcorequery.txt | wc -l ) ; c=$(cat corequery.txt| wc -l) ;echo $a+$b-$c | bc
	
"""


import sys
import re


def makedict_first_file(filename,bench):
	""" reads file to memory, compiles regex for two pattern, header and
	benchmark strain. reads file line by line
	"""
	with open (filename,"r") as handle:


		#lines = re.compile("==> pfam[0-9]*.csv <==").split(file)
		benchmark = re.compile(bench) 
		header = re.compile("==> pfam[0-9]*.csv <==")
		ndict = {}

		for line in handle:
			line = line.strip()

			if header.search(line):
				value = line.split(" ")[1]

			if benchmark.search(line):

				splitline = line.split("\t")
				seq =  splitline[3]
				ndict[seq] = value

	return ndict


def makedict_second_file(filename,bench):

	with open (filename,"r") as handle:


		#lines = re.compile("==> pfam[0-9]*.csv <==").split(file)
		benchmark = re.compile(bench)
		header = re.compile("==> pfam[0-9]*.csv <==")
		ndict = {}

		for line in handle:
			line = line.strip()

			if header.search(line):
				value = line.split(" ")[1]

			if benchmark.search(line):

				splitline = line.split("\t")
				seq =  splitline[3]
				ndict[seq] = value

	return ndict

comparisonflag = sys.argv[3]

if comparisonflag == "r":
	benchmark = "3655.faa"

elif comparisonflag == "q":
	benchmark = "11_S1.faa"
	
d1 = makedict_first_file(sys.argv[1], benchmark)
d2 = makedict_second_file(sys.argv[2], benchmark)

for i in d1:
	try:
		print("file1:\t%s\tfile2:\t%s" % (d1[i],d2[i]))
	except:
		print("file1:\t%s\tfile2:\t%s" % (d1[i],None))


