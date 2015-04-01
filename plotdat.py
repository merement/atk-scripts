#!/usr/bin/python
# For plotting dat-files
# The file is assumed to have several columns of data
# within each column the data is separated by commas (that is
# the file is taken as CSV)
#
# The first column may contain text, in this case it's ignored
# The first numerical column is taken as X, the remaining columns are Y
#
# TODO: option for plotting columns themselves

import numpy
import matplotlib.pyplot as plt

import pandas # for reading csv

import sys # for handling command line arguments

def is_number(s):
	# taken from http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-in-python
	try:
		float(s)
		return True
	except ValueError:
		return False
        
def main(argv) :
	
	if len(argv) < 1 :
		print "Plots comma separated dat files. Requires an input file"
		return -1
		
	print "Processing ", argv
	try:
		dat = pandas.read_csv(argv[0], delimiter = ',').values
	except IOError :
		print "File couldn't be open. Wrong name?"
		return -1
	
	# we want a bit of extra functionality
	if len(argv) > 2 :
		xcol = argv[2]
	else :
		xcol = 0
	if len(argv) > 1 :
		ycol = argv[1]
	else :
		ycol = 1
	
	if is_number(dat[0,0]) :
		xind = xcol
		yind = ycol
	else :
		xind = xcol + 1
		yind = ycol + 1
		
	plt.plot(dat[:,xind], dat[:,yind], '.-')
	plt.show()

if __name__ == "__main__" :
	main(sys.argv[1:])


