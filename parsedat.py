#!/usr/bin/env python
# Parses dat files produced by stanalysis.python
#
# TODO: field tests

import numpy
import matplotlib.pyplot as plt

import pandas # for reading csv

import sys # for handling command line arguments

if len(sys.argv) < 2 :
	sys.exit("Basic analysis of comma separated dat files. Requires an input file")

flagOld = True # for files with screwed R part
if len(sys.argv) > 2 :
	if sys.argv[2] == 'new' :
		flagOld = False
	elif sys.argv[2] == 'old' :
		flagOld = True


# The procedure is
# 1. Look for one of the three patterns in the name '_G_', '_K_', '_Kbar_'
# 2. Completes to the full set
# 3. Takes out the first letter in each name and substitutes 
#		'L' (angular momentum)
#		'D' (momentum)
#		'R' (position)

inFileName = sys.argv[1]

inPatt = ['_G_', '_K_', '_Kbar_']
repPrefix = ['L', 'D'] #, 'R'] # we exclude 'R' for a while due to the error in the output

for pat in inPatt :
	if inFileName.find(pat) > 0 :
		break
else :
	sys.exit("The input file doesn't look like produced by stanalysis.py")

# We need to look for data belonging to the same group

def parsedat(fName) :
	def adjust (val) :
		# takes a complex number as a string and rounds its parts to zero if 
		# they are too small
		tol = 1.0e-10

		num = complex(val)
		if abs(num.real) < tol :
			num = 0 + 1j*num.imag
		if abs(num.imag) < tol :
			num = num.real
		return num

	print "Processing ", fName
	try:
		dat = pandas.read_csv(fName, delimiter = ',').values
	except IOError :
		print "Warning! File couldn't be open. Skipping"
		return {}

	nlines, ncols = dat.shape
	# ncols should be 5 

	# Some parsing
	# The format of the line
	# band: 1 (KPoint: [ 0.333  0.333  0.   ] Energy: -12.7425419539), Up, X, Y, Z
	# we want to return the array of dictionaries
	# {	'band' : N, 'KPoint' : [x, y, z], 'Energy': E, 'Spin': Up, 'data': [X, Y, Z]}
	parsed_data = []
	curK = numpy.zeros(3)
	for i in range(nlines):
		out = {}
		place_start = dat[i,0].find(':') + 1
		place_end = dat[i,0].find('(') - 1
		out['band'] = int(dat[i,0][place_start:place_end].strip())

		place_start = dat[i,0].find('[') + 2
		place_end = dat[i, 0].find(']') - 1
		# we need to take care of possible multiple spaces and to add the terminating space
		# Why not just split????
		kstr = ' '.join(dat[i,0][place_start:place_end].split()) + ' '
		p_start = 0
		for j in range(3) :
			p_end = kstr.find(' ', p_start)
			#stri = kstr[p_start: p_end].strip()
			#print stri
			curK[j] = float(kstr[p_start: p_end].strip())
			p_start = p_end + 1
		out['KPoint'] = curK # since we assign to curK it's safe

		place_start = dat[i, 0].find(':', place_end) + 1
		place_end =  dat[i, 0].find(')') - 1
		out['Energy'] = float(dat[i,0][place_start: place_end].strip())

		out['Spin'] = dat[i,1].strip()

		out['data'] = numpy.array([adjust(val) for val in dat[i, 2:]])
		parsed_data.append(out.copy()) # but here we have to make a copy
	return parsed_data


for point in inPatt :
	curKPostfix = inFileName.replace(pat, point, 1)
	data = {}
	for quant in repPrefix :
		listFiles = quant + curKPostfix[1:]

		data[quant] = parsedat(listFiles)
		data[quant]['File'] = listFiles

	# Nah, this is stupid
	#if flagOld :
	#	coorName = 'RR' + curKPostfix[1:]
	#	with open(coorName, 'r') as cfile:


	print data