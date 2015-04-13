#!/usr/bin/env python
# Shows data collected in the json file by parselog/parsedat
import numpy
import matplotlib.pyplot as plt

import json

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--old', action='store_true', help = 'Set for flipping first two observers')

parser.add_argument('-c', default = 'NULL', help = 'The name of the .nc-file with the configuration')
parser.add_argument('--size_down', default = 17, help = 'The number of bands below the Fermi level')

parser.add_argument("jsonFileName", help = 'The name of the input json file')
args = parser.parse_args()

import sys # 

Tolerance = 1.0e-10
Tol2 = Tolerance**2

try:
	with open(args.jsonFileName, "rt") as fileJson :
		data = json.load(fileJson)
except IOError:
	raise

# data contains a list of dicts
# each dict is of the format
# {"Down": {"L": [..], "P": [...], "R": [...],  "fill": 73.03},
#   "Up": {....}
#   "band": 1,
#	"energy": -12.7425419539 }
# if "fill" == 0 then LPR are absent (can be taken as zero)
# the records 'Down' and 'Up' can be empty

# Energies vs band number
# Markers:
#	black - both spins
#	blue - spin down only
#	red - spin up only

kPoint = data[0]['kPoint']
data = data[1:]

bands = numpy.array([data[i]['band'] for i in range(len(data))])
energies = numpy.array([data[i]['energy'] for i in range(len(data))])

colors = []
flagBands = False 

for i in range(len(data)) :
	if len(data[i]['Up']) == 0 or len(data[i]['Down']) == 0 :
		colors.append('black')
		flagBands = True
		continue

	if data[i]['Up']['fill'] > Tolerance and data[i]['Down']['fill'] > Tolerance :
		colors.append('black')
	elif data[i]['Up']['fill'] > Tolerance :
		colors.append('red')
	else :
		colors.append('blue')

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.scatter(bands, energies, c = colors, lw = 0)
ax1.set_xlim([bands[0], bands[-1]])
ax1.set_ylim([energies[0], energies[-1]])
ax1.set_xlabel('Bands')
ax1.set_ylabel('Energies, eV')
ax1.set_title('Section at k = %s' % kPoint)

outFileName = args.jsonFileName[:args.jsonFileName.find('.json')] + ".png"
fig1.savefig("energy_" + outFileName, bbox_inches='tight')

if flagBands :
	sys.exit()

def getLadj (rec, spin) :
	l = numpy.array(rec[spin]['L'])
	r = numpy.array(rec[spin]['R'])
	p = numpy.array(rec[spin]['P'])
	return l - numpy.cross(r,p)

def getR (rec, spin) :
	r = numpy.array(rec[spin]['R'])
	return r

# adjusted angular momentum and coordinates
Lz = []
R = []

for i in range(len(data)):
	if data[i]['Up']['fill'] > Tolerance :
		l = getLadj(data[i], 'Up')
		r = getR(data[i], 'Up')
	else :
		l = getLadj(data[i], 'Down')
		r = getR(data[i], 'Down')

	if abs(l[0]) + abs(l[1]) > abs(l[2])/10 :
		print "Band %s has noticeable in-plane component" % data[i]['band']
	Lz.append(l[2])
	R.append(r.tolist())

R = numpy.array(R)
Rzav = (R[:,2]).sum()/len(R)

for i in range(len(data)) :
	if abs(R[i][2] - Rzav) > 0.1 :
		print "Significant vertical variation in band %s" % data[i]['band']

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

ax2.scatter(bands, Lz, c = colors, lw = 0)
ax1.set_xlim([bands[0], bands[-1]])
#ax2.set_ylim([energies[0], energies[-1]])
ax2.set_xlabel('Bands')
ax2.set_ylabel('Orbital angular momentum, eV')
ax2.set_title('Section at k = %s' % kPoint)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.scatter(R[:,0],R[:,1], c = colors, lw = 0)
#ax3.set_xlim([bands[0], bands[-1]])
#ax2.set_ylim([energies[0], energies[-1]])
ax3.set_xlabel('x coordinate of CM')
ax3.set_ylabel('y coordinate of CM')
ax3.set_title('Section at k = %s' % kPoint)

s = - 1
for band, x, y in zip(bands, R[:,0], R[:,1]) :
	x = float(x)
	y = float(y)
	ax3.annotate(band, xy = (x, y), xytext = (s*10, 5), 
		textcoords = 'offset points')
	s = -s

fig2.savefig("orbital_" + outFileName, bbox_inches='tight')
fig3.savefig("cm_" + outFileName, bbox_inches='tight')
#plt.show()
