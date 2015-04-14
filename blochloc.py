#!/usr/bin/env atkpython
# Analysis of the localization properties of Bloch states in complex structures
# NOT COMPLETED
#
# It deals with in-plane distributions.
#
# Requires:
#   1. File with DFT calculations
#   2. File with Bloch functions (TODO: make calculations from the scratch)
#   3. File with the bandstructure (TODO: the same as above)
#
# Calculates:
#
# 1. Integrated over z distribution
# 2. Center of mass of the in-plane distribution
# 3. Its dispersion
# 4. Its inverse participation ratio
#
# TODO: implement support for file management

Tolerance = 1e-10

from NanoLanguage import *

import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-c', default = 'NULL', 
		help = 'The name of the .nc-file with the configuration (results of the DFT calculations)')
parser.	add_argument('-b',
		help = 'The name of the file containing Bloch functions')
parser. add_argument('-s',
        help = 'The name of the file containing bandstructure')

args = parser.parse_args()

flagReadMain = True

if args.c == 'NULL' :
    print "The name of the file with DFT calculations is assumed to be in main.txt"
    # default file with the bulk configuration
    if flagReadMain :
    # We extract the name of file with the DFT results
        try :
            with open('main.txt') as file:
                confFileName =  file.readline()
        except IOError as e:
            print "Unable to open file 'main.txt' with the reference to DFT calculations"
            raise
    else :
        # Or use the given file name
        confFileName = "Mo_2_BS_DOS_OS_QAD.nc"
else :
    confFileName = args.c

namePostfix = confFileName[:confFileName.find(".nc")]

bulk_conf = nlread(confFileName, BulkConfiguration)[0]

# This requires a lot of memory
blochStates = nlread(args.b, BlochState)
numStates = len(blochStates)
print "Bloch states are read. Total %s records" % numStates

# DEBUG version
DEB_flag = False
if DEB_flag :
    print "Debug version is used: list of states is truncated"
    numStates = 1

# number of fields in the record
# |kpoint|, band, energy, IPR
numFields = 1+1+1+1 
listStates = numpy.zeros([numStates, numFields])

bandstructures = nlread(args.s)
numBStruct = len(bandstructures)
# There are few potential problems here.
# 1. It cannot be guaranteed that the first record in the bandstructure
#   file contains the bandstructure corresponding to the file with 
#   Bloch functions. Those files by default collect everyting written
#   into them. Since the bandstructure is needed for information 
#   purposes only, one shoudl allow for some slight discrepancy
#
# 2. We should also allow for possibility that Bloch states correspond
#   to different k points (TODO: implement this fully)

print "Bloch functions correspond to the following states"

for i in range(numStates):
    stateKpoint = blochStates[i].kPoint()
    bandNumber = blochStates[i].quantumNumber()
    listStates[i, 0] = numpy.sqrt((stateKpoint*stateKpoint).sum())
    listStates[i, 1] = bandNumber

    print "State #%s: k-point: %s, Band %s" % \
        (i, stateKpoint, bandNumber)

    # Now we look for respective bandstructures
    # The bandstructure is always calculated at least at two points
    # (the 'feature' of their Bandstructure calculator)

    for j in range(numBStruct) :
        for ind, k in enumerate(bandstructures[j].kpoints()) :
            df = k - stateKpoint
            if (df*df).sum() < Tolerance :
                break
        else :
            continue

        ebands = bandstructures[j].evaluate()[ind,:]
        if len(ebands) > bandNumber :
            listStates[i,2] = ebands[ind]
            print "Energy: %s" % ebands[ind]
            break
    else :
        print "ERROR: Couldn't find the bandstructure corresponding to this Bloch function"
        sys.exit(1)

outFileName = "IPR_%s_%s.dat" % (args.k, namePostfix)

spins = [Spin.Up, Spin.Down]
strSpins = ['Up', 'Down']

Norms = numpy.zeros(2)

numpy.set_printoptions(precision=3)

# Now we have all preliminaries satisfied

# taken from visBloch-red.py (TODO: modularize fcs)
def planarDensity(bState, ind) :
    """
    Calculates the density as the function of the in-plane coordinates.

    Returns 
        n_xy, x_cell, y_cell --- 2D numpy array, and x and y coordinates of the grid points
    inside the elementary cell (in Angstrom).

    If the density is too small, the function returns zero, None, None

    If 0 < ind < 1, it's taken as the relative z-coordinate and the function
    outputs the distribution of the density at this plane. Otherwise,
    it calculates the integral over the z axis.
    """
    # The main ideas are taken from
    # http://quantumwise.com/documents/manuals/ATK-13.8/ReferenceManual/ref.electrondifferencedensity.html

    ni, nj, nk = bState.shape()
    # Find the volume elements.
    dX, dY, dZ = bState.volumeElement()
    # Calculate the unit area in the y-z plane.
    length_unit = dX.unit()
    dAXY = numpy.linalg.norm( numpy.cross(dX,dY) )*length_unit**2
    dz = dZ.norm()

    psi = bState.toArray()
    density = (psi*psi.conj()).real

    if density.sum() < Tolerance :
        return 0, None, None

    if 0 < ind < 1 :
        n_xy = numpy.array(density[:,:,int(nk*ind)])
    else :
        n_xy = numpy.array(
            [ [ density[ii, jj, :].sum() 
                for ii in range(ni) ] for jj in range(nj) ])

    x_cell = numpy.array([[bState.gridCoordinate(ii, jj, 0).inUnitsOf(Angstrom)[0] 
                           for ii in range(ni)] for jj in range(nj)])
    y_cell = numpy.array([[bState.gridCoordinate(ii, jj, 0).inUnitsOf(Angstrom)[1] 
                           for ii in range(ni)] for jj in range(nj)])

    return n_xy, x_cell, y_cell


def processBloch(state, spinState) :
    spinIndex = 0 if spinState == Spin.Up else 1

    n, x, y = planarDensity(state.spinProjection(spinState), -1)

    if x == None :
        norm = 0
        ipr = 0
    else :
        norm = n.sum()
        n = n/norm
        pp = (n*n).sum()
        ipr = 1/pp
    return norm, ipr

print "Analysing Bloch states..."

for i in range(numStates) :
    print "Band %s (Energy: %s) is being processed" % (listStates[i,0], listStates[i,1])

    bState = blochStates[i]

    # This way is faster for some reason than turning the whole state into array
    norm, ipr = processBloch(i, bState, Spin.All)
    # norm, ipr = processBloch(i, bState, Spin.Down)
    print "Norm: %s, IPR: %s" % (norm, ipr)
    listStates[i,2] = ipr

with open(outFileName, "w") as iprfile :
    for i in range(numStates) :
        iprfile.write("band: %s, KPoint: %s, Energy: %s, IPR: %s " 
                % (listStates[i,0], KPoint, listStates[i,1], listStates[i,2]) )