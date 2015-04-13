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

bandstructure = nlread(args.s)[0]

# The bandstructure is always calculated at least at two points
# (the 'feature' of their Bandstructure calculator)

if len(bandstructure.kpoints()) > 2 :
    print "WARNING: possible inconsistency between bandstructure and Bloch functions"
    print "Energy assignment is not reliable"

KPoint = bandstructure.kpoints()[0]

print "Bands at the k-point: %s are read" % KPoint
print "The Fermi level is at %s eV" % bandstructure.fermiLevel().inUnitsOf(eV)
enarray = numpy.array(bandstructure.evaluate()[0])

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
# band, energy, IPR
numFields = 1+1+1 
listStates = numpy.zeros([numStates, numFields])

print "Bloch states correspond to the following bands"

for i in range(numStates):
    bandNumber = blochStates[i].quantumNumber()
    bandEnergy = enarray[bandNumber-1]
    listStates[i, 0] = bandNumber
    listStates[i, 1] = bandEnergy

    print "Band %s: %s" % (bandNumber, bandEnergy)

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