#!/usr/bin/env atkpython
# Analysis of Bloch states in complex structures
#
# For a given k-point it 
#   1. calculates all the bands in the provided range
#   2. For each band it calculates
#       2.1 Population for each spin projection (up/down)
#       2.2. Average coordinate, momentum, angular momentum
#
# TODO: implement arbitrary k-point to support encircling
# TODO: implement support for file management

from NanoLanguage import *

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-c', default = 'NULL', help = 'The name of the .nc-file with the configuration (results of the DFT calculations)')
parser.add_argument('-k', default = 'G', help = 'The k point (currently: G, K, Kbar only')
parser.add_argument('--bloch', action = 'store_true', help = 'Calculate only Bloch functions')
parser.add_argument('--size_up', type = int, default = 17, help = 'The number of bands above the Fermi level')
parser.add_argument('--size_down', type = int, default = 17, help = 'The number of bands below the Fermi level')
args = parser.parse_args()

sizeWindowUp = args.size_up
sizeWindowDown = args.size_down
sizeWindow = sizeWindowDown + sizeWindowUp

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

k_GPoint = numpy.array([0.0, 0.0, 0.0])
k_KPoint = numpy.array([0.33333333, 0.33333333, 0.0])

KPoint = k_GPoint if args.k == 'G' else \
    k_KPoint if args.k == 'K' else -k_KPoint

bulk_conf = nlread(confFileName, BulkConfiguration)[0]

print "Bands at the k-point: %s are calculated" % KPoint

bandstructure = Bandstructure(configuration = bulk_conf, 
    kpoints = [KPoint, KPoint], 
    bands_above_fermi_level = sizeWindowUp)

outBandFileName = "band_%s_%s.nc" % (args.k, namePostfix)
outFileName = "%s_%s.dat" % (args.k, namePostfix)

nlsave(outBandFileName, bandstructure)

print "Finding bands is completed. The Fermi level is at %s eV" % bandstructure.fermiLevel().inUnitsOf(eV)
enarray = numpy.array(bandstructure.evaluate()[0])

refEnergy = 0

listBelow = numpy.where(enarray < refEnergy)[0]
listAbove = numpy.where(enarray > refEnergy)[0]

b_list = numpy.concatenate([listBelow[-sizeWindowDown:], listAbove])
numStates = len(b_list)

print "Within the chosen range of bands the energies are distributed as follows"

for i in b_list :
    print "Band %s: %s" % (i, enarray[i])

# FileLog(Kpoint, b_list, energies[b_list])

def ccross(a, b) :
    """ Calculates cross product of vectors whose components are 
    complex numbers. The standard numpy.cross cannot deal with them.
    """
    return numpy.array([a[1]*b[2] - b[1]*a[2], \
        a[2]*b[0] - b[2]*a[0], \
        a[0]*b[1] - b[0]*a[1]])

BlochFileName = "Bloch_%s_%s.nc" % (args.k, namePostfix)

# arrays to store R, P, L
L = numpy.zeros([numStates, 3, 2], dtype = complex)
P = numpy.zeros([numStates, 3, 2], dtype = complex)
R = numpy.zeros([numStates, 3, 2])

spins = [Spin.Up, Spin.Down]

Norms = numpy.zeros(2)

numpy.set_printoptions(precision=3)

def addrecord(num, sp, Ns, l, p, r):
    l = l/(1.0*Ns)
    p = p/(1.0*Ns)
    r = r/(1.0*Ns)
    L[num, :, sp] = l
    P[num, :, sp] = p
    R[num, :, sp] = r
    
    print "(Lx, Ly, Lz) : ", l
    print "(px, py, pz) : ", p
    print "(rx, ry, rz) : ", r
    print " <r> x <p> :", ccross(r, p)

strSpins = ['Up', 'Down']

def processBloch(count, state, spinState) :
    # constant vector of zeros to pass instead of actual small data
    zers = numpy.zeros(3)

    spinIndex = 0 if spinState == Spin.Up else 1

    state  = state.spinProjection(spinState)
    psi = state.toArray()
    Norm = (psi*psi.conj()).sum()
    print "%s component has norm: %s" % (strSpins[spinIndex], Norm)
    if Norm < 1e-7 :
        addrecord(count, spinIndex, 1, zers, zers, zers)
    else :
        psi_r = (psi + psi.conj())/2.0
        psi_i = (psi - psi.conj())/2.0/1j

        print "The magnitudes of the real and imaginary parts: %s, %s" % ((psi_r*psi_r).sum(), (psi_i*psi_i).sum())
        # The gradient of the wave function and such
        rdpsi = numpy.zeros(psi.shape + (3,), dtype = complex)
        dpsi = numpy.zeros(psi.shape + (3,), dtype = complex)
        rpsi = numpy.zeros(psi.shape + (3,), dtype = complex)

        ni, nj, nk = state.shape()
        irange = range(ni)
        jrange = range(nj)
        krange = range(nk)
        for ii in irange:
            for jj in irange :
                for kk in krange :
                    rr = state.gridCoordinate(ii, jj, kk)
                    d = state.derivatives(rr[0], rr[1], rr[2]).tolist()
                    # to not have problems with complex multiplication
                    r = numpy.array(rr.tolist()) 
                    rdcross = ccross(r,d)
                    dpsi[ii,jj,kk,:] = d                    
                    rdpsi[ii,jj,kk,:] = rdcross
                    rpsi[ii,jj,kk,:] = psi[ii,jj,kk]*r
        # loop over the grid
        p = [(dpsi[:,:,:, k]*psi[:,:,:].conj()).sum() for k in range(3)]
        r = [(rpsi[:,:,:, k]*psi[:,:,:].conj()).sum() for k in range(3)]
        l = [(rdpsi[:,:,:, k]*psi[:,:,:].conj()).sum() for k in range(3)]
        addrecord(count, spinIndex, Norm, l, p, r)

print "Analysing Bloch states..."

for i in range(numStates) :
    print "Band %s (Energy: %s) is being processed" %(b_list[i], enarray[b_list[i]])
    bState =  BlochState(configuration=bulk_conf, quantum_number=int(b_list[i]),
        k_point=KPoint)
    nlsave(BlochFileName, bState, comment = str(b_list[i]) + 'th Bloch function')
    print "Bloch state is extracted"
    # This way is faster for some reason than turning the whole state into array
    if not args.bloch :
        processBloch(i, bState, Spin.Up)
        processBloch(i, bState, Spin.Down)

# a bit of ugly coding (TODO: correct this mess)        
if not args.bloch :
    with open("LS_"+outFileName, "w") as Lfile :
        for i in range(numStates) :
            Lfile.write("band: %s (KPoint: %s Energy: %s), Up, " 
                % (b_list[i], KPoint, enarray[b_list[i]]) )
            Lfile.write(", ".join(str(x) for x in L[i, :, 0]) )
            Lfile.write("\n")
            Lfile.write("band: %s (KPoint: %s Energy: %s), Down, " 
                % (b_list[i], KPoint, enarray[b_list[i]]) )
            Lfile.write(", ".join(str(x) for x in L[i, :, 1]) )
            Lfile.write("\n")

    with open("DS_" + outFileName, "w") as Lfile :
        for i in range(numStates) :
            Lfile.write("band: %s (KPoint: %s Energy: %s), Up, " 
                % (b_list[i], KPoint, enarray[b_list[i]]) )
            Lfile.write(", ".join(str(x) for x in P[i, :, 0]) )
            Lfile.write("\n")
            Lfile.write("band: %s (KPoint: %s Energy: %s), Down, " 
                % (b_list[i], KPoint, enarray[b_list[i]]) )
            Lfile.write(", ".join(str(x) for x in P[i, :, 1]) )
            Lfile.write("\n")

    with open("RS_" + outFileName, "w") as Lfile :
        for i in range(numStates) :
            Lfile.write("band: %s (KPoint: %s Energy: %s), Up, " 
                % (b_list[i], KPoint, enarray[b_list[i]]) )
            Lfile.write(", ".join(str(x) for x in R[i, :, 0]) )
            Lfile.write("\n")
            Lfile.write("band: %s (KPoint: %s Energy: %s), Down, " 
                % (b_list[i], KPoint, enarray[b_list[i]]) )
            Lfile.write(", ".join(str(x) for x in R[i, :, 1]) )
            Lfile.write("\n")
# end of saving results of processing Bloch states