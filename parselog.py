#!/usr/bin/python
# Parses dat files produced by stanalysis.python
import numpy
#import matplotlib.pyplot as plt

import json

# to dump numpy arrays
class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray) and obj.ndim == 1:
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

import sys # for handling command line arguments

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--bands', action = 'store_true', help = 'Process only basic bands information')
parser.add_argument('logFileName', help = 'The name of the log file to process')
args = parser.parse_args()

Tolerance = 1.0e-10
Tol2 = Tolerance**2

import re
def grep(pattern,fileObj):
    r=[]
    linenumber=0
    for line in fileObj:
        linenumber +=1
        if re.search(pattern,line):
            r.append((linenumber,line))
    fileObj.seek(0) # reset the cursor # redo with an array
    return r

def adjust(num) :
    tol = Tolerance
    return num if abs(num) > tol else 0.0*num

def adjustcomplex(num) :
    tol = Tolerance
    return adjust(num.real) + adjust(num.imag)*1j

def addcomplex (repart, impart) :
    # takes a complex number as a pair of strings 
    # combines them and rounds its parts to zero if they are too small

    out = adjustcomplex(complex(repart) + complex(impart))
    return out

def getVec(line) :
    # The format is
    # (...) :  [ -5.111e-15 +1.150e-45j   1.523e-15 +4.957e-45j  -3.970e-01 -2.342e-46j]
    # or
    # (...) :  [ 0.  0.  0.]
    #
    # or (notice spaces)
    # (...) :  [  1.987e-15+1.336j   9.039e-15-0.595j  -1.302e+00-0.33j ]
    place_start = line.find('[') + 1
    place_end = line.find(']')
    strPre = line[place_start:place_end]
    strNums = strPre.split()
    
    if complex(strNums[0]) == 0j :
        # we have case #2 but usually we don't get here
        return numpy.zeros(3)

    # A bit of traditional coding
    count = 0
    out = []
    while count < len(strNums) :
        try:
            partRe = float(strNums[count])
            count += 1
            partIm = complex(strNums[count])
            num = partRe + partIm
        except ValueError :
            num = complex(strNums[count])
        out.append(adjustcomplex(num))
        count += 1
    if len(out) != 3:
        sys.exit("Error! Parsing inconsistency")

    return numpy.array(out)

def getKVec(line) :
    # The format is
    # Bands at the k-point: [ 0.33333333  0.33333333  0.        ] are calculated
    
    place_start = line.find('[') + 1
    place_end = line.find(']')
    strNums = line[place_start:place_end].split()

    # Here everything is real
    return numpy.array([float(strNums[i]) for i in range(3)])
    
def addvalue(dataref, spin, upd) :
    if dataref[spin]['fill'] > Tolerance :
        dataref[spin].update(upd)

try :
    with open(args.logFileName,'r') as fileLog:
        # Extract k
        linesGet = grep('Bands at the k-point', fileLog)
        vec = getKVec(linesGet[0][1])
        data = [{'kPoint' : vec}]
        print "Data taken at %s is processed" % vec

        # First we extract energies
        linesGet = grep('is being processed', fileLog)
        # The output is a list of tuples
        # (line_number, 'Band 18 (Energy: 1.99234743534) is being processed\n')
        #data = []
        for line in linesGet :
            rec = []
            start = line[1].find(" ") +1
            end = line[1].find(" ", start)
            numBand = int(line[1][start:end].strip())
            start = line[1].find(":") + 2
            end = line[1].find(")") 
            energy = float(line[1][start:end].strip())
            add = {'band': numBand, 'energy' : energy, \
                   'Up' : {}, 'Down' : {}}
            # here making copy is not probably necessary but just to be on the safe side
            data.append(add.copy())

        totBands = len(data) - 1
        print "Totally %s bands are recovered" % totBands

        if not args.bands :

            # Then we have pairs Up/Down
            # Populations of the bands
            print "Processing populations"
            linesGet = grep('component has norm', fileLog)
            print "%s records are found" % len(linesGet)
            if len(linesGet) != 2*totBands :
                sys.exit('Mismatch of number of records. Log file is corrupted')

            flagUp = True
            bandCount = 1
            for line in linesGet :
                pop = adjust(complex(
                    line[1][
                        (line[1].find('(') + 1):
                        line[1].find(')')]).real)

                if flagUp :
                    data[bandCount]['Up'].update({'fill' : pop})
                else :
                    data[bandCount]['Down'].update({'fill' : pop})
                    bandCount += 1
                flagUp = not flagUp

            # Orbital momentum
            print "Processing orbital momenta"
            linesGet = grep('(Lx, Ly, Lz)', fileLog)
            print "%s records are found" % len(linesGet)
            if len(linesGet) != 2*totBands :
                sys.exit('Mismatch of number of records. Log file is corrupted')

            flagUp = True
            bandCount = 1
            for line in linesGet :
                vec = getVec(line[1])/1.0j
                # Here we may loose something
                vec = vec.real
                addvalue(data[bandCount], 
                         spin = 'Up' if flagUp else 'Down', upd = {'L' : vec})
                if not flagUp :
                    bandCount += 1

                flagUp = not flagUp

            # momentum
            print "Processing momenta"
            linesGet = grep('(px, py, pz)', fileLog)
            if len(linesGet) != 2*totBands :
                sys.exit('Mismatch of number of records. Log file is corrupted')

            flagUp = True
            bandCount = 1
            for line in linesGet :
                vec = getVec(line[1])/1.0j
                # Here we may loose something
                vec = vec.real
                addvalue(data[bandCount], 
                         spin = 'Up' if flagUp else 'Down', upd = {'P' : vec})
                if not flagUp :
                    bandCount += 1

                flagUp = not flagUp

            # coordinate
            print "Processing coordinates"
            linesGet = grep('(rx, ry, rz)', fileLog)
            if len(linesGet) != 2*totBands :
                sys.exit('Mismatch of number of records. Log file is corrupted')

            flagUp = True
            bandCount = 1
            for line in linesGet :
                vec = getVec(line[1]).real
                addvalue(data[bandCount], 
                         spin = 'Up' if flagUp else 'Down', upd = {'R' : vec})
                if not flagUp :
                    bandCount += 1

                flagUp = not flagUp
        # end of full processing

        # now we dump the processed data
        with open(args.logFileName + '.json', 'wt') as out:
            res = json.dump(data, out, sort_keys=True, indent=4, separators=(',', ': '), cls = NumpyAwareJSONEncoder)

except IOError :
    sys.exit("The log file couldn't be open")