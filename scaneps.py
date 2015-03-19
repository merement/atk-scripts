#!/usr/bin/env atkpython
# Script scaneps by Misha
#
# Calculates optical spectrum, \epsilon(\omega), from the provided nc-file
# with the bandstructure (TODO) and all
# and outputs nc-files with the susceptibility
# and dat-files with \epsilon(\omega) and n(\omega)

from NanoLanguage import *

import numpy
import cmath

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--upmin', default = 1)
parser.add_argument('--upmax', default = 13)
parser.add_argument('--downmin', default = 1)
parser.add_argument('--downmax', default = 13)
parser.add_argument('--fixed', action='store_true', help = "calculate spectrum for given --upmin and --donwmin")

args = parser.parse_args()

# We extract the name of file with the DFT results
try :
	with open('main.txt') as file:
		bandFileName =  file.readline()
except IOError as e:
	print "Unable to open file 'main.txt' with the reference to DFT calculations"
	raise

bandFileName_prefix = bandFileName[:(bandFileName.find("nc") - 1)]

specFilePrefix = "spec"
specFilePostfix = bandFileName

# nc-outputs
specFileName = specFilePrefix + "_" + specFilePostfix

# The results of calculations will be written to
epsPostfix = bandFileName_prefix + ".dat"

# First we read the configuration
bulk_configuration = nlread(bandFileName, BulkConfiguration)[0] 

# Next we calculate the optical spectrum

BandUpMin = int(args.upmin)

BandUpMax = BandUpMin if args.fixed else int(args.upmax)

BandDownMin = int(args.downmin)
BandDownMax = BandDownMin if args.fixed else int(args.downmax)

for down in range(BandDownMin, BandDownMax + 1):
	for up in range(BandUpMin, BandUpMax + 1):

		opt_spec = OpticalSpectrum(
    		configuration=bulk_configuration,
    		kpoints=MonkhorstPackGrid(15,15,15),
    		energies=numpy.linspace(0.1,2.5,200)*eV,
    		broadening=0.02*eV,
    		bands_below_fermi_level=down,
		    bands_above_fermi_level=up,
    		)
		specFileName = specFilePrefix +"_%s_%s_" % (down, up)
		specFileName += specFilePostfix
		nlsave(specFileName, opt_spec)


		# the array of energies
		enes = numpy.array(opt_spec.energies().tolist())

		# number of energy points
		nump = enes.shape[0]

		# optical parameters
		epsre = opt_spec.evaluateDielectricConstant()
		epsim = opt_spec.evaluateImaginaryDielectricConstant()

		for ii in range(3):
			for jj in range(ii, 3):

				epsFileName = "epsRE_%s_%s_%s_%s_" % (down, up, ii, jj) + epsPostfix
				with open(epsFileName, "w") as epsfile :
					for kk in range(nump) :
						epsfile.write(", ".join(str(x) 
							for x in (enes[kk], epsre[ii, jj, kk], epsim[ii, jj, kk])))
						epsfile.write("\n")
