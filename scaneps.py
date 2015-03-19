# Script scaneps by Misha
#
# Calculates optical spectrum, \epsilon(\omega), from the provided nc-file
# with the bandstructure and all
# and outputs dat-files with \epsilon(\omega) and n(\omega)

from NanoLanguage import *
#import NLEngine

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
bandFileName = bandFileName_prefix + ".nc"

specFilePrefix = "spec"
specFilePostfix = bandFileName

specFileName = "spec_" + bandFileName

# The results of calculations will be written to
epsPostfix = bandFileName_prefix + ".dat"

# First we read the configuration
bulk_configuration = nlread(bandFileName, BulkConfiguration)[0] 
			# , object_id='gID003'

# Next we calculate the optical spectrum

BandUpMin = int(args.upmin)
if args.fixed :
	BandUpMax = BandUpMin
else:
	BandUpMax = int(args.upmax)

BandDownMin = int(args.downmin)
if args.fixed :
	BandDownMax = BandDownMin
else:
	BandDownMax = int(args.downmax)

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

		# The dielectric tensor is pretty much diagonal and its in-plane
		# eigenvalues are the same

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
