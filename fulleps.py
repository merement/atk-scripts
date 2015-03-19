# Script fulleps by Misha
#
# Calculates optical spectrum, \epsilon(\omega), from the provided nc-file
# with the bandstructure and all
# and outputs dat-files with \epsilon(\omega) and n(\omega)
#
# TODO: level it up with the current ideology for such scripts
#	1) support of command line options
#	2) ...
#
# It should be replaced by a boosted version of scaneps.py

from NanoLanguage import *
#import NLEngine

import numpy
import cmath

# We extract the name of file with the DFT results
try :
	with open('main.txt') as file:
		bandFileName =  file.readline()
except IOError as e:
	print "Unable to open file 'main.txt' with the reference to DFT calculations"
	raise

bandFileName_prefix = bandFileName[:(bandFileName.find("nc") - 1)]
bandFileName = bandFileName_prefix + ".nc"
specFileName = "spec_" + bandFileName

# The results of calculations will be written to
epsxyFileName = "eps_xy_" + bandFileName_prefix + ".dat"
epszFileName = "eps_z_" + bandFileName_prefix + ".dat"
nFileName =  "n_" + bandFileName_prefix + ".dat"

# First we read the configuration
bulk_configuration = nlread(bandFileName, BulkConfiguration)[0] 
			# , object_id='gID003'

# Next we calculate the optical spectrum

opt_spec = OpticalSpectrum(
    configuration=bulk_configuration,
    kpoints=MonkhorstPackGrid(15,15,15),
    energies=numpy.linspace(0.1,4.5,900)*eV,
    broadening=0.025*eV,
    bands_below_fermi_level=13,
    bands_above_fermi_level=13,
    )
nlsave(specFileName, opt_spec)

# the array of energies
enes = numpy.array(opt_spec.energies().tolist())

# number of energy points
nump = enes.shape[0]

# The dielectric tensor is pretty much diagonal and its in-plane
# eigenvalues are the same

epsre = opt_spec.evaluateDielectricConstant()[0,0,:]
epsim = opt_spec.evaluateImaginaryDielectricConstant()[0,0,:]

ezre = opt_spec.evaluateDielectricConstant()[2,2,:]
ezim = opt_spec.evaluateImaginaryDielectricConstant()[2,2,:]


with open(epsxyFileName, "w") as epsfile :
	for ii in range(nump) :
		epsfile.write(", ".join(str(x) for x in (enes[ii], epsre[ii], epsim[ii])))
		epsfile.write("\n")

with open(epszFileName, "w") as epsfile :
	for ii in range(nump) :
		epsfile.write(", ".join(str(x) for x in (enes[ii], ezre[ii], ezim[ii])))
		epsfile.write("\n")

with open(nFileName, "w") as nfile :
	for ii in range(nump) :
		nind = cmath.sqrt(epsre[ii] + 1j*epsim[ii])
		nfile.write(", ".join(str(x) for x in (enes[ii], nind.real, nind.imag)))
		nfile.write("\n")
