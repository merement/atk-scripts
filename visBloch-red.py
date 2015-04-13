#!/usr/bin/env atkpython
# script visBloch by Misha
#
# Looks at all Bloch states in the respective nc-file and plots
# the spatial distribution of the density in the elementary cell.
#
# The plots are saved to png files.
#
# [2015-2-23] adopted for section files
#
# [2015-3-11] plots only selected spins
#
# TODO: define sections to plot, flag and options for the transversal
# [2015-3-24]
# Implemented the above TODO except options for the transversal.
#
# TODO: choice between bright and dark color schemes
#
# [2014-4-10]
# Implemented color choice
# TODO: options for the distribution along the z-axis

Tolerance = 1e-10

from NanoLanguage import *

import numpy

# to avoid keeping X-tunnel open
from matplotlib import use
use('Agg')

import matplotlib.pyplot as plt

# Common part for involved files
readmainFlag = True

import argparse
import json
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-s', default = 'All', help = 'Projection of the spin')
parser.add_argument('-c', default = 'NULL', help = 'The name of the .nc-file with the configuration')
parser.add_argument('blochFileName', help = 'The name of the .nc-file with Bloch functions')
parser.add_argument("--section", default = '[-1]', help="The list of sections to plot. \
    Expected format is --section [s1, s2, ...]. Parameters 0 < s_k < 1 specify the \
    relative z-coordinate of the (x, y)-plane. If s_k is outside of this range\
    the distribution integrated over z is plotted. ")
parser.add_argument("-t", action='store_true', help = "Provide to plot the distribution in the z-direction.")

# Arguments dealing with color schemes
color_select = parser.add_mutually_exclusive_group()
color_select.add_argument('--bright', action = 'store_true', 
        help = 'Provide to select "grey" scheme suitable for bw-printing')
color_select.add_argument('--dark', default = True, action = 'store_true', 
        help = 'Provide to select dark scheme (jet)')
color_select.add_argument('--color',  
        help = 'Specify the color scheme.')

args = parser.parse_args()

# validation of sections
try:
    listSections = json.loads(args.section)
except ValueError:
    print "The list of sections couldn't be recognized. \
    Expected format is --section [s1, s2, ...]"
    sys.exit(1)

# validation of the color choice
# http://matplotlib.org/examples/color/colormaps_reference.html

listcmaps = ['Blues', 'BuGn', 'BuPu','GnBu', 'Greens', 'Greys', 'Oranges', 
    'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 
    'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'bone', 'cool', 
    'copper', 'gist_heat', 'gray', 'hot', 'pink', 'spring', 'summer', 
    'winter', 'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
    'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'seismic', 'Accent', 
    'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3',
    'gist_earth', 'terrain', 'ocean', 'gist_stern', 'brg', 'CMRmap', 
    'cubehelix', 'gnuplot', 'gnuplot2', 'gist_ncar', 'nipy_spectral', 
    'jet', 'rainbow', 'gist_rainbow', 'hsv', 'flag', 'prism']

if args.bright :
    colChoice = 'gray'
elif args.dark :
    colChoice = 'jet'
elif args.color in listcmaps :
    colChoice = args.color
else:
    colChoice = 'jet' 

spinState = Spin.Up if args.s == 'Up' else \
    Spin.Down if args.s == 'Down' else Spin.All
    
strSpin = 'Up' if args.s == 'Up' else \
    'Down' if args.s == 'Down' else 'All'

if args.c == 'NULL' :
    # default action
    if readmainFlag :
    # We extract the name of file with the DFT results
        try :
            with open('main.txt') as file:
                bandFileName =  file.readline()
        except IOError as e:
            print "Unable to open file 'main.txt' with the reference to DFT calculations"
            raise
    else :
        # Or use the given file name
        bandFileName = "Mo_2_BS_DOS_OS_QAD.nc"
else :
    bandFileName = args.c

FileName_prefix = args.blochFileName[:(args.blochFileName.find(".nc"))]

def planarDensity(bState, ind) :
    """
    Calculates the density as the function of the in-plane coordinates.

    Returns 2D numpy array, and x and y coordinates of the grid points
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

def crossDensity(bState) :
    """
    Calculates the density as the function of the out-of-plane coordinate.

    Returns 
        two 1D numpy arrays - with \rho(z) and \psi'(z), where
                \psi'(z) = int dx dy real(psi(x,y,z))
        1D numpy array with z coordinates (in Angstroms)

    If the density is too small, the function returns zero, None, None

    If 0 < ind < 1, it's taken as the relative z-coordinate and the function
    outputs the distribution of the density at this plane. Otherwise,
    it calculates the integral over the z axis.
    """
    ni, nj, nk = bState.shape()
    # Find the volume elements.
    dX, dY, dZ = bState.volumeElement()
    # Calculate the unit area in the y-z plane.
    length_unit = dX.unit()
    dAXY = numpy.linalg.norm( numpy.cross(dX,dY) )*length_unit**2

    dz = dZ.norm()

    psi = bState.toArray()
    density = (psi*psi.conj()).real
    reals = psi.real

    norm = density.sum()
    if norm < Tolerance :
        return 0, None, None

    density = density/norm
#    reals = reals/3.0

    rhoz = numpy.array(
        [ density[:,:,ii].sum() for ii in range(nk) ])
    psiz = numpy.array(
        [ reals[:,:,ii].sum() for ii in range(nk) ])

    psiz = psiz/numpy.sqrt((psiz*psiz).sum())/3.0
    z_cell = numpy.array([bState.gridCoordinate(0, 0, ii).inUnitsOf(Angstrom)[2] 
                           for ii in range(nk)])

    return rhoz, psiz, z_cell

def show_density(state, lat, colors, ind) :
    """
    INPUT:
    	state - the instance of the BlochState class
    			describes the state we want to visualize
    	lat - the array of atom's coordinates
    	colors - the array of colors to depict atoms
    	ind - specifies which density we want to look at (see above)

    OUTPUT:
    	fig, ax - figure and axes handlers
    """
    print "Prosessing shift: ", ind
    n_xy, x, y = planarDensity(state, ind)

    if x == None :
        return None, None

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = fig.add_subplot(111, aspect = 'equal')
    #ax.plot_surface(x, y, n_xy)
    z_min, z_max = 0, numpy.abs(n_xy).max()
    #ax.pcolor(x, y, n_xy, cmap = colChoice, vmin = -z_max, vmax = z_max)
    ax.pcolor(x, y, -n_xy, cmap = colChoice)
    ax.scatter(lat[:,0], lat[:,1], c = colors, edgecolor = colors, marker = 'o', s = 3.0)
    ax.autoscale(tight=True)

    return fig, ax

def show_cross(state) :
    """
    INPUT:
        state - the instance of the BlochState class
                describes the state we want to visualize

    OUTPUT:
        fig, ax - figure and axes handlers
    """
    print "Prosessing cross densities... "
    rhoz, psiz, z = crossDensity(state)

    if z == None :
        return None, None

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = fig.add_subplot(111)
    ax.plot(z, rhoz, 'k', z, psiz, 'b')
    ax.autoscale(tight=True)

    return fig, ax

def show_density_fig(state, coords, colors, ind) :
    # These wrappers return only figure handlers
    fig, ax = show_density(state, coords, colors, ind)
    return fig

def show_cross_fig(state) :
    fig, ax = show_cross(state)
    return fig

def save_fig_list(fig_list, identity) :
    for ii in range(len(fig_list)) :
        if fig_list[ii] == None :
            continue
        print "Saving figure " + str(ii)
        outFileName = 'fig_' + str(ii) + '_' + identity + '.png'
        fig_list[ii].savefig(outFileName, bbox_inches='tight')

def buildColors(listE, colScheme) :
    """
    Creates the array with colors corresponding to different atoms.
    INPUT:
    	listE - array with names of the atoms
        colScheme - chosen colormap

    OUTPUT:
    	listC - array with colors: 'k' or 'w'
    """
    ## Ought to be done using map but this will do for now
    cDict = {
        'Sulfur' : 'w', 
        'Molybdenum' : 'k', 
        'Other' : 'r'}

    if colScheme == 'gray' :
        # for 'print' color scheme we want gray Sulfurs
        cDict['Sulfur'] = '0.75'
    # The key 'Other' is enacted when the element is absent in the 'database'

    listC = [cDict[elem] if elem in cDict else cDict['Other'] 
             for elem in listE ]

    return listC

# We need the configuration in order to get the coordinates of atoms
conf = nlread(bandFileName, BulkConfiguration, read_state=False)[0]  
coords = conf.cartesianCoordinates().inUnitsOf(Angstrom)
elems = conf.elements()
listElems = [elems[i].name() for i in range(len(elems))]
colorAtom = buildColors(listElems, colChoice)

# This requires a lot of memory
blochStates = nlread(args.blochFileName, BlochState)

numStates = len(blochStates)

# DEBUG version
DEB_flag = False
if DEB_flag :
    print "Debug version is used: list of states is truncated"
    numStates = 1

for i in range(numStates) :
    band = blochStates[i].quantumNumber()
    kPoint = blochStates[i].kPoint()
    print "Processing state %s (band # %s) out of %s ... " % (i, band, numStates)

    figs = [show_density_fig(blochStates[i].spinProjection(spinState), 
            coords, colorAtom, ind) for ind in listSections]
    if args.t :
        figs.append(show_cross_fig(blochStates[i].spinProjection(spinState)))
        
    ident = str(band) + "_" + strSpin +"_" + FileName_prefix
    save_fig_list(figs, ident)

