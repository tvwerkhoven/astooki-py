#!/usr/bin/env python
# encoding: utf-8
"""
* Testing atmosphere model coordinate system, TvW, 20090520

The coordinate systems of the telescope and the model must match. The two
coordinate systems involved are the one of the lenslets (1) and the one of the
subfield angles (2).

Coordinate system 1 is a spatial coordinate system which is not model
dependent, if there is a rotation or mirroring in this system in the model wrt
the telescope, this does matter for the results.

Coordinate system 2 is an angular coordinate system, which defines the angles
of the subfields in the WFWFS images. This coordinate system *is* model
dependent. If the system is mirrored wrt the real system, this means that
we're trying to find correlations between beams that do not exist.

The four conversions that are possible are:

(2a) angle(x,y) = pos(x,y) * ccdres
(2b) angle(x,y) = pos(-x,y) * ccdres
(2c) angle(x,y) = pos(x,-y) * ccdres
(2d) angle(x,y) = pos(-x,-y) * ccdres

with ccdres the angular resolution of the ccd (arcsec/pixel).

For example, if we assume the relation between the pixel position and the
angle we look at is (2a), but in fact it is (2d), we project the subfield
beams on the wrong atmospheric layer cells such that two subfields that *do*
sense the same atmosphere are attributed to different cells.

This script attempts to find the correct setup.
"""

# Libraries
# ==================================

import numpy as N
import libfile as lf
import liblog as log
log.VERBOSITY +=2
import libsh
import libtomo as lt
import pylab

# Functions
# ==================================

def rms(data, remdc=False):
	if remdc:
		return N.sqrt(N.mean((data-N.mean(data))**2.0))
	else:
		return N.sqrt(N.mean(data**2.0))


# Settings
# ==================================

# Files
ddir = '/Users/tim/workdocs/data/wfwfs/sst/2009.04.28-run05/proc/'
shfile = ddir + 'subshift-20x20/2009.04.28-run05_wfwfs_test_im28Apr2009_1-999-shifts.npy'
safile = ddir + 'samask/2009.04.28-run05-samask-ll-centroid.csv'
sffile = ddir + 'sfmask/2009.04.28-run05-sfmask-20x20.csv'

# Settings
aptr = 0.49
ccdres = 0.35 * N.pi /60./60./180.
nlay = N.array([1, 2])
lh = N.array([[0,0], [5000, 10000]])
lcells = N.array([8,8])

# Load data
(nsa, sapos, sasize) = libsh.loadSaSfConf(safile)
(nsf, sfccdpos, sfsize) = libsh.loadSaSfConf(sffile)

# Setup configuration
geoms = N.zeros((N.product(nlay), len(nlay)))
origs = N.zeros((N.product(nlay), len(nlay)))
sizes = N.zeros((N.product(nlay), len(nlay)))
fov = (N.max(sfccdpos + sfsize, axis=0) - N.min(sfccdpos, axis=0)) * ccdres 
sffov = sfsize * ccdres
sfang = (sfccdpos + sfsize/2.0) * ccdres
sfang -= N.mean(sfang, axis=0)

# Setup all layer configurations
for l in range(len(nlay)):
	thisHeights = N.linspace(lh[l,0], lh[l,1], nlay[l])
	geoms[:, l] = N.tile(\
		N.repeat(thisHeights, N.product(nlay[l+1:])), N.product(nlay[:l]))

# Layer origin and sizes
telang = [0.0, 0.0]
lorigs = geoms.reshape(N.product(nlay), \
	 	len(nlay), 1) * N.tan(telang).reshape(1,1,2)
lsizes = aptr + \
		geoms.reshape(N.product(nlay), len(nlay), 1)*\
		N.tan(0.5 * fov).reshape(1,1,2)

# Setup matrices
# ==================================

svdCache = lt.cacheSvd(geoms, lsizes, lorigs, lcells, sasize, sapos, sfang, \
 	sffov, matroot='/Users/tim/workdocs/data/matrices/')

modmats = {}
modmats['inv'] = []
modmats['fwd'] = []
log.prNot(log.NOTICE, "Pre-computing inversion and forward matrices...")

for geom in range(N.product(nlay)):
	# Setup inversion model matrix
	modmats['inv'].append(\
		N.dot(\
			svdCache[geom]['vh'].T, \
			N.dot(\
				svdCache[geom]['s_inv'] *
					N.identity(len(svdCache[geom]['s_inv'])), \
					svdCache[geom]['u'].T)))
	# Setup forward model matrix
	modmats['fwd'].append(\
		N.dot(\
			svdCache[geom]['u'], \
			N.dot(\
				svdCache[geom]['s'] *\
				N.identity(len(svdCache[geom]['s'])), \
				svdCache[geom]['vh'])))


# Process data
# ==================================
shifts = lf.loadData(shfile, asnpy=True)

# Average over number of references
#shifts_a = shifts.mean(axis=1)
# Remove average over all frames
# shifts_aa = shifts_a - N.mean(shifts_a, axis=0).reshape(1, shifts_a.shape[1], shifts_a.shape[2], 2)

# Process one shift set
it = 400
sh = shifts[it]
# Average over reference
sh = sh.mean(0)
# Set sa 81 and 84 to zero, bad subaps
sh[81] = 0
sh[84] = 0
# Set average over all subaps to zero
sh -= sh.mean(0).reshape(1, -1, 2)
# SIMULATE TIP-TILT CORRECTION (DO NOT USE IN REAL ANALYSIS)
#sh -= sh.mean(1).reshape(-1, 1, 2)
log.prNot(log.NOTICE, "Data frame %d/%d" % (it+1, shifts.shape[0]))

inrms = N.r_[\
	rms(sh[...,0]), \
	rms(sh[...,1]), \
	rms(sh[...,0], remdc=True), \
	rms(sh[...,1], remdc=True)]

wfwfsX = sh[...,0].flatten()
wfwfsY = sh[...,1].flatten()

recatm = []
wfwfsRec = []
recrms = []
diffrms = []

for geom in range(N.product(nlay)):
	recatm.append( \
		N.array([\
			(N.dot(modmats['inv'][geom], wfwfsX)).reshape( \
				len(nlay), lcells[0], lcells[1]), \
			(N.dot(modmats['inv'][geom], wfwfsY)).reshape( \
						len(nlay), lcells[0], lcells[1])\
		]))
	print recatm[-1].shape
	wfwfsRec.append( \
		N.array([\
			N.dot(modmats['fwd'][geom], \
				recatm[-1][0].flatten()), \
			N.dot(modmats['fwd'][geom], \
				recatm[-1][1].flatten()) \
		]))
	print wfwfsRec[-1].shape
	recrms.append(N.r_[\
		rms(wfwfsRec[-1][0]), \
		rms(wfwfsRec[-1][1]), \
		rms(wfwfsRec[-1][0], remdc=True), \
		rms(wfwfsRec[-1][1], remdc=True)])
	diffrms.append(N.r_[\
		rms(wfwfsRec[-1][0]-wfwfsX),\
		rms(wfwfsRec[-1][1]-wfwfsY),\
		rms(wfwfsRec[-1][0]-wfwfsX,remdc=True),\
		rms(wfwfsRec[-1][1]-wfwfsY,remdc=True)])

wfwfsRec = N.array(wfwfsRec)
recatm = N.array(recatm)
