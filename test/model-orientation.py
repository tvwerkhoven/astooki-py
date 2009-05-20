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

* Results

- Without 'tip-tilt correction'

[ 1.  1.]
[ 2.56105151  2.64738488  2.56105151  2.64738488]
[array([ 2.41957  ,  2.5457506,  2.41957  ,  2.5457506]), 
array([ 2.41942144,  2.54575878,  2.41942144,  2.54575878])]
[array([ 0.83944372,  0.72649886,  0.83944372,  0.72649886]), 
array([ 0.83987182,  0.72647032,  0.83987182,  0.72647032])]

[ 1. -1.]
[ 2.56105151  2.64738488  2.56105151  2.64738488]
[array([ 2.41964949,  2.54604396,  2.41964948,  2.54604396]), 
array([ 2.41946316,  2.54603927,  2.41946316,  2.54603927])]
[array([ 0.83921454,  0.72547021,  0.83921454,  0.72547021]), 
array([ 0.8397516 ,  0.72548666,  0.8397516 ,  0.72548666])]

[-1.  1.]
[ 2.56105151  2.64738488  2.56105151  2.64738488]
[array([ 2.42000366,  2.54567882,  2.42000366,  2.54567882]), 
array([ 2.41970324,  2.5457816 ,  2.41970324,  2.5457816 ])]
[array([ 0.83819277,  0.72675039,  0.83819277,  0.72675038]), 
array([ 0.83905963,  0.72639029,  0.83905963,  0.72639029])]

[-1. -1.]
[ 2.56105151  2.64738488  2.56105151  2.64738488]
[array([ 2.42080883,  2.54691589,  2.42080883,  2.54691589]), 
array([ 2.42014299,  2.54663038,  2.42014299,  2.54663038])]
[array([ 0.8358644 ,  0.72240311,  0.8358644 ,  0.72240311]), 
array([ 0.83779048,  0.72340906,  0.83779048,  0.72340906])]

- With 'tip-tilt corretion'

[ 1.  1.]
[ 0.30162596  0.38384209  0.30162596  0.38384209]
[array([ 0.05644125,  0.04671362,  0.05644125,  0.04671362]), 
array([ 0.04917027,  0.04716309,  0.04917027,  0.04716309])]
[array([ 0.29629817,  0.38098896,  0.29629817,  0.38098896]), 
array([ 0.29759117,  0.38093358,  0.29759117,  0.38093358])]
0.0160326649117 8.2087244348 0.16192476782

[ 1. -1.]
[ 0.30162596  0.38384209  0.30162596  0.38384209]
[array([ 0.0591062 ,  0.06079126,  0.0591062 ,  0.06079126]), 
array([ 0.0499766 ,  0.06202045,  0.0499766 ,  0.06202045])]
[array([ 0.29577809,  0.37899759,  0.29577809,  0.37899759]), 
array([ 0.29745682,  0.37879838,  0.29745682,  0.37879838])]
0.00813757557978 4.16643869685 0.123269427078


[-1.  1.]
[ 0.30162596  0.38384209  0.30162596  0.38384209]
[array([ 0.07193514,  0.04386062,  0.07193514,  0.04386062]), 
array([ 0.06129477,  0.04826221,  0.06129477,  0.04826221])]
[array([ 0.29292244,  0.38132794,  0.29292244,  0.38132794]), 
array([ 0.29533231,  0.38079589,  0.29533231,  0.38079589])]
0.0128423154853 6.57526552845 0.155788323111

[-1. -1.]
[ 0.30162596  0.38384209  0.30162596  0.38384209]
[array([ 0.0945659 ,  0.08940066,  0.0945659 ,  0.08940066]), 
array([ 0.07630392,  0.08286977,  0.07630392,  0.08286977])]
[array([ 0.28641842,  0.37328578,  0.28641842,  0.37328578]), 
array([ 0.29181489,  0.37478975,  0.29181489,  0.37478974])]
-0.00619474457005 -3.17170921987 0.197274495195

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
nlay = N.array([1,1])
lh = N.array([[5000, 5000], [25000, 25000]])
lcells = N.array([10,10])

# Load data
_shifts = lf.loadData(shfile, asnpy=True)
# Process one shift set
it = 400
shifts = _shifts[it]
del _shifts
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
lsizes = (aptr + \
		geoms.reshape(N.product(nlay), len(nlay), 1)*\
		N.tan(0.5 * fov).reshape(1,1,2)) * 0.9

# Setup matrices
# ==================================

orientations = N.array([[1,-1], [1,1], [-1,1], [-1,-1]], dtype=N.float)
data = {}
for _or in orientations:
	orid = str(_or[0])+str(_or[1])
	data[orid] = {}
	svdCache = lt.cacheSvd(geoms, lsizes, lorigs, lcells, sasize, sapos, \
		sfang * _or, sffov, matroot='/Users/tim/workdocs/data/matrices/')
		
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


	# Average over reference
	sh = shifts.mean(0)
	# Set sa 81 and 84 to zero, bad subaps
	#sh[81] = 0
	#sh[84] = 0
	# Set average over all subaps to zero
	sh -= sh.mean(0).reshape(1, -1, 2)
	# SIMULATE TIP-TILT CORRECTION (DO NOT USE IN REAL ANALYSIS)
	sh -= sh.mean(1).reshape(-1, 1, 2)
	#log.prNot(log.NOTICE, "Data frame %d/%d" % (it+1, shifts.shape[0]))

	inrms = N.r_[\
		rms(sh[...,0]), \
		rms(sh[...,1]), \
		rms(sh[...,0], remdc=True), \
		rms(sh[...,1], remdc=True)]
	
	data[orid]['inrms'] = inrms
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
		wfwfsRec.append( \
			N.array([\
				N.dot(modmats['fwd'][geom], \
					recatm[-1][0].flatten()), \
				N.dot(modmats['fwd'][geom], \
					recatm[-1][1].flatten()) \
			]))
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

	data[orid]['recatm'] = N.array(recatm)
	data[orid]['wfwfsRec'] = N.array(wfwfsRec)
	data[orid]['recrms'] = N.array(recrms)
	data[orid]['diffrms'] = N.array(diffrms)
	
	#print _or
	#print inrms
	#print recrms
	#print diffrms
	#print recatm.mean(), recatm.sum(), recatm.std()

print "Done"

dk = data.keys()
print dk
print data[dk[0]].keys()
for k in data.keys():
	print "==============", k, "=============="
	for (quant, val) in data[k].items():
		if (quant in ['inrms', 'diffrms', 'recrms']):
			print quant, "----------------"
			print val
		if (quant == 'diffrms'):
			pylab.figure(1)
			#pylab.plot(val[:,2], val[:,3], '-')
			pylab.plot([val[0,2]], [val[0,3]], 'o')
			#pylab.plot([val[1,2]], [val[1,3]], 'x')
		elif (quant == 'recrms'):
			pylab.figure(2)
			#pylab.plot(val[:,2], val[:,3], '-')
			pylab.plot([val[0,2]], [val[0,3]], 'o')
			#pylab.plot([val[1,2]], [val[1,3]], 'x')
		elif (quant == 'recatm'):
			print 'fig3'
			val[0,0].flatten().shape
			pylab.figure(4)
			#pylab.plot(val[:,2], val[:,3], '-')
			pylab.plot(val[0,0].flatten())
			pylab.plot(val[0,1].flatten())
