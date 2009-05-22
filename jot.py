#!/usr/bin/env python2.5
# encoding: utf-8

### Test libshifts-c library
import sys
import pyana
import numpy as N
import scipy
import astooki.libsh as libsh
import astooki.libshifts as libshifts
import astooki.clibshifts as clibshifts
import time

def main():
	pdir = 'data/2009.04.28-run01/proc/'
	ddir = 'data/2009.04.28-run01/'
	imfile = 'wfwfs_test_im28Apr2009.0000100'
	
	# load SA / SF config
	(nsa, saccdpos, saccdsize) = libsh.loadSaSfConf(pdir+'2009.04.28-mask.csv')
	(nsf, sfccdpos, sfccdsize) = \
		libsh.loadSaSfConf(pdir+'2009.04.28-subfield-big.csv')
	
	#saccdpos = saccdpos[20:30]
	nsa = len(saccdpos)
	# load image 
	img = pyana.getdata(ddir+imfile)
	img = img.astype(N.float32)
	
	# # Make fake image with gaussian
	# im1 = mk2dgauss((512,512), (153.5, 150.6))
	# im2 = mk2dgauss((512,512), (300, 150.))
	# im3 = mk2dgauss((512,512), (151.2, 303.3))
	# im4 = mk2dgauss((512,512), (302.4, 300.3))
	# img = (im1+im2+im3+im4).astype(N.float32)
	# saccdpos = N.array([[150, 150], [300, 150], [150, 300], [300, 300]], dtype=N.int32) - N.int32([32,32])
	# saccdsize = N.array([64,64], dtype=N.int32)
	# sfccdpos = N.array([[40, 40]], dtype=N.int32)
	# sfccdsize = N.array([12, 12], dtype=N.int32)
	# #saccdpos = saccdpos[:2]
	# nsa = len(saccdpos)
	# Pass to C library
	beg1 = time.time()
	datr = _libshifts.calcShifts(img, saccdpos, saccdsize, sfccdpos, sfccdsize, N.array([7,7]))
	print datr.keys()
	print datr['shifts'].shape
	# for sa in xrange(nsa):
	# 	print datr['shifts'][0][sa][0]
	print datr['refapts'].shape
	
	refsa = datr['refapts'][0]
	ref = img[saccdpos[refsa][1]:saccdpos[refsa][1]+saccdsize[1], saccdpos[refsa][0]:saccdpos[refsa][0]+saccdsize[0]].astype(N.float32)
	ref = ref/N.float32(ref.mean())
	pyshift = []
	beg2 = time.time()
	for sa in xrange(nsa):
		pos = saccdpos[sa]
		c = img[pos[1]:pos[1]+saccdsize[1], pos[0]:pos[0]+saccdsize[0]]
		
		#rms = (N.sum((c - c.mean())**2.0)/(N.product(saccdsize)))**0.5
		#print "sa #%d @ (%d,%d) mean: %g, rms: %g, %g" % (sa, pos[0], pos[1], c.mean(), rms, 100*rms/c.mean()),
		c = c / c.mean()
		_subfield = c[sfccdpos[0][1]:sfccdpos[0][1]+sfccdsize[1], \
			sfccdpos[0][0]:sfccdpos[0][0]+sfccdsize[0]].astype(N.float32)
		diffmap = libshifts.sqDiffWeave(_subfield, ref, sfccdpos[0], N.array([7,7]))
		#print "sf #%d mean: %g" % (0, _subfield.mean()),
		#print diffmap.max(), diffmap.min(), diffmap.mean()
		diffmap = diffmap.astype(N.float32)
		pyshift.append(libshifts.quadInt2dWeave(diffmap, range=N.array([7,7]), limit=N.array([7,7])))
	pyshift = N.array(pyshift)
	
	end = time.time()
	for sa in xrange(nsa):
		print "diff @ sa %d:" % (sa), datr['shifts'][0][sa][0] - pyshift[sa]
	
	print "C took: %g sec, Python took %g." % (beg2-beg1, end-beg2)
	


def mk2dgauss(size, orig):
	im = N.indices(size)
	dat1 = N.exp(-(im[0]-orig[1])**2.0/300.)
	dat2 = N.exp(-(im[1]-orig[0])**2.0/300.)
	dat = (dat1*dat2)/(dat1*dat2).max()
	return dat


	
def fill(l=None, m='e'):
	if l is not None:
		if m == 'e':
			print 'extending'
			l.extend([1,2,3,4])
		elif m == 'a':
			print 'appending'
			l.append([1,2,3,4])


if __name__ == "__main__":
	sys.exit(main())

pdir = '../data/2009.04.28-run01/proc/'
ddir = '../data/2009.04.28-run01/'
imfile = 'wfwfs_test_im28Apr2009.0000100'

# load SA / SF config
(nsa, saccdpos, saccdsize) = libsh.loadSaSfConf(pdir+'2009.04.28-mask.csv')
(nsf, sfccdpos, sfccdsize) = \
 	libsh.loadSaSfConf(pdir+'2009.04.28-subfield-big.csv')

# load image 
img = pyana.getdata(ddir+imfile)
img = img.astype(N.float32)

# Pass to C library

_libshifts.calcShifts(img, saccdpos, saccdsize, sfccdpos, sfccdsize, N.array([4,4]))

# done

### Process subshifts

import numpy as N

shfile = 'subshift/2009.04.28-run01_wfwfs_test_im28Apr2009_3-1002-shifts.npy'
subsh = N.load(shfile)

# Find non-finite entries
notfin = N.argwhere(N.isfinite(subsh) == False)
notfin_perc = notfin.shape[0]*100./subsh.size

# Find bad frames
# finframes = range(subsh.shape[0])
# for i in (N.unique(notfin[:,0]))[::-1]: finframes.pop(i)
# subshfin = subsh[finframes]
# NB: Almost all frames are bad...

# Set non-finite to 0 (bad solution)
notfin = N.argwhere(N.isfinite(subsh) == False)
for nfidx in notfin:
	subsh[tuple(nfidx)] = 0.0

# Subtract mean from all subfield shifts per frame and per reference subapt
for f in xrange(subsh.shape[0]):
	for r in xrange(subsh.shape[1]):
		avg = N.mean(subsh[f, r, :, :, :].reshape(-1,2), axis=0)
		subsh[f, r, :, :, :] -= avg.reshape(1,1,2)
		


# Compare different references
sfs = [1,50,51,61,71,100,150,269]
sas = [1,10,15,50,60,70,71]
for sa in sas:
	for sf in sfs:
		for i in xrange(3):
			print "%d,%d: (%.3g, %.3g) +- (%.3g, %.3g)" % \
				((sa, sf) + tuple(N.mean(subsh[:,i,sa,sf,:], 0)) + \
				tuple(N.var(subsh[:,i,sa,sf,:], 0)))


# Process shifts, average over frames and references

for r in xrange(subsh.shape[1]):
	subsh[:, r, :, :, :]
	avg = N.mean(subsh[f, r, :, :, :].reshape(-1,2), axis=0)
	subsh[f, r, :, :, :] -= avg.reshape(1,1,2)

### ==========================================================================
### Analyze shift data
### ==========================================================================


import numpy as N
import pylab

newsh = N.load('2009.04.28-run05_wfwfs_test_im28Apr2009_1-999-shifts.npy')
oldsh = N.load('2009.04.28-run05_wfwfs_test_im28Apr2009_1-999-shifts.npy.old0')

print newsh.shape, oldsh.shape

diff = oldsh[:,0,40,0] - newsh[:,0,40,0]
pylab.figure()
pylab.plot(diff[0], diff[1], 'o')

pylab.figure()
sa = 42
pylab.plot(newsh[:,0,sa,0,0], newsh[:,0,sa,0,1], 'o')
pylab.plot(newsh[:,1,sa,0,0], newsh[:,1,sa,0,1], 'o')
pylab.plot(newsh[:,2,sa,0,0], newsh[:,2,sa,0,1], 'o')
pylab.plot(newsh[:,3,sa,0,0], newsh[:,3,sa,0,1], 'o')

pylab.figure()
sa = 42
r1 = 0
r2 = 3
diffx = (newsh[:,r1,sa,0,0] - N.mean(newsh[:,r1,sa,0,0])) - \
	(newsh[:,r2,sa,0,0] - N.mean(newsh[:,r2,sa,0,0]))
diffy = (newsh[:,r1,sa,0,1] - N.mean(newsh[:,r1,sa,0,1])) - \
	(newsh[:,r2,sa,0,1] - N.mean(newsh[:,r2,sa,0,1]))
pylab.plot(diffx, diffy, 'o')

# Now average over the different references
nref = newsh.shape[1]
for i in range(nref):
	newsh[:,i,:,0,:] -= newsh[:,i,:,0,:].mean(1).reshape(-1,1,2)
newsh_a = newsh.mean(1)
print newsh_a.shape

pylab.figure()
sa = 42
pylab.plot(newsh_a[:,sa,0,0], newsh_a[:,sa,0,1], 'o')


### ==========================================================================
### Inspect tomographic inversion data
### ==========================================================================

# Settings
# ==================================

import numpy as N
import astooki.libfile as lf
import astooki.liblog as log
log.VERBOSITY +=2
import astooki.libsh as libsh
import astooki.libtomo as lt
import pylab

def rms(data, remdc=False):
if remdc:
return N.sqrt(N.mean((data-N.mean(data))**2.0))
else:
return N.sqrt(N.mean(data**2.0))


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
shifts = lf.loadData(shfile, asnpy=True)
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

# Allocate data for reconstruction
#recatm = N.zeros((shifts.shape[0], N.product(nlay), \
# 	len(nlay), lcells[0], lcells[1], 2), dtype=N.float32)
#inrms = N.zeros((shifts.shape[0], 4))
#recrms = N.zeros((shifts.shape[0], N.product(nlay), 4))
#diffrms = N.zeros((shifts.shape[0], N.product(nlay), 4))

# Process data
# ==================================

svdCache = lt.cacheSvd(geoms, lsizes, lorigs, lcells, sasize, sapos, sfang, \
 	sffov, matroot='/Users/tim/workdocs/data/matrices/')

# notfin = N.argwhere(N.isfinite(shifts) == False)
# if (notfin.shape[0] > 0):
# 	log.prNot(log.WARNING, "Some measurements are non-finite, check configuration!")			
# 	# TODO: setting to zero is a poor solution to NaNs
# 	for nfidx in notfin:
# 		shifts[tuple(nfidx)] = 0.0

# Setup inversion and forward matrices from SVD components
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
# APPLY TIP-TILT CORRECTION (DO NOT USE IN REAL ANALYSIS)
sh -= sh.mean(1).reshape(-1, 1, 2)
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

# Plot data
pylab.figure()
pylab.plot(sh[:,:,0].flatten(), sh[:,:,1].flatten(), 'x')
pylab.plot(sh[80,:,0].flatten(), sh[80,:,1].flatten(), 'x')
pylab.plot(sh[30,:,0].flatten(), sh[30,:,1].flatten(), 'x')
pylab.plot(sh[81,:,0].flatten(), sh[81,:,1].flatten(), 'x')

pylab.figure()
pylab.plot(wfwfsX, wfwfsY, 'o')
pylab.plot(wfwfsRec[0,0], wfwfsRec[0,1], 'o')
pylab.plot(wfwfsRec[1,0], wfwfsRec[1,1], 'o')
pylab.plot(wfwfsRec[2,0], wfwfsRec[2,1], 'o')
pylab.plot(wfwfsRec[3,0], wfwfsRec[3,1], 'o')

pylab.plot(wfwfsX-wfwfsRec[0,0], wfwfsY-wfwfsRec[0,1], 'o')

pylab.figure()
pylab.plot(wfwfsX)
pylab.plot(wfwfsRec[0,0]-wfwfsX)
pylab.plot(wfwfsRec[1,0]-wfwfsX)
pylab.plot(wfwfsRec[2,0]-wfwfsX)
pylab.plot(wfwfsRec[3,0]-wfwfsX)

pylab.figure()
pylab.plot(wfwfsX)
pylab.plot(wfwfsRec[0,0])
pylab.plot(wfwfsRec[5,0])
pylab.plot(wfwfsRec[9,0])

pylab.plot(wfwfsX-wfwfsRec[1,0])
pylab.plot(wfwfsX-wfwfsRec[2,0])
pylab.plot(wfwfsX-wfwfsRec[3,0])
pylab.plot(wfwfsX-wfwfsRec[4,0])
pylab.plot(wfwfsX-wfwfsRec[5,0])

pylab.figure()
pylab.plot(recatm[0,0].flatten())
pylab.plot(recatm[5,0].flatten())
pylab.plot(recatm[9,0].flatten())

### ==========================================================================
### Simulate tomographic analysis
### ==========================================================================

import numpy as N
import astooki.libfile as lf
import astooki.liblog as log
log.VERBOSITY +=2
import astooki.libsh as libsh
import astooki.libtomo as lt
import pylab

def rms(data, remdc=False):
if remdc:
return N.sqrt(N.mean((data-N.mean(data))**2.0))
else:
return N.sqrt(N.mean(data**2.0))



# Files
ddir = '/Users/tim/workdocs/data/wfwfs/sst/2009.04.28-run05/proc/'
safile = ddir + 'samask/2009.04.28-run05-samask-ll-centroid.csv'
sffile = ddir + 'sfmask/2009.04.28-run05-sfmask-20x20.csv'
# Settings
aptr = 0.49
ccdres = 0.35 * N.pi /60./60./180.
nlay = N.array([1, 2])
lh = N.array([[0, 0], [5000, 10000]])
lcells = N.array([8,8])

# Load data
(nsa, sapos, sasize) = libsh.loadSaSfConf(safile)
(nsf, sfccdpos, sfsize) = libsh.loadSaSfConf(sffile)

# Setup configuration
geoms = N.zeros((N.product(nlay), len(nlay)))
fov = (N.max(sfccdpos + sfsize, axis=0) - N.min(sfccdpos, axis=0)) * ccdres 
sffov = sfsize * ccdres
sfang = (sfccdpos + sfsize/2.0) * ccdres
sfang -= N.mean(sfang, axis=0)


for l in range(len(nlay)):
	thisHeights = N.linspace(lh[l,0], lh[l,1], nlay[l])
	geoms[:, l] = N.tile(N.repeat(thisHeights, N.product(nlay[l+1:])), N.product(nlay[:l]))


# Layer origin and sizes
telang = [0.0, 0.0]
lorigs = geoms.reshape(N.product(nlay), \
	 	len(nlay), 1) * N.tan(telang).reshape(1,1,2)
lsizes = aptr + \
		geoms.reshape(N.product(nlay), len(nlay), 1)*\
		N.tan(0.5 * fov).reshape(1,1,2)

# Generate one forward matrix first
# (mat, matc, mattag) = lt.computeFwdMatrix(geoms[1], lsizes[1], lorigs[1], lcells, sasize, sapos, sfang, sffov, matroot='/Users/tim/workdocs/data/matrices/')
# matd = mat-matc

# Generate whole SVD cache
svdCache = lt.cacheSvd(geoms, lsizes, lorigs, lcells, sasize, sapos, sfang, \
 	sffov, matroot='/Users/tim/workdocs/data/matrices/')


# Setup inversion and forward matrices from SVD components
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

fakeatm = []
fakedata = []
fakerecatm = []
fakerecdata = []

for geom in range(N.product(nlay)):
	fakeatm.append(N.random.random((len(nlay), lcells[0], lcells[1])))
	print fakeatm[-1].shape
	fakedata.append(N.dot(modmats['fwd'][geom], fakeatm[-1].flatten()))
	print fakedata[-1].shape
	fakerecatm.append(N.dot(modmats['inv'][geom], fakedata[-1]).reshape( \
		len(nlay), lcells[0], lcells[1]))
	print fakerecatm[-1].shape
	fakerecdata.append(N.dot(modmats['fwd'][geom], fakerecatm[-1].flatten()))

fakeatm = N.array(fakeatm)
fakedata = N.array(fakedata)
fakerecatm = N.array(fakerecatm)
fakerecdata = N.array(fakerecdata)

# Plot
atm = 5
pylab.figure()
pylab.plot(fakeatm[atm].flatten())
pylab.plot(fakerecatm[atm].flatten())
pylab.plot((fakeatm[atm] - fakerecatm[atm]).flatten())

pylab.figure()
pylab.plot(fakedata[atm].flatten())
pylab.plot(fakerecdata[atm].flatten())
pylab.plot((fakedata[atm] - fakerecdata[atm]).flatten())


### ==========================================================================
### Analyze tomographic inversion data
### ==========================================================================

import astooki.libfile as lf
import pylab

dat, meta = lf.restoreData('astooki-meta-data.pickle')

inrms = dat['inrms'][:750]
diffrms = dat['diffrms'][:750]
recrms = dat['recrms'][:750]

pylab.figure()
pylab.plot(inrms[:,0])
pylab.plot(recrms[:,:,0].mean(1))
pylab.plot(diffrms[:,:,0].mean(1))

pylab.figure()
pylab.plot(inrms[:,0], inrms[:,1], 'o')
pylab.plot(recrms[:,:,0].mean(1), recrms[:,:,1].mean(1), 'o')
pylab.plot(diffrms[:,:,0].mean(1), diffrms[:,:,1].mean(1), 'o')

# Average over all frames, exposing the inter-geometry differences
diffrms_a = diffrms.mean(0)
diffrms_s = diffrms.std(0)
recrms_a = recrms.mean(0)
recrms_s = recrms.std(0)

pylab.figure()
pylab.plot(diffrms_a[:,0])
pylab.plot(diffrms_a[:,1])
# pylab.errorbar(range(diffrms.shape[1]), diffrms_a[:,0], yerr=diffrms_s[:,0], fmt='o')
# pylab.errorbar(range(diffrms.shape[1]), diffrms_a[:,1], yerr=diffrms_s[:,1], fmt='o')

pylab.figure()
pylab.plot(recrms_a[:,0])
pylab.plot(recrms_a[:,1])
# pylab.errorbar(range(recrms.shape[1]), recrms_a[:,0], yerr=recrms_a[:,0], fmt='o')
# pylab.errorbar(range(recrms.shape[1]), recrms_a[:,1], yerr=recrms_a[:,1], fmt='o')

pylab.cla()
pylab.plot(diffrms_a[:,0])
pylab.plot(recrms_a[:,0])

### ==========================================================================
### Subfield mask problems
### ==========================================================================
import numpy as N

sasize = N.array([88,92])
border = N.array([10,10])
sfsize = N.array([68,72])
overlap = N.array([0,0])


sasize = N.array([88,92])
border = N.array([6,6])
sfsize = N.array([16,16])
overlap = N.array([0.5,0.5])

effsize = sasize - 2*border
print effsize
pitch = sfsize * (1-overlap)
print pitch
nsf = N.floor((effsize-sfsize+pitch) / pitch)
print nsf
effpitch = (effsize-sfsize)/(nsf-1)
effpitch[(nsf == 1)] = 0
print effpitch

sfpos = border + \
	N.indices(nsf, dtype=N.float).reshape(2,-1).T * effpitch
sfpos = N.floor(sfpos).astype(N.int32)

print sasize
print sfpos.min(0)
print sfpos.max(0), sfpos.max(0) + sfsize + border