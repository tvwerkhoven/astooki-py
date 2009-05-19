#!/usr/bin/env python2.5
# encoding: utf-8

### Test libshifts-c library
import sys
import pyana
import libsh
import numpy as N
import scipy
import _libshifts
import libshifts
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


### Analyze tomographic inversion data
### ==================================

import libfile as lf
import pylab

dat, meta = lf.restoreData('astooki-meta-data.pickle')

inrms = dat['inrms'][:750]
diffrms = dat['diffrms'][:750]
recrms = dat['recrms'][:750]

pylab.cla()
pylab.plot(inrms[:,0])
pylab.plot(recrms[:,:,0].mean(1))
pylab.plot(diffrms[:,:,0].mean(1))

# take only until frame 750
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
