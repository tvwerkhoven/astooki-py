#!/usr/bin/env python2.5
# encoding: utf-8

### Test libshifts-c library

import pyana
import libsh
import numpy
import scipy
#import Clibshifts

pdir = './data/2009.04.27/proc/'
ddir = './data/2009.04.27/'
imfile = 'wfwfs_test_im27Apr2009.0001000'

# load SA / SF config
(nsa, saccdpos, saccdsize) = libsh.loadSaSfConf(pdir+'2009.04.22-mask.csv')
(nsf, sfccdpos, sfccdsize) = \
 	libsh.loadSaSfConf(pdir+'2009.04.22-subfield-big.csv')

# load image 
img = pyana.getdata(ddir+imfile)

# Pass to C library

#Clibshifts.calcShifts(img, saccdpos, saccdsize, sfccdpos, sfccdsize, N.array([4,4]))

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