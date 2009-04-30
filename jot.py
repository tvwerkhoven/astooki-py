
import pyana
import libsh
import Clibshifts

pdir = './data/2009.04.27/proc/'
ddir = './data/2009.04.27/'
imfile = 'wfwfs_test_im27Apr2009.0001000'

# load SA / SF config
(nsa, saccdpos, saccdsize) = libsh.loadSaSfConf(pdir+'2009.04.22-mask.csv')
(nsf, sfccdpos, sfccdsize) = \
 	libsh.loadSaSfConf(pdir+'2009.04.22-subfield-big.csv')

# load image 
img = pyana.load(ddir+imfile)

# Pass to C library

Clibshifts.calcShifts(img, saccdpos, saccdsize, sfccdpos, sfccdsize, N.array([4,4]))

# done