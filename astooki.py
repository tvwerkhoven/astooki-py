#!/usr/bin/env /sw/bin/python2.5
# encoding: utf-8
"""
@filename astooki.py
@author Tim van Werkhoven (tim@astro.su.se)
@date 20090422 14:34
@brief AsTooki: the astronomical toolkit - imagemagick for astronomers.

Created by Tim van Werkhoven on 2009-04-22.
Copyright (c) 2009 Tim van Werkhoven (tim@astro.su.se)

This file is licensed under the Creative Commons Attribution-Share Alike
license versions 3.0 or higher, see
http://creativecommons.org/licenses/by-sa/3.0/
"""

import sys
import getopt
import os
import numpy as N
import scipy as S
import libsh
import libfile
import libplot
import liblog as log

VERSION = "0.0.1"
AUTHOR = "Tim van Werkhoven (tim@astro.su.se)"
DATE = "2009-04-24"

help_message = '''astooki version %s (%s) by %s.
Usage: astooki <TOOL> [OPTIONS] [FILES]

Tools
 convert                     Convert files to another format
 stats                       Get statistics on files
 samask                      Make a subaperture mask
 sfmask                      Make a subfield mask
 saopt                       Optimze a subaperture mask with a flat
 saupd                       Update a subaperture mask with an offset
 shifts                      Measure image shifts in various subfields and 
                               subimages

Input formats supported
 ana                         ANA format
 fits                        Flexible Image Transport System
 npy                         NumPy binary data files

Output formats supported
 ana                         ANA format
 fits                        Flexible Image Transport System
 png                         Portable Network Graphics
 npy                         NumPy binary data files

Common options
 -v, --verbose               increase verbosity
 -h, --help                  show this help
 -s, --stats                 show RMS of each file processed [False]
 -i, --informat=FORMAT       input file-format [ana]
     --ff=FILEPATH           use this flatfield when correcting files [none]
     --fm=N                  flatfield consists of N summed frames
     --df=FILEPATH           use this darkfield when correcting files [none]
     --dm=N                  darkfield consists of N summed frames
     --mf=FILEPATH           use this configuration as mask file [none]
     --[no]norm              normalize pixels when masking [True]
     --crop=X,Y,W,H          crop region at X,Y size W,H

Convert options
 -f, --file=FILEPATH         if single file, store converted file here 
     --scale=FACTOR          scale output resolution by this factor [1.0]
 -o, --outformat=FORMAT      output file-format [fits]
     --intclip=LOW,HIGH      clip intensity to this range

Stats options
 -f, --file=FILEPATH         file to store statistics to

Samask options
 -f, --file=FILEPATH         file to store subaperture configuration to
     --rad=RAD               radius of the subaperture pattern [1024]
     --shape=SHAPE           pattern shape (circular or square) [circular]
     --size=X,Y              size of the subapertures
     --pitch=X,Y             pitch of the subapertures
     --xoff=EVENOFF,ODDOFF   x-offset for subapertures in even, odd rows in
                               units of --size [0, 0.5]
     --disp=X,Y              global pattern offset [0,0]
     --scale=SCALE           global pattern scaling factor [1.0]
     --[no]plot              make a plot of the subaperture mask [yes]

Sfmask options
 -f, --file=FILEPATH         file to store subaperture configuration to
     --sfsize=X,Y            subfield size [16, 16]
     --sasize=X,Y            subaperture size to fit things in
     --overlap=X,Y           how many overlap to allow between subfields in X 
                               and Y direction [0.5, 0.5]
     --border=X,Y            add a border around the subfields

Saopt options
 -f, --file=FILEPATH         file to store optimized configuration to
     --saifac=FLOAT          intensity drop-off factor considered 'dark' [0.7]
     --rad=RAD               radius of the subaperture pattern, used for 
                               plotting [1024]
Saupd options
     --offsets=FILE          file holding offset vectors for all subapertures

Shifts options
 -r, --range=INT             shift range to use for cross-correlation [7]
     --safile=FILEPATH       subaperture locations, same format as maskfile
     --sffile=FILEPATH       subfield locations w.r.t. subaperture locations
 -n, --nref=INT              number of references to use [4]

Examples
 To calculate stats for a series of files in one directory, using dark- and
 flatfielding and cropping using a pre-calculated mask:
   astooki.py stats -vvv --ff *ff*1001 --fm=500 --df ../2009.04.22/*dd*1995 \\
   --dm=500 --mf=proc/2009.04.22-mask-crop3.csv --file=proc/2009.04.27-stats \\  
   wfwfs_test_im27Apr2009.0000???

 To convert one file to fits format while dark/flat-fielding and using a 
 cropmask without normalizing:
   astooki.py convert -vvv --ff *ff*1001 --fm=500 --df \\
   ../2009.04.22/*dd*1995 --dm=500 --mf=proc/2009.04.22-mask-crop3.csv \\
   --nonorm wfwfs_test_im27Apr2009.0000042

 To crop a portion of the image and store it as png:
   astooki.py convert -vvv --ff *ff*1001 --fm=500 --df \\
   ../2009.04.22/*dd*1995 --dm=500 --crop 652,1108,168,146 --mf \\
   proc/2009.04.22-mask.csv --intclip 0.9,1.1 --outformat png \\
   wfwfs_test_im27Apr2009.0000001

 To measure image shifts over the whole subaperture, use a subfield file with 
 only one subfield with the size of the complete subimage except for a guard 
 range, combined with a regular subimage config file:
   astooki.py shifts -vvv --ff *ff*1001 --fm=500 --df \\
   ../2009.04.22/*dd*1995 --dm=500 --safile proc/2009.04.22-mask.csv \\
   --sffile proc/2009.04.22-subfield-big.csv --range 7 --nref 5 \\
   --file wfwfs_test_im27Apr2009-shifts wfwfs_test_im27Apr2009.0000{1,2}??

 Updating subaperture positions:
 astooki.py saupd -vv --plot --mf 2009.04.22-mask.csv --offset offset-csv.csv

 Making a subfield mask with 16x16 pix subfields and 0.5 overlap on each side:
 astooki.py sfmask -vv --plot --file 'subfieldmask.csv' --sfsize=16,16 \\
   --sasize=165,139 --overlap=0.5,0.5 --border=6,6

 Subfield shift measurement using the above corrected subimage positions and
 subfields:
 astooki.py shifts -vv --ff *ff*1001 --fm=500 --df ../2009.04.22/*dd*1995 \\
   --dm=500 --safile proc/2009.04.22-mask-updated.csv --sffile \\
   proc/2009.04.27-subfieldmask.csv --range 6 --nref 2 \\
   wfwfs_test_im27Apr2009.000020?

 Convert numpy data to fits format:
 astooki.py convert -vvv --informat npy --outformat fits \\
   2009.04.27_wfwfs_test_im27Apr2009_200-209-shifts.npy
''' % (VERSION, DATE, AUTHOR)

### Supported formats
_FORMAT_ANA = 'ana'
_FORMAT_FITS = 'fits'
_FORMAT_PNG = 'png'
_FORMAT_NPY = 'npy'
_INFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_NPY)
_OUTFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_PNG, _FORMAT_NPY)
# Tools available
_TOOLS = ('convert', 'stats', 'samask', 'sfmask', 'saopt', 'saupd', 'shifts')

### ==========================================================================
### Startup functions
### ==========================================================================

def main(argv=None):
	# Parse command-line options from sys.argv
	(tool, params, files) = parse_options()
	# Sanity check on parameters
	check_params(tool, params)
	# Perform action requested
	log.prNot(log.INFO, "Tool: %s." % tool)
	if (tool == 'convert'): ConvertTool(files,params)
	elif (tool == 'stats'): StatsTool(files, params)
	elif (tool == 'samask'): SubaptConfTool(files, params)
	elif (tool == 'sfmask'): SubfieldConfTool(files, params)
	elif (tool == 'saopt'): SubaptOptTool(files, params)
	elif (tool == 'saupd'): SubaptUpdateTool(files, params)
	elif (tool == 'shifts'): ShiftTool(files, params)
	# done
	log.prNot(log.INFO, "Complete.")
	return 0


def parse_options():
	"""
	Process command-line options. See help_message for possible options
	"""
	
	argv = sys.argv
	# First check whether argv[1] is present (could be -h, --help or a tool)
	try:
		tool = argv[1]
		if tool in ["-h", "--help"]: print_help()
		elif tool not in _TOOLS:
			raise Exception
	except:
		print >> sys.stderr, os.path.basename(sys.argv[0]) + ": Syntax incorrect."
		print >> sys.stderr, "\t for help use --help"
		sys.exit(2)
	
		
	
	# Parse common options first
	# ==========================
	
	params = get_defaults(tool)
	
	opts, args = getopt.getopt(argv[2:], "vhsi:o:f:r:n:", ["verbose", "help", "stats", "informat=", "ff=", "fm=", "df=", "dm=", "mf=", "outformat=", "intclip=", "crop=", "file=", "scale=", "rad=", "shape=", "size=", "pitch=", "xoff=", "disp=", "plot", "noplot", "norm", "nonorm", "saifac=", "range=", "sffile=", "safile=", "nref=", "sfsize=", "sasize=", "overlap=", "border=", "offsets="])
	# Remaining 'args' must be files
	files = args
	
	for option, value in opts:
		log.prNot(log.DEBUG, 'Parsing common: %s:%s' % (option, value))
		if option in ["-v", "--verbose"]: log.VERBOSITY += 1
		if option in ["-h", "--help"]: print_help(tool)
		if option in ["-s", "--stats"]: params['stats'] = True
		if option in ["-i", "--informat"]: params['informat'] = value
		if option in ["--ff"]: params['flatfield'] = os.path.realpath(value)
		if option in ["--fm"]: params['flatmulti'] = float(value)
		if option in ["--df"]: params['darkfield'] = os.path.realpath(value)
		if option in ["--dm"]: params['darkmulti'] = float(value)
		if option in ["--mf"]: params['maskfile'] = os.path.realpath(value)
		if option in ["--norm"]: params['norm'] = True
		if option in ["--nonorm"]: params['norm'] = False
		if option in ["--plot"]: params['plot'] = True
		if option in ["--noplot"]: params['plot'] = False
		if option in ["-f", "--file"]: params['file'] = os.path.realpath(value)
		if option in ["--crop"]: params['crop'] = value.split(',')
	
	# Parse tool-specific parameters now
	# ==================================
	
	if (tool == 'convert'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing convert: %s:%s' % (option, value))
			if option in ["-o", "--outformat"]: params['outformat'] = value
			if option in ["--scale"]: params['scale'] = float(value)
			if option in ["--intclip"]: params['intclip'] = value.split(',')
	elif (tool == 'samask'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing samask: %s:%s' % (option, value))
			if option in ["--rad"]: params['rad'] = float(value)
			if option in ["--shape"]: params['shape'] = value
			if option in ["--size"]: params['size'] = value.split(',')
			if option in ["--pitch"]: params['pitch'] = value.split(',')
			if option in ["--xoff"]: params['xoff'] = value.split(',')
			if option in ["--disp"]: params['disp'] = value.split(',')
			if option in ["--scale"]: params['scale'] = float(value)
	elif (tool == 'sfmask'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing sfmask: %s:%s' % (option, value))
			if option in ["--sfsize"]: params['sfsize'] = value.split(',')
			if option in ["--sasize"]: params['sasize'] = value.split(',')
			if option in ["--overlap"]: params['overlap'] = value.split(',')
			if option in ["--border"]: params['border'] = value.split(',')
	elif (tool == 'saopt'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing saopt: %s:%s' % (option, value))
			if option in ["--saifac"]: params['saifac'] = float(value)
			if option in ["--rad"]: params['rad'] = float(value)
	elif (tool == 'saupd'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing saupd: %s:%s' % (option, value))
			if option in ["--offsets"]: params['offsets'] = os.path.realpath(value)
	elif (tool == 'shifts'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing shift: %s:%s' % (option, value))
			if option in ["-r", "--range"]: params['shrange'] = N.int32(value)
			if option in ["--safile"]: params['safile'] = os.path.realpath(value)
			if option in ["--sffile"]: params['sffile'] = os.path.realpath(value)
			if option in ["-n", "--nref"]: params['nref'] = N.int32(value)
	
	return (tool, params, files)


def get_defaults(tool):
	"""
	Set default paramets for 'tool' in a dict, and return it .
	"""
	default = {}
	# Common defaults
	default['flatfield'] = False
	default['flatmulti'] = 1
	default['darkfield'] = False
	default['darkmulti'] = 1
	default['maskfile'] = False
	default['norm'] = True
	default['stats'] = False
	default['informat'] = _FORMAT_ANA
	default['file'] = False
	default['plot'] = False
	# Convert defaults:
	default['scale'] = 1.0
	default['intclip'] = False
	default['crop'] = False
	default['outformat'] = _FORMAT_FITS
	# Samask defaults:
	default['rad'] = 1024
	default['shape'] = 'circular'
	default['size'] = ['128','128']
	default['pitch'] = ['164','164']
	default['xoff'] = ['0', '0.5']
	default['disp'] = ['0','0']
	# sfmask defaults
	default['sfsize'] = ['16', '16']
	default['sasize'] = ['0', '0']
	default['overlap'] = ['0.5', '0.5']
	default['border'] = ['6', '6']
	# Saopt defaults
	default['saifac'] = 0.7
	# Saupd defaults
	default['offsets'] = False
	
	# Shift defaults
	default['shrange'] = 7
	default['safile'] = False
	default['sffile'] = False
	default['nref'] = 4
	
	# Internal configuration, cannot be changed command-line
	# pixel offset for any operation within a masked image
	default['_pix'] = 2
	
	return default


def check_params(tool, params):
	"""
	Check whether tool 'action' exists and whether parameters 'params' are sane 
	(files exists, options are allowed, etc.)
	"""
	
	# Tool must be valid
	# ===================================================================
	
	if tool not in _TOOLS:
		print >> sys.stderr, os.path.basename(sys.argv[0]) + ": Tool '%s' not available." % (tool)
		print >> sys.stderr, "\t for help use --help"
		sys.exit(2)
	
	# Common options requirements
	# ===================================================================
	
	# Check outformat
	if params['outformat'] not in _OUTFORMATS:
		log.prNot(log.ERROR, "Unsupported output format '%s'" %(params['outformat']))
	# Check informat
	if params['informat'] not in _INFORMATS:
		log.prNot(log.ERROR, "Unsupported input format '%s'" % (params['informat']))
	# Intclip should be (float, float)
	if (params['intclip'] is not False):
		try: 
			tmp = params['intclip'][1]
			params['intclip'] = N.array(params['intclip'][0:2]).astype(N.float)
		except: log.prNot(log.ERROR, "intclip invalid, should be <float>,<float>.")
	# Crop should be (int, int, int, int)
	if (params['crop'] is not False):
		try: 
			tmp = params['crop'][3]
			params['crop'] = N.array(params['crop'][0:4]).astype(N.float)
		except: log.prNot(log.ERROR, "crop invalid, should be 4 floats.")
	# Flatfield must exist
	if (params['flatfield']) and \
		(not os.path.exists(params['flatfield'])):
		log.prNot(log.ERROR, "flatfield '%s' does not exist." % \
		 	(params['flatfield']))
	# Darkfield must exist
	if (params['darkfield']) and \
		(not os.path.exists(params['darkfield'])):
		log.prNot(log.ERROR, "darkfield '%s' does not exist." % \
		 	(params['darkfield']))
	# Makfile needs to exist
	if (params['maskfile']) and \
		(not os.path.exists(params['maskfile'])):
		log.prNot(log.ERROR, "maskfile '%s' does not exist." % \
		 	(params['maskfile']))
	# File should not exist, find a new file if it does
	if (params['file']):
		libfile.saveOldFile(params['file'], postfix='.old', maxold=5)
	# Shape should be 'circular' or 'square'
	if (params['shape'] not in ['square', 'circular']):
		log.prNot(log.ERROR, "shape invalid, should be 'circular' or 'square'")
	# Size should be (int, int)
	try: 
		tmp = params['size'][1]
		params['size'] = N.array(params['size'][0:2]).astype(N.int32)
	except: log.prNot(log.ERROR, "size invalid, should be <int>,<int>.")
	# pitch should be int, int
	try: 
		tmp = params['pitch'][1]
		params['pitch'] = N.array(params['pitch'][0:2]).astype(N.int32)
	except: log.prNot(log.ERROR, "pitch invalid, should be <int>,<int>.")
	# xoff should be floats
	try: 
		tmp = params['xoff'][1]
		params['xoff'] = N.array(params['xoff'][0:2]).astype(N.float)
	except: log.prNot(log.ERROR, "xoff invalid, should be <float>,<float>.")
	# disp should be int, int
	try: params['disp'] = N.array(params['disp'][0:2]).astype(N.int32)
	except: log.prNot(log.ERROR, "disp invalid, should be <int>,<int>.")
	
	try: 
		tmp = params['sfsize'][1]
		params['sfsize'] = N.array(params['sfsize'][0:2]).astype(N.int32)
	except: log.prNot(log.ERROR, "sfsize invalid, should be <int>,<int>.")
	try: 
		tmp = params['sasize'][1]
		params['sasize'] = N.array(params['sasize'][0:2]).astype(N.int32)
	except: log.prNot(log.ERROR, "sasize invalid, should be <int>,<int>.")
	try: 
		tmp = params['overlap'][1]
		params['overlap'] = N.array(params['overlap'][0:2]).astype(N.float32)
	except:
		log.prNot(log.ERROR, "overlap invalid, should be <float>,<float>.")
	try: 
		tmp = params['border'][1]
		params['border'] = N.array(params['border'][0:2]).astype(N.int32)
	except:
		log.prNot(log.ERROR, "border invalid, should be <int>,<int>.")
	# Offset file needs to exist
	if (params['offsets']) and \
		(not os.path.exists(params['offsets'])):
		log.prNot(log.ERROR, "offset file '%s' does not exist." % \
		 	(params['offsets']))

	# Requirements depending on tools (where defaults are not sufficient)
	# ===================================================================
	if (tool == 'saopt'):
		# need flatfield, maskfile
		if (params['flatfield']) and \
			(not os.path.exists(params['flatfield'])):
			log.prNot(log.ERROR, "Tool 'saopt' requires flatfield.")
		if (params['maskfile']) and \
			(not os.path.exists(params['maskfile'])):
			log.prNot(log.ERROR, "Tool 'saopt' requires maskfile.")
	elif (tool == 'shifts'):
		# need safile and sffile
		if (params['safile']) and \
			(not os.path.exists(params['safile'])):
			log.prNot(log.ERROR, "safile '%s' does not exist." % \
			 	(params['safile']))
		if (params['sffile']) and \
			(not os.path.exists(params['sffile'])):
			log.prNot(log.ERROR, "sffile '%s' does not exist." % \
			 	(params['sffile']))
	elif (tool == 'sfmask'):
		# sasize > sfsize, both must be int, int
		if (params['sasize'] < params['sfsize']).any():
			log.prNot(log.ERROR, "sasize must be bigger than sfsize.")
		if (params['overlap'] > 1).any() or (params['overlap'] < 0).any():
			log.prNot(log.ERROR, "overlap must be between 0 and 1.")
	elif (tool == 'saupd'):
		# Saupd needs maskfile, offset file
		if (params['maskfile']) and \
			(not os.path.exists(params['maskfile'])):
			log.prNot(log.ERROR, "Tool 'saupd' requires maskfile.")
		if (params['offsets']) and \
			(not os.path.exists(params['offsets'])):
			log.prNot(log.ERROR, "Tool 'saupd' requires offsets file.")
	
	# Done


### ==========================================================================
### Tool classes
### ==========================================================================

class Tool(object):
	"""Generic Tool class with common functions."""
	
	def __init__(self, files, params):
		self.files = files
		self.params = params
		# Save some options common for all tools
		self.informat = params['informat']
		self.flatfield = params['flatfield']
		self.flatmulti = params['flatmulti']
		self.flatdata = None
		self.darkfield = params['darkfield']
		self.darkmulti = params['darkmulti']
		self.darkdata = None
		self.gaindata = None
		self.maskfile = params['maskfile']
		self.mask = None
		self.norm = params['norm']
		self.crop = params['crop']
		self.plot = params['plot']
		self.file = params['file']
		if (self.maskfile):
			log.prNot(log.DEBUG, "Loading mask files.")
			(self.nsa, self.saccdpos, self.saccdsize) = \
			 	libsh.loadSaSfConf(self.maskfile)
		
		# This is the dataid, containing the direct parent directory, the base of 
		# the file, and the range of extensions we have:
		
		log.prNot(log.INFO, "Processing %d files" % \
			(len(self.files)))
	
	
	def load(self, filename):
		log.prNot(log.DEBUG, "Loading '%s'" % (filename))
		if not os.path.exists(filename):
			log.prNot(log.WARN, "File '%s' does not exist" % (filename))
			return None
		# Load data, using different methods depending on type
		if (self.informat == _FORMAT_ANA):
			data = self.__anaload(filename)
		elif (self.informat == _FORMAT_FITS):
			data = self.__fitsload(filename)	
		elif (self.informat == _FORMAT_NPY):
			data = self.__npyload(filename)	
		else:
			log.prNot(log.WARN, "Filetype unsupported." % (filename))			
			return None
		self.origres = data.shape
		if (self.crop is not False):
			data = data[self.crop[1]:self.crop[1] + self.crop[3], \
				self.crop[0]:self.crop[0] + self.crop[2]]

		return data
	
	
	def __anaload(self, filename):
		import pyana
		return pyana.getdata(filename)
	
	
	def __fitsload(self, filename):
		import pyfits
		return pyana.getdata(filename)
	
	
	def __npyload(self, filename):
		import numpy
		return numpy.load(filename)
	
	
	def darkflat(self, data):
		# Init dark-flat fields
		if (self.flatfield or self.darkfield):
			self.__initdarkflat()
			log.prNot(log.DEBUG, "Dark-flatfielding data. Dark avg: %.4g, gain avg: %.4g" % (N.mean(self.darkdata), N.mean(self.gaindata)))		
			# Now process the frame
			return (data-self.darkdata) * self.gaindata
		return data
	
	
	def __initdarkflat(self):
		# Get flats and darks, if not already present
		if (self.flatdata is None):
			if (self.flatfield):
				log.prNot(log.DEBUG, "Loading flatfield...")
				self.flatdata = self.load(self.flatfield)
				self.flatdata = self.flatdata.astype(N.float32)
				self.flatdata /= 1.0*self.flatmulti
				log.prNot(log.DEBUG, "Flatfield average: %.6g" % N.mean(self.flatdata))
			else:
				# Maybe we don't want flatfielding, in that case set it to 1.0
				log.prNot(log.DEBUG, "Not flatfielding, setting to 1.0")
				self.gaindata = N.float32(1.0)
		if (self.darkdata is None):
			if (self.darkfield):
				log.prNot(log.DEBUG, "Loading darkfield...")
				self.darkdata = self.load(self.darkfield)
				self.darkdata = self.darkdata.astype(N.float32)
				self.darkdata /= 1.0*self.darkmulti
				log.prNot(log.DEBUG, "Darkfield average: %.6g" % N.mean(self.darkdata))
			else:
				# Maybe we don't want darkfielding, in that case set it to 1.0
				log.prNot(log.DEBUG, "Not darkfielding, setting to 0.0")
				self.darkdata = N.float32(0.0)
		if (self.gaindata is None):
			# Make a gain for faster processing
			self.gaindata = 1.0/(self.flatdata - self.darkdata)
			#self.gaindata /= N.mean(self.gaindata)
	
	
	def fitssave(self, data, filepath, overwrite=True):
		"""
		Save 'data' as FITS file to 'filepath'.
		"""
		import pyfits
		log.prNot(log.DEBUG, "Tool.fitssave(): Saving data to '%s'." % (filepath))
		pyfits.writeto(filepath, data, clobber=overwrite)
	
	
	def anasave(self, data, filepath, compressed=1):
		"""
		Save 'data' as ANA file to 'filepath'. Can be compressed (default: yes).
		"""
		import pyana
		log.prNot(log.DEBUG, "Tool.anasave(): Saving data to '%s'." % (filepath))
		pyana.fzwrite(filepath, data, compressed)
	
	
	def npysave(self, data, filepath):
		"""
		Save 'data' as npy file to 'filepath'.
		"""
		log.prNot(log.DEBUG, "Tool.npysave(): Saving data to '%s'." % (filepath))
		N.save(data, filepath)
	
	
	def pngsave(self, data, filepath, scale=True):
		"""
		Save 'data' as PNG file to 'filepath'. Data can be scaled to full range
		(default: yes).
		"""
		log.prNot(log.DEBUG, "Tool.pngsave(): Saving data to '%s'." % (filepath))
		if (scale):
			# Scale the values to 0-255
			maxval = N.max(data)
			minval = N.min(data)
			scdata = (255.0*(data - minval)/(maxval - minval))
			scdata = (scdata.astype(N.uint8))
		else: scdata = data
		if (scdata.shape[1] % 4 != 0): 
			raise RuntimeError("Cannot save PNG files with horizontal resolution not a multiple of 4 (cairo bug).")
		
		import cairo
		# Create surface from data
		surf = cairo.ImageSurface.create_for_data(scdata, \
		 	cairo.FORMAT_A8, scdata.shape[1], scdata.shape[0])
		# Create context for mirror matrix
		ctx = cairo.Context(surf)
		mat = cairo.Matrix(1, 0, 0, -1, 0, surf.get_height())
		ctx.transform(mat)
		# Save to disk
		cairo.ImageSurface.write_to_png(surf, filepath)
	
	
	def maskimg(self, data):
		"""
		Apply a mask on an image, set all values outside the mask to the minimum 
		value inside the mask, so it will appear black.
		"""
		if (self.maskfile is False):
			return data
		log.prNot(log.DEBUG, "Masking image if necessary, res: %d,%d" % \
		 	(tuple(self.origres)))
		self.__initmask(self.origres)
		data[self.mask == False] = N.min(data[self.mask])
		if (self.norm):
			# Normalize data per subapt
			for p in self.saccdpos:
				size = self.saccdsize
				if (self.crop is not False):
					p = p - self.crop[0:2]
					# Skip subapertures that lie outside the cropped field of view
					if (p + self.saccdsize < 0).any(): continue
					if (p > self.crop[2:]).any(): continue
					# Make sure positions are positive, and adapt size to that
					poff = N.clip(p, 0, p.max()) - p
					size = size - poff
					p = N.clip(p, 0, p.max())
				avg = N.mean(data[\
					p[1]:p[1]+size[1], \
					p[0]:p[0]+size[0]])
				data[\
					p[1]:p[1]+size[1], \
					p[0]:p[0]+size[0]] /= avg
				log.prNot(log.DEBUG, "maskimg(): normalizing, avg: %.3g" % (avg))
		
		return data
	
	
	def __initmask(self, res):
		if (self.mask is None) or (res != self.maskres):
			# We need to make a new mask here
			log.prNot(log.DEBUG, "maskimg(): (re-)initializing mask with resolution %d,%d" % (tuple(res)))
			(self.mask, self.maskborder) = \
				libsh.makeSubaptMask(self.saccdpos, self.saccdsize, res)
			self.maskres = self.mask.shape
			if (self.crop is not False):
				self.mask = self.mask[self.crop[1]:self.crop[1] + self.crop[3], \
					self.crop[0]:self.crop[0] + self.crop[2]]
				self.maskborder = self.maskborder[self.crop[1]:self.crop[1] + \
				 	self.crop[3], self.crop[0]:self.crop[0] + self.crop[2]]
				
	



class SubaptConfTool(Tool):
	"""Calculate configurations."""
	def __init__(self, files, params):
		super(SubaptConfTool, self).__init__(files, params)
		# Subaperture pattern radius
		self.rad = params['rad']
		# Subaperture size
		self.size = params['size']
		# Subaperture pitch
		self.pitch = params['pitch']
		# Aperture shape, either circular or square
		self.shape = params['shape']
		# Odd-row offset
		self.xoff = params['xoff']
		# Scale factor for complete pattern
		self.scale = params['scale']
		# Displacement vector for complete pattern
		self.disp = params['disp']
		# Output file
		if params['file']: self.file = params['file']
		else: self.file = './astooki-subaptconf.csv'
				
		self.run()
	
	
	def run(self):
		# Generate pattern
		(nsa, saccdpos, saccdsize) = \
			libsh.calcSubaptConf(self.rad, self.size, self.pitch, self.shape, \
			self.xoff, self.disp, self.scale)
		# Save to file
		libsh.saveSaSfConf(self.file, nsa, [-1,-1], saccdsize, saccdpos)
		if (self.plot):
			plfile = os.path.splitext(self.file)[0]+'-plot.eps'
			libplot.showSaSfLayout(plfile, saccdpos, saccdsize, \
				plrange=[[0, 2*self.rad]]*2)
		# Done
		
	


class SubfieldConfTool(Tool):
	"""Calculate subfield configurations."""
	def __init__(self, files, params):
		super(SubfieldConfTool, self).__init__(files, params)
		# Subaperture size
		self.sasize = params['sasize']
		# Subfield size
		self.sfsize = params['sfsize']
		# Overlap
		self.overlap = params['overlap']
		# Border 
		self.border = params['border']
		
		if params['file']: self.file = params['file']
		else: self.file = './astooki-subfieldconf.csv'
		
		self.run()
	
	
	def run(self):
		# Generate subfield positions
		effsize = self.sasize - 2*self.border
		pitch = self.sfsize * (1-self.overlap)
		nsf = effsize / pitch
		nsf = N.round(nsf-1)
		effpitch = effsize/(nsf+1)
		
		sfpos = self.border + \
			N.indices(nsf, dtype=N.float32).reshape(2,-1).T * effpitch
		sfpos = N.round(sfpos).astype(N.int32)
		totnsf = N.product(nsf).astype(N.int32)
		
		log.prNot(log.INFO, "Found %d x %d subfields." % tuple(nsf))
		log.prNot(log.INFO, "Size %d,%d" % tuple(self.sfsize))
		libsh.saveSaSfConf(self.file, totnsf, [-1,-1], self.sfsize, sfpos)
		if (self.plot):
			plfile = os.path.splitext(self.file)[0]+'-plot.eps'
			libplot.showSaSfLayout(plfile, sfpos, self.sfsize, \
				plrange=[[0, self.sasize[0]], [0, self.sasize[1]]])
		# Done
		
	


class SubaptUpdateTool(Tool):
	"""Update subaperture configurations with an offset."""
	def __init__(self, files, params):
		super(SubaptUpdateTool, self).__init__(files, params)
		# Output file
		tmp = os.path.splitext(self.maskfile)
		self.file = tmp[0]+'-updated'+tmp[1]
		self.offsets = params['offsets']
		
		self.run()
		
	
	
	def run(self):
		# Load mask file
		(nsa, pos, size) = \
		 	libsh.loadSaSfConf(self.maskfile)
		# load offsets
		off = libfile.loadData(self.offsets, ascsv=True)
		
		# Compare
		if (off.shape[0] != nsa):
			log.prNot(log.ERROR, "SubaptUpdateTool(): offsets not the same size as subaperture positions.")
		
		# Offset the new positions
		newpos = (pos + off).astype(N.int32)
		# Crop the subaperture size by twice the maximum offset
		newsize = (size - (N.max(off, axis=0)*2)).astype(N.int32)
		# Store
		libsh.saveSaSfConf(self.file, nsa, [-1,-1], newsize, newpos)
		if (self.plot):
			plfile = os.path.splitext(self.file)[0]+'-plot.eps'
			plran = [[0, N.ceil(1600./512)*512]]*2
			libplot.showSaSfLayout(plfile, newpos, newsize, \
				plrange=plran)
		# Done
	


class SubaptOptTool(Tool):
	"""Optimize subaperture configurations."""
	def __init__(self, files, params):
		super(SubaptOptTool, self).__init__(files, params)
		# Output file
		tmp = os.path.splitext(self.maskfile)
		if params['file']: self.file = params['file']
		else: self.file = tmp[0]+'-optimized'+tmp[1]
		self.saifac = params['saifac']
		self.rad = params['rad']
		
		self.run()
	
	
	def run(self):
		# Load flatfield
		flatimg = self.load(self.flatfield)
		# Load old configuration
		(nsa, pos, size) = \
		 	libsh.loadSaSfConf(self.maskfile)
		# Optimize pattern
		(onsa, opos, osize) = libsh.optSubapConf(flatimg, pos, size, self.saifac)
		# Save to file
		libsh.saveSaSfConf(self.file, onsa, [-1,-1], osize, opos)
		if (self.plot):
			plfile = os.path.splitext(self.file)[0]+'-plot.eps'
			libplot.showSaSfLayout(plfile, opos, osize, \
				plrange=[[0, 2*self.rad]]*2)
		# Done
	


class ShiftTool(Tool):
	"""Calculate 'static' shifts between subapetures."""
	def __init__(self, files, params):
		super(ShiftTool, self).__init__(files, params)
		self.shrange = params['shrange']
		self.safile = params['safile']
		self.sffile = params['sffile']
		self.nref = params['nref']
		# Load safile and sffile
		(self.nsa, self.saccdpos, self.saccdsize) = \
			libsh.loadSaSfConf(self.safile)
		(self.nsf, self.sfccdpos, self.sfccdsize) = \
			libsh.loadSaSfConf(self.sffile)
		if (len(self.files) > 0):
			pardir = \
			 	os.path.basename(os.path.dirname(os.path.realpath(self.files[0])))
			basefile = os.path.splitext(os.path.basename(self.files[0]))[0]
			begid = int((os.path.splitext(os.path.basename(self.files[0]))[1])[1:])
			endid = int((os.path.splitext(os.path.basename(self.files[-1]))[1])[1:])
			self.dataid = "%s_%s_%d-%d" % (pardir, basefile, begid, endid)
		else:
			self.dataid = "astooki-generic"			
		# Make sure we have a file to save results to
		if (self.file is False):
			self.file = os.path.realpath(self.dataid)
		# Run analysis
		self.run()
	
	
	def run(self):
		import libshifts as ls
		# Process files
		allshifts = []
		allfiles = []
		for f in self.files:
			base = os.path.basename(f)
			log.prNot(log.INFO, "Measuring shifts for %s." % (base))
			# Load file
			img = self.load(f)
			if (img is None): 
				log.prNot(log.DEBUG, "Skipping %s, could not read file." % (base))
				continue
			allfiles.append(base)
			# Dark-flat file if needed
			dfimg = self.darkflat(img)
			
			# Measure shift
			imgshifts = ls.calcShifts(dfimg, self.saccdpos, self.saccdsize, \
			 	self.sfccdpos, self.sfccdsize, method=ls.COMPARE_ABSDIFFSQ, \
			 	extremum=ls.EXTREMUM_2D9PTSQ, refmode=ls.REF_BESTRMS, \
			 	refopt=self.nref, shrange=[self.shrange, self.shrange], \
			 	subfields=None, corrmaps=None)
			
			allshifts.append(imgshifts)
		
		# Process results, store to disk
		log.prNot(log.INFO, "Done, saving results to disk @ '%s'." % (self.file))
		allshifts = N.array(allshifts)
		# Store the list of files where we save data to
		files = {}
		files['shifts'] = libfile.saveData(self.file + '-shifts', \
			allshifts, asnpy=True, asfits=True)
		files['saccdpos'] = libfile.saveData(self.file + '-saccdpos', \
		 	self.saccdpos, asnpy=True)
		files['sfccdpos'] = libfile.saveData(self.file + '-sfccdpos', \
		 	self.sfccdpos, asnpy=True)
		files['saccdsize'] = libfile.saveData(self.file + '-saccdsize', \
		 	self.saccdsize, asnpy=True)
		files['sfccdsize'] = libfile.saveData(self.file + '-sfccdsize', \
		 	self.sfccdsize, asnpy=True)
		cpos = self.saccdpos.reshape(-1,1,2) + self.sfccdpos.reshape(1,-1,2) + \
			self.sfccdsize.reshape(1,1,2)/2.0
		files['sasfpos-c'] = libfile.saveData(self.file + '-sasfpos-c', \
		 	cpos, asnpy=True, asfits=True)
		files['files'] = libfile.saveData(self.file + '-files', \
		 	allfiles, asnpy=True, ascsv=True, csvfmt='%s')
		
		# If we have only one subfield, also calculate 'static' shifts
		if (len(self.sfccdpos) == 1):
			log.prNot(log.INFO, "Calculating static offsets.")
			(soff, sofferr) =libsh.procStatShift(allshifts[:,:,:,0,:])
			files['files'] = libfile.saveData(self.file + '-offset', \
			 	soff, asnpy=True, ascsv=True)
			files['files'] = libfile.saveData(self.file + '-offset-err', \
			 	sofferr, asnpy=True)
			libplot.plotShifts(self.file + '-offset-plot', allshifts, \
				self.saccdpos, self.saccdsize, self.sfccdpos, self.sfccdsize, \
				plorigin=(0,0), plrange=(2048, 2048), mag=7.0, allsh=True, \
				 title='Static offsets for' + self.dataid,  legend=True)
		
		libfile.saveData(self.file + '-meta', files, aspickle=True)

	


class StatsTool(Tool):
	"""Calculate statistics on files"""
	def __init__(self, files, params):
		super(StatsTool, self).__init__(files, params)
		if (len(self.files) > 0):
			pardir = \
			 	os.path.basename(os.path.dirname(os.path.realpath(self.files[0])))
			basefile = os.path.splitext(os.path.basename(self.files[0]))[0]
			begid = int((os.path.splitext(os.path.basename(self.files[0]))[1])[1:])
			endid = int((os.path.splitext(os.path.basename(self.files[-1]))[1])[1:])
			self.dataid = "%s_%s_%d-%d" % (pardir, basefile, begid, endid)
		else:
			self.dataid = "astooki-generic"
		self.file = os.path.realpath(self.dataid)
		self.run()
	
	
	def run(self):
		# Store all stats and all files
		allstats = []
		allfiles = []
		for f in self.files:
			base = os.path.basename(f)
			allfiles.append(base)
			allstats.append(self.dowork(f))
		
		# Process stats for all files, display average
		allstats = N.array(allstats)
		allfiles = (N.array(allfiles)).reshape(-1,1)
		all_avg = N.mean(allstats, axis=0)
		all_std = (N.var(allstats, axis=0))**0.5
		log.prNot(log.INFO, "all %d: mean: %.4g+-%.4g rms: %.4g+-%.4g (%.3g%%+-%.4g)" % \
			(allstats.shape[0], all_avg[0], all_std[0], all_avg[2], all_std[2], all_avg[3], all_std[3]))
		# Save results if requested
		nf = libfile.saveData(self.file + '-stats', allstats, asnpy=True)
		hdr = ['filename, path=%s' % \
		 	(os.path.dirname(os.path.realpath(self.files[0]))), 'avg', \
		 	'std', 'rms', 'fractional rms']
		cf = libfile.saveData(self.file + '-stats', N.concatenate((allfiles, \
		 	allstats), axis=1), ascsv=True, csvhdr=hdr, csvfmt='%s')
	
	
	def dowork(file):
		# Load file
		img = self.load(f)
		if (img is None): 
			log.prNot(log.DEBUG, "Skipping %s, could not read file." % (base))
			return False
		# Dark-flat file if needed
		dfimg = self.darkflat(img)
		data = dfimg
		# If a maskfile is given, only calculate stats within the subaperture.
		if (self.maskfile):
			size = self.saccdsize
			substat = []
			for p in self.saccdpos:
				_sub = data[\
					p[1]:p[1]+size[1], \
					p[0]:p[0]+size[0]]
				avg = N.mean(_sub)
				std = N.var(_sub)**0.5
				rms = (N.sum((_sub-avg)**2.0)/_sub.size)**0.5
				rmsrat = 100.0*rms/avg
				substat.append([avg, std, rms, rmsrat])
			substat = N.mean(N.array(substat), axis=0)
			log.prNot(log.INFO, "%s: mean: %.4g std: %.4g rms: %.4g (%.3g%%)" % \
				(base, substat[0], substat[1], substat[2], substat[3]))
			return list(substat)
		else:
			# If no maskfile given, calculate stats for all pixels
			r = (N.min(data), N.max(data))
			avg = N.mean(data)
			std = N.var(data)**0.5
			rms = (N.sum((data-avg)**2.0)/data.size)**0.5
			rmsrat = 100.0*rms/avg	
			log.prNot(log.INFO, "%s: mean: %.4g std: %.4g range: %.3g--%.3g rms: %.4g (%.3g%%)" % \
				(base, avg, std, r[0], r[1], rms, rmsrat))
			return [avg, std, rms, rmsrat]
		
		# done
	



class ConvertTool(Tool):
	"""Convert files to different format"""
	def __init__(self, files, params):
		super(ConvertTool, self).__init__(files, params)
		self.file = params['file']
		self.outformat = params['outformat']
		self.scale = params['scale']
		self.intclip = params['intclip']
		self.run()
	
	
	def run(self):
		# Process files
		for f in self.files:
			base = os.path.basename(f)
			log.prNot(log.INFO, "Converting file '%s' to %s" % (base, \
			 	self.outformat))
			# Load file
			img = self.load(f)
			if (img is None): continue
			# Dark-flat file if needed
			dfimg = self.darkflat(img)
			# Mask if needed
			data = self.maskimg(dfimg)
			if (self.scale != 1.0):
				import scipy.ndimage
				sc = self.scale
				orig = data.shape[0]
				nsc = (N.round(orig * sc / 8)*8)/orig
				log.prNot(log.INFO, "Scaling image by approximately %g." % self.scale)
				data = S.ndimage.zoom(data, nsc, mode='wrap')
			# Crop intensity if necessary
			if (self.intclip is not False):
				log.prNot(log.INFO, "Clipping intensity to %g--%g." % \
				 	tuple(self.intclip))
				data = N.clip(data, self.intclip[0], self.intclip[1])
			# Save again
			if (len(self.files) == 1 and self.file): savefile = self.file
			else: savefile = f+'.'+self.outformat
			log.prNot(log.INFO, "Saving '%s' as %s in '%s'." % \
				(base, self.outformat, os.path.basename(savefile)))
				
			if (self.outformat == _FORMAT_PNG): self.pngsave(data, savefile)
			elif (self.outformat == _FORMAT_FITS): self.fitssave(data, savefile)
			elif (self.outformat == _FORMAT_ANA): self.anasave(data, savefile)
			elif (self.outformat == _FORMAT_NPY): self.npysave(data, savefile)
		
	


### ==========================================================================
### Helper functions
### ==========================================================================

def print_help(tool='common', out=sys.stdout):
	if tool == 'common':
		print >> out, help_message
		sys.exit(0)
	else:
		print >> out, 'No specific help for this tool.'
		sys.exit(0)


if __name__ == "__main__":
	sys.exit(main())
