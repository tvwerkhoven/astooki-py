#!/usr/bin/env /sw/bin/python2.5
# encoding: utf-8
"""
@filename astooki.py
@author Tim van Werkhoven (tim@astro.su.se)
@date 20090422 14:34
@brief AsTooki: the astronomical toolkit - imagemagick for astronomers.

Created by Tim on 2009-04-22.
Copyright (c) 2009 Tim van Werkhoven. All rights reserved.

$Id$
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
 saopt                       Optimze a subaperture mask with a flat
 shifts                      Measure image shifts in various subfields and 
                               subimages

Input formats supported
 ana                         ANA format
 fits                        Flexible Image Transport System

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

Saopt options
 -f, --file=FILEPATH         file to store optimized configuration to
     --saifac=FLOAT          intensity drop-off factor considered 'dark' [0.7]
     --rad=RAD               radius of the subaperture pattern, used for 
                               plotting [1024]

Shifts options
 -r, --range=INT             shift range to use for cross-correlation [7]
     --safile=FILEPATH       subaperture locations, same format as maskfile
     --sffile=FILEPATH       subfield locations w.r.t. subaperture locations
 -n, --nref=INT              number of references to use [4]
''' % (VERSION, DATE, AUTHOR)

### Supported formats
_FORMAT_ANA = 'ana'
_FORMAT_FITS = 'fits'
_FORMAT_PNG = 'png'
_FORMAT_NPY = 'npy'
_INFORMATS = (_FORMAT_ANA, _FORMAT_FITS)
_OUTFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_PNG, _FORMAT_NPY)
# Tools available
_TOOLS = ('convert', 'stats', 'samask', 'saopt', 'shifts')

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
	elif (tool == 'saopt'): SubaptOptTool(files, params)
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
		if tool not in _TOOLS:
			raise Exception
	except:
		print >> sys.stderr, os.path.basename(sys.argv[0]) + ": Syntax incorrect."
		print >> sys.stderr, "\t for help use --help"
		sys.exit(2)
	if tool in ["-h", "--help"]: print_help()
		
	
	# Parse common options first
	# ==========================
	
	params = get_defaults(tool)
	
	opts, args = getopt.getopt(argv[2:], "vhsi:o:f:r:n:", ["verbose", "help", "stats", "informat=", "ff=", "fm=", "df=", "dm=", "mf=", "outformat=", "intclip=", "crop=", "file=", "scale=", "rad=", "shape=", "size=", "pitch=", "xoff=", "disp=", "plot", "noplot", "norm", "nonorm", "saifac=", "range=", "sffile=", "safile=", "nref="])
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
	elif (tool == 'saopt'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing saopt: %s:%s' % (option, value))
			if option in ["--saifac"]: params['saifac'] = float(value)
			if option in ["--rad"]: params['rad'] = float(value)
	elif (tool == 'shifts'):
		for option, value in opts:
			log.prNot(log.DEBUG, 'Parsing shift: %s:%s' % (option, value))
			if option in ["-r", "--range"]: params['shrange'] = float(value)
			if option in ["--safile"]: params['safile'] = value
			if option in ["--sffile"]: params['sffile'] = value
			if option in ["-n", "--nref"]: params['nref'] = int(value)
	
	return (tool, params, files)


def get_defaults(tool):
	"""
	Set default paramets for 'tool' in a dict, and return it .
	"""
	default = {}
	# Common defaults
	default['flatfield'] = 'none'
	default['flatmulti'] = 1
	default['darkfield'] = 'none'
	default['darkmulti'] = 1
	default['maskfile'] = 'none'
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
	# Saopt defaults
	default['saifac'] = 0.7
	# Shift defaults
	default['shrange'] = 7
	default['safile'] = 'none'
	default['sffile'] = 'none'
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
	if (params['flatfield'] is not 'none') and \
		(not os.path.exists(params['flatfield'])):
		log.prNot(log.ERROR, "Flatfield does not exist.")
	# Darkfield must exist
	if (params['darkfield'] is not 'none') and \
		(not os.path.exists(params['darkfield'])):
		log.prNot(log.ERROR, "darkfield does not exist.")
	# Makfile needs to exist
	if (params['maskfile'] is not 'none') and \
		(not os.path.exists(params['maskfile'])):
		log.prNot(log.ERROR, "maskfile '%s' does not exist." % (params['maskfile']))
	# File should not exist, find a new file if it does
	if (params['file']):
		libfile.saveOldFile(params['file'], postfix='.old', maxold=5)
	# Shape should be 'circular' or 'square'
	if (params['shape'] not in ['square', 'circular']):
		log.prNot(log.ERROR, "shape invalid, should be 'circular' or 'square'")
	# Size should be (int, int)
	try: 
		tmp = params['size'][1]
		params['size'] = N.array(params['size'][0:2]).astype(N.int)
	except: log.prNot(log.ERROR, "size invalid, should be <int>,<int>.")
	# pitch should be int, int
	try: 
		tmp = params['pitch'][1]
		params['pitch'] = N.array(params['pitch'][0:2]).astype(N.int)
	except: log.prNot(log.ERROR, "pitch invalid, should be <int>,<int>.")
	# xoff should be floats
	try: 
		tmp = params['xoff'][1]
		params['xoff'] = N.array(params['xoff'][0:2]).astype(N.float)
	except: log.prNot(log.ERROR, "xoff invalid, should be <float>,<float>.")
	# disp should be int, int
	try: params['disp'] = N.array(params['disp'][0:2]).astype(N.int)
	except: log.prNot(log.ERROR, "disp invalid, should be <int>,<int>.")
	
	# Requirements depending on tools (where defaults are not sufficient)
	# ===================================================================
	if (tool == 'saopt'):
		# need flatfield, maskfile
		if (not os.path.exists(params['flatfield'])):
			log.prNot(log.ERROR, "Tool 'saopt' requires flatfield.")
		if (not os.path.exists(params['maskfile'])):
			log.prNot(log.ERROR, "Tool 'saopt' requires maskfile.")
	elif (tool == 'shifts'):
		# need safile and sffile
		params['safile'] = os.path.realpath(params['safile'])
		if (not os.path.exists(params['safile'])):
			log.prNot(log.ERROR, "Tool 'shifts' requires safile (%s)." % \
			 	(params['safile']))
		params['sffile'] = os.path.realpath(params['sffile'])
		if (not os.path.exists(params['sffile'])):
			log.prNot(log.ERROR, "Tool 'shifts' requires sffile (%s)." % \
			 	(params['sffile']))
	
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
		if (self.maskfile is not 'none'):
			log.prNot(log.DEBUG, "Loading mask files.")
			(self.nsa, self.saccdpos, self.saccdsize) = \
			 	libsh.loadSaSfConf(self.maskfile)
		
		log.prNot(log.INFO, "Processing %d files." % (len(self.files)))
	
	
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
		else:
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
	
	
	def darkflat(self, data):
		# Init dark-flat fields
		self.__initdarkflat()
		# Now process the frame
		log.prNot(log.DEBUG, "Dark-flatfielding data. Dark avg: %.4g, gain avg: %.4g" % (N.mean(self.darkdata), N.mean(self.gaindata)))
		return (data.astype(N.float)-self.darkdata) * self.gaindata
	
	
	def __initdarkflat(self):
		# Get flats and darks, if not already present
		if (self.flatdata is None):
			if (self.flatfield is not 'none'):
				log.prNot(log.DEBUG, "Loading flatfield...")
				self.flatdata = self.load(self.flatfield)
				self.flatdata = self.flatdata.astype(N.float)
				self.flatdata /= 1.0*self.flatmulti
				log.prNot(log.DEBUG, "Flatfield average: %.6g" % N.mean(self.flatdata))
			else:
				# Maybe we don't want flatfielding, in that case set it to 1.0
				log.prNot(log.DEBUG, "Not flatfielding, setting to 1.0")
				self.gaindata = 1.0
		if (self.darkdata is None):
			if (self.darkfield is not 'none'):
				log.prNot(log.DEBUG, "Loading darkfield...")
				self.darkdata = self.load(self.darkfield)
				self.darkdata = self.darkdata.astype(N.float)
				self.darkdata /= 1.0*self.darkmulti
				log.prNot(log.DEBUG, "Darkfield average: %.6g" % N.mean(self.darkdata))
			else:
				# Maybe we don't want darkfielding, in that case set it to 1.0
				log.prNot(log.DEBUG, "Not darkfielding, setting to 0.0")
				self.darkdata = 0.0
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
		log.prNot(log.DEBUG, "Masking image.")
		if (self.maskfile is 'none'):
			return data
		self.__initmask(self.origres)
		data[self.mask == False] = N.min(data[self.mask])
		if (self.norm):
			# Normalize data per subapt
			for p in self.saccdpos:
				if (self.crop is not False):
					p -= self.crop[0:2]
					if (p < 0).any(): continue
					if (p > self.crop[2:]).any(): continue
				avg = N.mean(data[\
					p[1]:p[1]+self.saccdsize[1], \
					p[0]:p[0]+self.saccdsize[0]])
				data[\
					p[1]:p[1]+self.saccdsize[1], \
					p[0]:p[0]+self.saccdsize[0]] /= avg
				log.prNot(log.DEBUG, "maskimg(): normalizing, avg: %.3g" % (avg))
		
		return data
	
	
	def __initmask(self, res):
		if (self.mask is None) or (res != self.mask.shape):
			# We need to make a new mask here
				(self.mask, self.maskborder) = \
					libsh.makeSubaptMask(self.saccdpos, self.saccdsize, res)
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
		# Run analysis
		self.run()
	
	
	def run(self):
		import libshifts as ls
		# Process files
		allshifts = []
		for f in self.files:
			base = os.path.basename(f)
			# Load file
			img = self.load(f)
			if (img is None): 
				log.prNot(log.DEBUG, "Skipping %s, could not read file." % (base))
				continue
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
		log.prNot(log.INFO, "Done, saving results to disk.")
		allshifts = N.array(allshifts)
		if (self.file is not False):
			nf = libfile.saveData(self.file + '-shifts', allshifts, asnpy=True, \
			 	ascsv=True)

	


class StatsTool(Tool):
	"""Calculate statistics on files"""
	def __init__(self, files, params):
		super(StatsTool, self).__init__(files, params)
		self.run()
	
	
	def run(self):
		# Store all stats and all files
		allstats = []
		allfiles = []
		for f in self.files:
			base = os.path.basename(f)
			allfiles.append(base)
			# Load file
			img = self.load(f)
			if (img is None): 
				log.prNot(log.DEBUG, "Skipping %s, could not read file." % (base))
				continue
			# Dark-flat file if needed
			dfimg = self.darkflat(img)
			data = dfimg
			# If a maskfile is given, only calculate stats within the subaperture.
			if (self.maskfile is not 'none'):
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
				allstats.append(list(substat))
				log.prNot(log.INFO, "%s: mean: %.4g std: %.4g rms: %.4g (%.3g%%)" % \
					(base, substat[0], substat[1], substat[2], substat[3]))
			else:
				# If no maskfile given, calculate stats for all pixels
				r = (N.min(data), N.max(data))
				avg = N.mean(data)
				std = N.var(data)**0.5
				rms = (N.sum((data-avg)**2.0)/data.size)**0.5
				rmsrat = 100.0*rms/avg	
				log.prNot(log.INFO, "%s: mean: %.4g std: %.4g range: %.3g--%.3g rms: %.4g (%.3g%%)" % \
					(base, avg, std, r[0], r[1], rms, rmsrat))
				allstats.append([avg, std, rms, rmsrat])
		# Process stats for all files, display average
		allstats = N.array(allstats)
		allfiles = (N.array(allfiles)).reshape(-1,1)
		all_avg = N.mean(allstats, axis=0)
		all_std = (N.var(allstats, axis=0))**0.5
		log.prNot(log.INFO, "all %d: mean: %.4g+-%.4g rms: %.4g+-%.4g (%.3g%%+-%.4g)" % \
			(allstats.shape[0], all_avg[0], all_std[0], all_avg[2], all_std[2], all_avg[3], all_std[3]))
		# Save results if requested
		if (self.file is not False):
			nf = libfile.saveData(self.file, allstats, asnpy=True)
			hdr = ['filename, path=%s' % \
			 	(os.path.dirname(os.path.realpath(self.files[0]))), 'avg', \
			 	'std', 'rms', 'fractional rms']
			cf = libfile.saveData(self.file, N.concatenate((allfiles, allstats), \
			 	axis=1), ascsv=True, csvhdr=hdr, csvfmt='%s')

	


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
			log.prNot(log.INFO, "Convertintg file '%s' to %s" % (base, \
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
