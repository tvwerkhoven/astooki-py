#!/usr/bin/env python2.5
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

import sys, os, time
import getopt
import numpy as N
import scipy as S
import libsh
import libfile as lf
import liblog as log
#import libplot

VERSION = "0.0.2"
AUTHOR = "Tim van Werkhoven (tim@astro.su.se)"
DATE = "2009-04-24"

help_message = '''astooki version %s (%s) by %s.
Usage: astooki <TOOL> [OPTIONS] [FILES]

Tools
 convert                     Convert files to another format
 stats                       Get statistics on files
 shiftoverlay                Overlay image shifts on the raw images
 samask                      Make a subaperture mask
 sfmask                      Make a subfield mask
 saopt                       Optimze a subaperture mask with a flat
 saupd                       Update a subaperture mask with an offset
 shifts                      Measure image shifts in various subfields and 
                               subimages
 procshifts                  Process shifts measured with the shifts tool

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
 -l, --log=FILE              log messages to this file as well
 -h, --help                  show this help
 -d, --dir=DIR               use DIR as outputdir for all files
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

Shiftoverlay options
     --intclip=LOW,HIGH      clip intensity to this range
     --subap=SA1,SA2,...     subapertures to process, set to -1 for all.
     --safile=FILEPATH       subaperture positions
     --sffile=FILEPATH       subfield positions relative to subap (csv)
     --shifts=FILEPATH       image shifts file
     --shape=[box,dot]       shape to indicate the shifts with
     --skip=INT              skip this many entries in shifts file

Samask options
 -f, --file=FILEPATH         file to store subaperture configuration to
     --rad=RAD               radius of the subaperture pattern [1024]
     --shape=SHAPE           pattern shape (circular or square) [circular]
     --sasize=X,Y            size of the subapertures
     --pitch=X,Y             pitch of the subapertures
     --xoff=EVENOFF,ODDOFF   x-offset for subapertures in even, odd rows in
                               units of --sasize [0, 0.5]
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
   astooki.py saupd -vv --plot --mf 2009.04.22-mask.csv --offset \\
   offset-csv.csv

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
 
 Overlay shift vectors on raw data:
 N.B. Make sure the order of files correspond with the order of shifts, i.e. 
 entry N in the shifts file must correspond with file N in the file list!
   astooki.py shiftoverlay -vv --shape box --scale 1.471 --outformat png \\
   --intclip=0.85,1.15 --subap 61  --safile \\
   ../2009.04.28-run01/proc/2009.04.28-mask-updated.csv --sffile \\
   ../2009.04.28-run01/proc/2009.04.28-subfield-24x24.csv --shifts \\
   proc/subshift2/2009.04.28-run01_wfwfs_test_im28Apr2009_3-1002-shifts.npy \\
   --ff ../2009.04.28-flats/wfwfs_test_ff28Apr2009.0000002 --fm 500 --df \\
   ../2009.04.28-darks/wfwfs_test_dd28Apr2009.0000002 --dm 500 \\
   wfwfs_test_im28Apr2009.0000*
''' % (VERSION, DATE, AUTHOR)

### Supported formats
_FORMAT_ANA = 'ana'
_FORMAT_FITS = 'fits'
_FORMAT_PNG = 'png'
_FORMAT_NPY = 'npy'
_INFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_NPY)
_OUTFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_PNG, _FORMAT_NPY)
### Tools available
_TOOLS = ('convert', 'stats', 'shiftoverlay', 'samask', 'sfmask', 'saopt', 'saupd', 'shifts', 'procshifts')
### Default types to use
_ftype = N.float64
_itype = N.int32

### ==========================================================================
### Startup functions
### ==========================================================================

def main(argv=None):
	beg = time.time()
	# Parse command-line options from sys.argv
	(tool, params, files) = parse_options()
	# Sanity check on parameters
	check_params(tool, params)
	log.prNot(log.NOTICE, "Parameters: '%s'" % (str(params)))
	log.prNot(log.NOTICE, "Files prefix: '%s'" % (os.path.commonprefix(files)))
	# Perform action requested
	log.prNot(log.NOTICE, "Tool: %s." % tool)
	if (tool == 'convert'): ConvertTool(files,params)
	elif (tool == 'stats'): StatsTool(files, params)
	elif (tool == 'shiftoverlay'): ShiftOverlayTool(files, params)
	elif (tool == 'samask'): SubaptConfTool(files, params)
	elif (tool == 'sfmask'): SubfieldConfTool(files, params)
	elif (tool == 'saopt'): SubaptOptTool(files, params)
	elif (tool == 'saupd'): SubaptUpdateTool(files, params)
	elif (tool == 'shifts'): ShiftTool(files, params)
	elif (tool == 'procshifts'): ProcShiftsTool(files, params)
	
	# done
	dur = time.time() - beg
	log.prNot(log.NOTICE, "Complete, processed %d files in %d seconds (%.3gfps)" % (len(files), dur, len(files)/dur))
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
	
	opts, args = getopt.getopt(argv[2:], "vhsi:o:d:f:r:n:l:", ["verbose", "help", "stats", "informat=", "ff=", "fm=", "df=", "dm=", "mf=", "outformat=", "intclip=", "crop=", "file=", "dir=", "scale=", "rad=", "shape=", "pitch=", "xoff=", "disp=", "plot", "noplot", "norm", "nonorm", "saifac=", "range=", "sffile=", "safile=", "nref=", "sfsize=", "sasize=", "overlap=", "border=", "offsets=", "subap=", "shifts=", "log=", "skip="])
	# Remaining 'args' must be files
	files = args
	
	for option, value in opts:
		log.prNot(log.INFO, 'Parsing: %s:%s' % (option, value))
		if option in ["-v", "--verbose"]: log.VERBOSITY += 1
		if option in ["-l", "--log"]: params['logfile'] = os.path.realpath(value)
	for option, value in opts:
		log.prNot(log.INFO, 'Parsing: %s:%s' % (option, value))
		if option in ["-d", "--dir"]: params['outdir'] = os.path.realpath(value)
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
		if option in ["-f", "--file"]: params['file'] = value
		if option in ["--crop"]: params['crop'] = value.split(',')
		# Convert
		if option in ["-o", "--outformat"]: params['outformat'] = value
		if option in ["--scale"]: params['scale'] = float(value)
		if option in ["--intclip"]: params['intclip'] = value.split(',')
		# ShiftOverlay
		if option in ["--subap"]: params['subap'] = value.split(',')
		if option in ["--shifts"]: params['shifts'] = os.path.realpath(value)
		if option in ["--shape"]: params['shape'] = value
		if option in ["--skip"]: params['skip'] = int(value)
		# Samask
		if option in ["--rad"]: params['rad'] = float(value)
		if option in ["--shape"]: params['shape'] = value
		#if option in ["--size"]: params['size'] = value.split(',')
		if option in ["--pitch"]: params['pitch'] = value.split(',')
		if option in ["--xoff"]: params['xoff'] = value.split(',')
		if option in ["--disp"]: params['disp'] = value.split(',')
		if option in ["--scale"]: params['scale'] = float(value)
		# Sfmask
		if option in ["--sfsize"]: params['sfsize'] = value.split(',')
		if option in ["--sasize"]: params['sasize'] = value.split(',')
		if option in ["--overlap"]: params['overlap'] = value.split(',')
		if option in ["--border"]: params['border'] = value.split(',')
		# Saopt
		if option in ["--saifac"]: params['saifac'] = float(value)
		#if option in ["--rad"]: params['rad'] = float(value)
		# Saupd
		if option in ["--offsets"]: params['offsets'] = os.path.realpath(value)
		# Shifts
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
	default['logfile'] = 'astooki-log'
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
	default['outdir'] = os.path.realpath('astooki-out/')
	default['crop'] = False
	if (tool == 'convert'):
		default['scale'] = 1.0
		default['intclip'] = False
		default['outformat'] = _FORMAT_FITS
	elif (tool == 'shiftoverlay'):
		default['subap'] = ['-1']
		default['shifts'] = False
		default['shape'] = 'box'
		default['skip'] = 0
	elif (tool == 'samask'):
		default['rad'] = 1024
		default['shape'] = 'circular'
		default['sasize'] = ['128','128']
		default['pitch'] = ['164','164']
		default['xoff'] = ['0', '0.5']
		default['disp'] = ['0','0']
	elif (tool == 'sfmask'):
		default['sfsize'] = ['16', '16']
		default['sasize'] = ['0', '0']
		default['overlap'] = ['0.5', '0.5']
		default['border'] = ['6', '6']
	elif (tool == 'saopt'):
		default['saifac'] = 0.7
	elif (tool == 'saupd'):
		default['offsets'] = False
	elif (tool == 'shifts'):
		default['shrange'] = 7
		default['safile'] = False
		default['sffile'] = False
		default['nref'] = 4
	
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
		
	# First fix log-related things
	# ===================================================================
	if (not os.path.isdir(params['outdir'])): os.makedirs(params['outdir'])
	log.initLogFile(os.path.join(params['outdir'], params['logfile']))
	
	# Common options requirements
	# ===================================================================
	if params['informat'] not in _INFORMATS:
		log.prNot(log.ERR, "Unsupported input format '%s'" % (params['informat']))
	
	if (params.has_key('intclip') and params['intclip'] is not False):
		try: params['intclip'] = N.array(params['intclip']).astype(N.float)[[0,1]]
		except: log.prNot(log.ERR, "intclip invalid, should be <float>,<float>.")
	
	if (params.has_key('crop') and params['crop'] is not False):
		try: params['crop'] = N.array(params['crop']).astype(N.float)[[0,1,2,3]]
		except: log.prNot(log.ERR, "crop invalid, should be 4 floats.")
	
	if (params['flatfield']) and (not os.path.exists(params['flatfield'])):
			log.prNot(log.ERR, "flatfield '%s' does not exist." % \
		 		(params['flatfield']))
	
	if (params['darkfield']) and (not os.path.exists(params['darkfield'])):
		log.prNot(log.ERR, "darkfield '%s' does not exist." % \
		 	(params['darkfield']))
	
	if (params['maskfile']) and (not os.path.exists(params['maskfile'])):
		log.prNot(log.ERR, "maskfile '%s' does not exist." % \
		 	(params['maskfile']))
	
	if (params.has_key('file') and params['file']):
		lf.saveOldFile(params['file'], postfix='.old', maxold=5)
	
	if params.has_key('pitch'):
		try: params['pitch'] = N.array(params['pitch']).astype(N.int)[[0,1]]
		except: log.prNot(log.ERR, "pitch invalid, should be <int>,<int>.")
	
	if params.has_key('xoff'):
		try: params['xoff'] = N.array(params['xoff']).astype(N.float)[[0,1]]
		except: log.prNot(log.ERR, "xoff invalid, should be <float>,<float>.")
	
	if params.has_key('disp'):
		try: params['disp'] = N.array(params['disp']).astype(N.int)[[0,1]]
		except: log.prNot(log.ERR, "disp invalid, should be <int>,<int>.")
	
	if params.has_key('sfsize'):
		try: params['sfsize'] = N.array(params['sfsize']).astype(N.int)[[0,1]]
		except: log.prNot(log.ERR, "sfsize invalid, should be <int>,<int>.")
	
	if params.has_key('sasize'):
		try: params['sasize'] = N.array(params['sasize']).astype(N.int)[[0,1]]
		except: log.prNot(log.ERR, "sasize invalid, should be <int>,<int>.")
	
	if params.has_key('overlap'):
		try: params['overlap'] = N.array(params['overlap']).astype(N.float)[[0,1]]
		except: log.prNot(log.ERR, "overlap invalid, should be <float>,<float>.")
	
	if params.has_key('border'):
		try: params['border'] = N.array(params['border']).astype(N.int)[[0,1]]
		except: log.prNot(log.ERR, "border invalid, should be <int>,<int>.")
	
	if params.has_key('subap'):
		try: params['subap'] = N.array(params['subap']).astype(N.int)
		except: log.prNot(log.ERR, "subap invalid, should be list of ints.")
	
	if params.has_key('offsets'):	
		if params['offsets'] and not os.path.exists(params['offsets']):
			log.prNot(log.ERR, "offset file '%s' does not exist." % \
		 		(params['offsets']))
	
	if params.has_key('outformat') and (params['outformat'] not in _OUTFORMATS):
		log.prNot(log.ERR, "Unsupported output '%s'" % (params['outformat']))
	
	# Requirements depending on tools (where defaults are not sufficient)
	# ===================================================================
	if (tool == 'saopt'):
		# need flatfield, maskfile
		if (not params['flatfield']) or (not os.path.exists(params['flatfield'])):
			log.prNot(log.ERR, "Tool 'saopt' requires flatfield.")
		if (not params['maskfile']) or (not os.path.exists(params['maskfile'])):
			log.prNot(log.ERR, "Tool 'saopt' requires maskfile.")
	elif (tool == 'shiftoverlay'):
		if (params['scale'] <= 1.0):
			log.prNot(log.WARNING, "Recommend to set scale to higher than 1.0")
		if (not params['sffile']) or (not os.path.exists(params['sffile'])):
			log.prNot(log.ERR, "Tool 'shiftoverlay' requires sffile.")
		if (not params['safile']) or (not os.path.exists(params['safile'])):
			log.prNot(log.ERR, "Tool 'shiftoverlay' requires safile.")
		if (not params['shifts']) or (not os.path.exists(params['shifts'])):
			log.prNot(log.ERR, "Tool 'shiftoverlay' requires shifts file.")
		if (params['shape'] not in ['box', 'dot']):
			log.prNot(log.ERR, "'shape' should be in ['box', 'dot'].")
	elif (tool == 'shifts'):
		# need safile and sffile
		if (not params['safile']) or (not os.path.exists(params['safile'])):
			log.prNot(log.ERR, "Tool 'shifts' requires safile.")
		if (not params['sffile']) or (not os.path.exists(params['sffile'])):
			log.prNot(log.ERR, "Tool 'shifts' requires sffile.")
	elif (tool == 'samask'):
		# Shape should be 'circular' or 'square'
		if (params['shape'] not in ['square', 'circular']):
			log.prNot(log.ERR, "shape invalid, should be 'circular' or 'square'")
	elif (tool == 'sfmask'):
		# sasize > sfsize, both must be int, int
		if (params['sasize'] < params['sfsize']).any():
			log.prNot(log.ERR, "sasize must be bigger than sfsize.")
		if (params['overlap'] > 1).any() or (params['overlap'] < 0).any():
			log.prNot(log.ERR, "overlap must be between 0 and 1.")
	elif (tool == 'saupd'):
		# Saupd needs maskfile, offset file
		if (params['maskfile']) and \
			(not os.path.exists(params['maskfile'])):
			log.prNot(log.ERR, "Tool 'saupd' requires maskfile.")
		if (params['offsets']) and \
			(not os.path.exists(params['offsets'])):
			log.prNot(log.ERR, "Tool 'saupd' requires offsets file.")
	
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
		self.outdir = params['outdir']
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
			log.prNot(log.INFO, "Loading mask files.")
			(self.nsa, self.saccdpos, self.saccdsize) = \
			 	libsh.loadSaSfConf(self.maskfile)
		# Save all files used to outdir
		for (key, val) in params.items():
			if (val.__class__ == 'string'.__class__):
				if os.path.isfile(val):
					f = os.path.basename(val)
					uri = os.path.join(self.outdir, f)
					lf.saveOldFile(uri)
					os.link(val, uri)
		log.prNot(log.NOTICE, "Processing %d files" % (len(self.files)))
	
	
	def load(self, filename):
		log.prNot(log.INFO, "Loading '%s'" % (filename))
		if not os.path.exists(filename):
			log.prNot(log.WARNING, "File '%s' does not exist" % (filename))
			return None
		# Load data, using different methods depending on type
		if (self.informat == _FORMAT_ANA):
			data = self.__anaload(filename)
		elif (self.informat == _FORMAT_FITS):
			data = self.__fitsload(filename)	
		elif (self.informat == _FORMAT_NPY):
			data = self.__npyload(filename)	
		else:
			log.prNot(log.WARNING, "Filetype unsupported." % (filename))			
			return None
		self.origres = data.shape
		if (self.crop is not False):
			data = data[self.crop[1]:self.crop[1] + self.crop[3], \
				self.crop[0]:self.crop[0] + self.crop[2]]

		return data.astype(_ftype)
	
	
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
			log.prNot(log.INFO, "Dark-flatfielding data. Dark avg: %.4g, gain avg: %.4g" % (N.mean(self.darkdata), N.mean(self.gaindata)))		
			# Now process the frame
			return (data-self.darkdata) * self.gaindata
		return data
	
	
	def __initdarkflat(self):
		# Get flats and darks, if not already present
		if (self.flatdata is None):
			if (self.flatfield):
				log.prNot(log.INFO, "Loading flatfield...")
				self.flatdata = self.load(self.flatfield)
				self.flatdata /= 1.0*self.flatmulti
				log.prNot(log.INFO, "Flatfield average: %.6g" % N.mean(self.flatdata))
			else:
				# Maybe we don't want flatfielding, in that case set it to 1.0
				log.prNot(log.INFO, "Not flatfielding, setting to 1.0")
				self.gaindata = 1.0
		if (self.darkdata is None):
			if (self.darkfield):
				log.prNot(log.INFO, "Loading darkfield...")
				self.darkdata = self.load(self.darkfield)
				self.darkdata /= 1.0*self.darkmulti
				log.prNot(log.INFO, "Darkfield average: %.6g" % N.mean(self.darkdata))
			else:
				# Maybe we don't want darkfielding, in that case set it to 1.0
				log.prNot(log.INFO, "Not darkfielding, setting to 0.0")
				self.darkdata = 0.0
		if (self.gaindata is None):
			# Make a gain for faster processing
			invgain = (self.flatdata - self.darkdata)
			# Prevent infinity
			invgain[invgain <= 0] = 1.0
			self.gaindata = 1.0/invgain
			#self.gaindata /= N.mean(self.gaindata)
	
	
	def fitssave(self, data, filepath, overwrite=True):
		"""
		Save 'data' as FITS file to 'filepath'.
		"""
		import pyfits
		log.prNot(log.INFO, "Tool.fitssave(): Saving data to '%s'." % (filepath))
		pyfits.writeto(filepath, data, clobber=overwrite)
	
	
	def anasave(self, data, filepath, compressed=1):
		"""
		Save 'data' as ANA file to 'filepath'. Can be compressed (default: yes).
		"""
		import pyana
		log.prNot(log.INFO, "Tool.anasave(): Saving data to '%s'." % (filepath))
		pyana.fzwrite(filepath, data, compressed)
	
	
	def npysave(self, data, filepath):
		"""
		Save 'data' as npy file to 'filepath'.
		"""
		log.prNot(log.INFO, "Tool.npysave(): Saving data to '%s'." % (filepath))
		N.save(data, filepath)
	
	
	def pngsave(self, data, filepath, scale=True):
		"""
		Save 'data' as PNG file to 'filepath'. Data can be scaled to full range
		(default: yes).
		"""
		log.prNot(log.INFO, "Tool.pngsave(): Saving data to '%s'." % (filepath))
		if (scale):
			# Scale the values to 0-255
			maxval = N.max(data)
			minval = N.min(data)
			scdata = (255.0*(data - minval)/(maxval - minval))
			scdata = (scdata.astype(N.uint8))
		else: scdata = data
		if (scdata.shape[1] % 4 != 0): 
			raise RuntimeError("Cannot save PNG with width not a multiple of 4 (cairo bug, width was %d)." % (scdata.shape[1]))
		
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
	
	
	def mkuri(self, path):
		"""
		"""
		_path = os.path.basename(path)
		if (_path != path):
			log.prNot(log.WARNING, "mkuri(): got path instead of filename.")
		
		return os.path.join(self.outdir, _path)
	
	
	def maskimg(self, data):
		"""
		Apply a mask on an image, set all values outside the mask to the minimum 
		value inside the mask, so it will appear black.
		"""
		if (self.maskfile is False):
			return data
		log.prNot(log.INFO, "Masking image if necessary, res: %d,%d" % \
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
				log.prNot(log.INFO, "maskimg(): normalizing, avg: %.3g" % (avg))
		
		return data
	
	
	def __initmask(self, res):
		if (self.mask is None) or (res != self.maskres):
			# We need to make a new mask here
			log.prNot(log.INFO, "maskimg(): (re-)initializing mask with resolution %d,%d" % (tuple(res)))
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
		self.sasize = params['sasize']
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
			libsh.calcSubaptConf(self.rad, self.sasize, self.pitch, self.shape, \
			self.xoff, self.disp, self.scale)
		# Save to file
		libsh.saveSaSfConf(self.file, nsa, [-1,-1], saccdsize, saccdpos)
		if (self.plot):
			import libplot
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
			N.indices(nsf, dtype=N.float).reshape(2,-1).T * effpitch
		sfpos = N.round(sfpos).astype(N.int32)
		totnsf = N.product(nsf).astype(N.int32)
		
		log.prNot(log.NOTICE, "Found %d x %d subfields." % tuple(nsf))
		log.prNot(log.NOTICE, "Size %d,%d" % tuple(self.sfsize))
		libsh.saveSaSfConf(self.file, totnsf, [-1,-1], self.sfsize, sfpos)
		if (self.plot):
			import libplot
			# TODO: fix this, splitting does not work correctly when files have 
			# multiple dots (2009.04.05-mask.csv)
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
		off = lf.loadData(self.offsets, ascsv=True)
		
		# Compare
		if (off.shape[0] != nsa):
			log.prNot(log.ERR, "SubaptUpdateTool(): offsets not the same size as subaperture positions.")
		
		# Offset the new positions
		maxsh = N.ceil(N.max(abs(off), axis=0))
		newpos = (pos + off + maxsh).astype(N.int32)
		# Crop the subaperture size by twice the maximum offset
		newsize = (size - maxsh*2).astype(N.int32)
		# Store
		libsh.saveSaSfConf(self.file, nsa, [-1,-1], newsize, newpos)
		if (self.plot):
			import libplot
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
			import libplot
			plfile = os.path.splitext(self.file)[0]+'-plot.eps'
			libplot.showSaSfLayout(plfile, opos, osize, \
				plrange=[[0, 2*self.rad]]*2)
		# Done
	


class ShiftTool(Tool):
	"""Calculate shifts between different subfields and subapertures."""
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
		allrefs = []
		for f in self.files:
			base = os.path.basename(f)
			log.prNot(log.NOTICE, "Measuring shifts for %s." % (base))
			# Load file
			img = self.load(f)
			if (img is None): 
				log.prNot(log.INFO, "Skipping %s, could not read file." % (base))
				continue
			allfiles.append(base)
			# Dark-flat file if needed
			dfimg = self.darkflat(img)
			
			# Measure shift
			refaps = []
			imgshifts = ls.calcShifts(dfimg, self.saccdpos, self.saccdsize, \
			 	self.sfccdpos, self.sfccdsize, method=ls.COMPARE_ABSDIFFSQ, \
			 	extremum=ls.EXTREMUM_2D9PTSQ, refmode=ls.REF_BESTRMS, \
			 	refopt=self.nref, shrange=[self.shrange, self.shrange], \
			 	subfields=None, corrmaps=None, refaps=refaps)
			allrefs.append(refaps)
			allshifts.append(imgshifts)
		
		# Process results, store to disk
		log.prNot(log.NOTICE, "Done, saving results to disk @ '%s'." % (self.file))
		allshifts = N.array(allshifts)
		# Store the list of files where we save data to
		files = {}
		files['shifts'] = lf.saveData(self.mkuri(self.file + '-shifts'), \
		 	allshifts, asnpy=True, asfits=True)
		files['refaps'] = lf.saveData(self.mkuri(self.file + '-refaps'), \
			allrefs, asnpy=True, ascsv=True)
		files['saccdpos'] = lf.saveData(self.mkuri(self.file + '-saccdpos'), \
			self.saccdpos, asnpy=True)
		files['sfccdpos'] = lf.saveData(self.mkuri(self.file + '-sfccdpos'), \
		 	self.sfccdpos, asnpy=True)
		files['saccdsize'] = lf.saveData(self.mkuri(self.file + '-saccdsize'), \
		 	self.saccdsize, asnpy=True)
		files['sfccdsize'] = lf.saveData(self.mkuri(self.file + '-sfccdsize'), \
		 	self.sfccdsize, asnpy=True)
		cpos = self.saccdpos.reshape(-1,1,2) + self.sfccdpos.reshape(1,-1,2) + \
			self.sfccdsize.reshape(1,1,2)/2.0
		files['sasfpos-c'] = lf.saveData(self.mkuri(self.file + '-sasfpos-c'), \
		 	cpos, asnpy=True, asfits=True)
		files['files'] = lf.saveData(self.mkuri(self.file + '-files'), \
		 	allfiles, asnpy=True, ascsv=True, csvfmt='%s')
		# Add meta info
		files['path'] = os.path.dirname(os.path.realpath(self.file))
		files['base'] = os.path.basename(self.file)
		
		metafile = lf.saveData(self.mkuri(self.file + '-meta'), \
			files, aspickle=True)
		# If we have only one subfield, also calculate 'static' shifts


	


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
			# Load file
			img = self.load(f)
			if (img is None): 
				log.prNot(log.INFO, "Skipping %s, could not read file." % (base))
				return False
			allfiles.append(base)
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
				log.prNot(log.NOTICE, "%s: mean: %.4g std: %.4g rms: %.4g (%.3g%%)" % \
					(base, substat[0], substat[1], substat[2], substat[3]))
				allstats.append(list(substat))
			else:
				# If no maskfile given, calculate stats for all pixels
				r = (N.min(data), N.max(data))
				avg = N.mean(data)
				std = N.var(data)**0.5
				rms = (N.sum((data-avg)**2.0)/data.size)**0.5
				rmsrat = 100.0*rms/avg	
				log.prNot(log.NOTICE, "%s: mean: %.4g std: %.4g range: %.3g--%.3g rms: %.4g (%.3g%%)" % \
					(base, avg, std, r[0], r[1], rms, rmsrat))
				allstats.append([avg, std, rms, rmsrat])
				
		
		# Process stats for all files, display average
		allstats = N.array(allstats)
		allfiles = (N.array(allfiles)).reshape(-1,1)
		all_avg = N.mean(allstats, axis=0)
		all_std = (N.var(allstats, axis=0))**0.5
		log.prNot(log.NOTICE, "all %d: mean: %.4g+-%.4g rms: %.4g+-%.4g (%.3g%%+-%.4g)" % \
			(allstats.shape[0], all_avg[0], all_std[0], all_avg[2], all_std[2], all_avg[3], all_std[3]))
		# Save results if requested
		nf = lf.saveData(self.mkuri(self.file + '-stats'), allstats, asnpy=True)
		hdr = ['filename, path=%s' % \
		 	(os.path.dirname(os.path.realpath(self.files[0]))), 'avg', \
		 	'std', 'rms', 'fractional rms']
		cf = lf.saveData(self.mkuri(self.file + '-stats'), \
		 	N.concatenate((allfiles, allstats), axis=1), ascsv=True, csvhdr=hdr, \
		 	csvfmt='%s')
	



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
			log.prNot(log.NOTICE, "Converting file '%s' to %s" % (base, \
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
				orig = data.shape
				nsc = (N.round(orig[1] * sc / 8) * 8)/orig[1]
				data = S.ndimage.zoom(data, nsc, mode='wrap')
				log.prNot(log.NOTICE, "Scaling image by %g (from %d,%d to %d,%d)." % \
				 	(nsc, orig[0], orig[1], data.shape[0], data.shape[1]))
			# Crop intensity if necessary
			if (self.intclip is not False):
				log.prNot(log.NOTICE, "Clipping intensity to %g--%g." % \
				 	tuple(self.intclip))
				data = N.clip(data, self.intclip[0], self.intclip[1])
			# Save again
			if (len(self.files) == 1 and self.file): 
				savefile = self.mkuri(self.file)
			else: 
				savefile = self.mkuri(f+'.'+self.outformat)
			log.prNot(log.NOTICE, "Saving '%s' as %s in '%s'." % \
				(base, self.outformat, os.path.basename(savefile)))
				
			if (self.outformat == _FORMAT_PNG): self.pngsave(data, savefile)
			elif (self.outformat == _FORMAT_FITS): self.fitssave(data, savefile)
			elif (self.outformat == _FORMAT_ANA): self.anasave(data, savefile)
			elif (self.outformat == _FORMAT_NPY): self.npysave(data, savefile)
		
	


class ShiftOverlayTool(Tool):
	"""Crop out one subaperture, overlay shift vectors, output PNGs"""
	def __init__(self, files, params):
		super(ShiftOverlayTool, self).__init__(files, params)
		self.file = params['file']
		self.scale = params['scale']
		self.intclip = params['intclip']
		self.outformat = params['outformat']
		self.shape = params['shape']
		self.skip = params['skip']
		sffile = params['sffile']
		safile = params['safile']
		(self.nsa, self.saccdpos, self.saccdsize) = \
			libsh.loadSaSfConf(safile)
		(self.nsf, self.sfccdpos, self.sfccdsize) = \
			libsh.loadSaSfConf(sffile)
		
		if params['subap'][0] == -1: self.subaps = range(self.nsa)
		else: self.subaps = params['subap']
		self.sfccdposc = self.sfccdpos + self.sfccdsize/2.0
		self.shifts = params['shifts']
		self.run()
	
	
	def run(self):
		# Do not plot overlapping subfields when using boxes
		if self.shape == 'box':
			plsf = []
			allpos = []
			sz = self.sfccdsize
			for sf in xrange(self.nsf):
				shpos = self.sfccdpos[sf] * self.scale
				# Check if subfields don't overlap
				if (sf != 0 and (abs(N.array(allpos)-shpos.reshape(1,2)) < sz*self.scale).all(axis=1).any()): continue
				oldpos = shpos			
				allpos.append(oldpos)
				plsf.append(sf)
			log.prNot(log.NOTICE, "Using %d subfields: %s" % (len(plsf), str(plsf)))
		elif self.shape == 'dot':
			plsf = xrange(self.nsf)
		
		allshifts = lf.loadData(self.shifts, asnpy=True)
		log.prNot(log.NOTICE, "Using subapertures %s" % (str(self.subaps)))
		
		# Loop over files to process
		for fidx in xrange(len(self.files)):
			f = self.files[fidx]			
			base = os.path.basename(f)
			img = self.load(f)
			if (img is None): continue
			dfimg = self.darkflat(img)
			log.prNot(log.NOTICE, "Processing file %d/%d, '%s'" % \
				(fidx+1, len(self.files), base))
			
			# Filter out non-finite values
			log.prNot(log.NOTICE, "Pre-processing shifts.")
			shifts = allshifts[fidx+self.skip,:,:,:,:]
			notfin = N.argwhere(N.isfinite(shifts) == False)
			for nfidx in notfin: shifts[tuple(nfidx)] = 0.0
			log.prNot(log.NOTICE, "Setting %d non-finite values to 0." % \
				(len(notfin)))
			
			# Average over different references
			shifts = N.mean(shifts, axis=0)
			# Subtract average over different subapertures
			#for sa in xrange(shifts.shape[0]):
			#	avg = N.mean(shifts[sa,:,:], axis=0)
			#	shifts[:,sf,:] -= avg.reshape(1,2)
			#print shifts.shape
			# Loop over subaps to process
			for sa in self.subaps:
				# self.crop = N.array([self.saccdpos[sa, 0], self.saccdpos[sa, 1], \
				# 	self.saccdsize[0], self.saccdsize[1]])
				# Crop out subaperture, divide by mean
				data = dfimg[\
					self.saccdpos[sa, 1]: \
					self.saccdpos[sa, 1] + self.saccdsize[1], \
					self.saccdpos[sa, 0]: \
					self.saccdpos[sa, 0] + self.saccdsize[0]]
				data = data / data.mean()
				
				if (self.scale != 1.0):
					import scipy.ndimage
					sc = self.scale
					orig = data.shape
					nsc = (N.round(orig[1] * sc / 8) * 8)/orig[1]
					data = S.ndimage.zoom(data, nsc, mode='wrap')
					log.prNot(log.NOTICE, "Scaling image by %g (from %d,%d to %d,%d)."%\
					 	(nsc, orig[1], orig[0], data.shape[1], data.shape[0]))
					nsc = (N.array(data.shape)*1.0/N.array(orig))[::-1]
				# Clip intensity if necessary
				if (self.intclip is not False):
					log.prNot(log.NOTICE, "Clipping intensity to %g--%g." % \
					 	tuple(self.intclip))
					data = N.clip(data, self.intclip[0], self.intclip[1])
				# Add shift vectors, set value at shift vector to max
				dmax = data.max()
				dmin = data.min()
				for sf in plsf:
					if self.shape == 'box':
						shpos = (self.sfccdpos[sf] - shifts[sa, sf]) * nsc
						data[shpos[1]:shpos[1]+(sz[1]*nsc[1]), \
							shpos[0]] = dmax
						data[shpos[1]:shpos[1]+(sz[1]*nsc[1]), \
							shpos[0]+(sz[0]*nsc[0])-1] = dmax
						data[shpos[1]+(sz[1]*nsc[1])-1, \
							shpos[0]:shpos[0]+(sz[0]*nsc[0])] = dmax
						data[shpos[1], \
							shpos[0]:shpos[0]+(sz[0]*nsc[1])] = dmax
					elif self.shape == 'dot':
					# This gives a 9 pixel block:
						vec = (self.sfccdposc[sf] - shifts[sa, sf]) * nsc
						data[vec[1]-1:vec[1]+2, vec[0]-1:vec[0]+2] = dmax
						data[vec[1], vec[0]] = dmin
				
				# Save again
				savefile = f+'-subap%d-shifts.%s' % (sa, self.outformat)
				savefile = self.mkuri(savefile)
				log.prNot(log.NOTICE, "Saving '%s' to '%s'." % \
					(base, os.path.basename(savefile)))
				
				if (self.outformat == _FORMAT_PNG): self.pngsave(data, savefile)
				elif (self.outformat == _FORMAT_FITS): self.fitssave(data, savefile)
				elif (self.outformat == _FORMAT_ANA): self.anasave(data, savefile)
				elif (self.outformat == _FORMAT_NPY): self.npysave(data, savefile)
		
	


class ProcShiftsTool(Tool):
	"""Process shift data"""
	def __init__(self, files, params):
		super(ProcShiftsTool, self).__init__(files, params)
		self.run()
	
	
	def run(self):
		for f in self.files:
			# Try to load pickle file
			log.prNot(log.NOTICE, "Processing meta file %s" % (f))
			(data, metafiles) = lf.restoreData(f)
			# Make sure we have the right data
			try: sfccdpos = data['sfccdpos']
			except: log.prNot(log.ERR, "'sfccdpos' not found in data.")
			try: saccdpos = data['saccdpos']
			except: log.prNot(log.ERR, "'saccdpos' not found in data.")
			try: sfccdsize = data['sfccdsize']
			except: log.prNot(log.ERR, "'sfccdsize' not found in data.")
			try: saccdsize = data['saccdsize']
			except: log.prNot(log.ERR, "'saccdsize' not found in data.")
			try: allshifts = data['shifts']
			except: log.prNot(log.ERR, "'allshifts' not found in data.")
			
			# Process NaNs and other non-finite numbers
			notfin = N.argwhere(N.isfinite(allshifts) == False)
			notfin_perc = notfin.shape[0]*100./allshifts.size
			log.prNot(log.NOTICE, "%d (%.2g%%) non-finite entries, spread over %d frames, %d subaps, %d subfields." % \
			 	(notfin.shape[0], notfin_perc, N.unique(notfin[:,0]).size, \
			 	N.unique(notfin[:,2]).size, N.unique(notfin[:,3]).size))
			#log.prNot(log.NOTICE, "Worst frame: %d with %d non-finite entries." % \
			#	N.bincount(notfin[:,0])
			if (notfin_perc > 0.5):
				log.prNot(log.WARNING, "Percentage of non-finite entries very high!")
			metafiles['notfinite'] = lf.saveData(self.mkuri(data['base'] + \
			 	'-notfinite'), notfin, ascsv=True, asnpy=True)
			
			# Make a list of all finite frames by excluding all non-finite frames
			# NB: does not work properly, throws away too much data. How to repair
			# NaNs? Why do we get NaNs in the first place?
			finframes = range(allshifts.shape[0])
			for i in (N.unique(notfin[:,0]))[::-1]: finframes.pop(i)
			allshifts_fin = allshifts[finframes]
			
			# If we have one subfield, treat it as static shift data:
			if (len(sfccdpos) == 1):
				log.prNot(log.NOTICE, "Calculating static offsets.")
				(soff, sofferr) = libsh.procStatShift(allshifts_fin[:,:,:,0,:])
				metafiles['offsets'] = lf.saveData(self.mkuri(data['base'] + \
				 	'-offset'), soff, asnpy=True, ascsv=True)
				metafiles['offset-err'] = lf.saveData(self.mkuri(data['base'] + \
				 	'-offset-err'), sofferr, asnpy=True, ascsv=True)
				if (self.plot):
					import libplot
					libplot.plotShifts(data['base'] + '-offset-plot', allshifts_fin, \
						saccdpos, saccdsize, sfccdpos, sfccdsize, \
						plorigin=(0,0), plrange=(2048, 2048), mag=7.0, allsh=False, \
					 	title='Static offsets for' + data['base'] +', mag=7', legend=True)	
			# We have subfield data here, process
			else:
				pass
			
			# Store the new meta file
			lf.saveData(f, metafiles, aspickle=True, explicit=True)
			
	

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
