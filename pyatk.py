#!/usr/bin/env python2.5
# encoding: utf-8
"""
@filename pyatk.py
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
import astooki.libfile as lf
import astooki.liblog as log
import astooki.libsh as libsh
import getopt
import numpy as N
import scipy as S

VERSION = "0.0.2"
AUTHOR = "Tim van Werkhoven (tim@astro.su.se)"
DATE = "2009-04-24"

help_message = {}
help_message['preamble'] = """astooki version %s (%s) by %s.
Usage: pyatk <TOOL> [OPTIONS] [FILES]""" % (VERSION, DATE, AUTHOR)

help_message['syntaxerr'] = """pyatk.py: Syntax incorrect.
Usage: 
   pyatk <TOOL> [OPTIONS] [FILES]
or
   pyatk [TOOL] --help
for help."""

help_message['common'] = """Tools
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
 tomo                        Tomographically analyze differential image shifts
 sdimm                       SDIMM+ analysis of image shifts

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
     --config=FILE           python file holding configuration, instead of 
                               using command-line arguments [None]"""

help_message['convert'] = """Convert options
 -f, --file=FILEPATH         if single file, store converted file here 
     --scale=FACTOR          scale output resolution by this factor [1.0]
 -o, --outformat=FORMAT      output file-format [fits]
     --intclip=LOW,HIGH      clip intensity to this range"""

help_message['stats'] = """Stats options
 -f, --file=FILEPATH         file to store statistics to"""

help_message['shiftoverlay'] = """Shiftoverlay options
     --intclip=LOW,HIGH      clip intensity to this range
     --subap=SA1,SA2,...     subapertures to process, set to -1 for all.
     --safile=FILEPATH       subaperture positions (csv)
     --sffile=FILEPATH       subfield positions relative to subap (csv)
     --shifts=FILEPATH       image shifts file
     --shape=[box,dot]       shape to indicate the shifts with
     --skip=INT              skip this many entries in shifts file"""

help_message['samask'] = """Samask options
 -f, --file=FILEPATH         file to store subaperture configuration to
     --rad=RAD               radius of the subaperture pattern [1024]
     --shape=SHAPE           pattern shape (circular or square) [circular]
     --sasize=X,Y            size of the subapertures
     --pitch=X,Y             pitch of the subapertures
     --xoff=EVENOFF,ODDOFF   x-offset for subapertures in even, odd rows in
                               units of --sasize [0, 0.5]
     --disp=X,Y              global pattern offset [0,0]
     --scale=SCALE           global pattern scaling factor [1.0]
     --[no]plot              make a plot of the subaperture mask [no]"""

help_message['sfmask'] = """Sfmask options
 -f, --file=FILEPATH         file to store subaperture configuration to
     --sfsize=X,Y            subfield size [16, 16]
     --sasize=X,Y            subaperture size to fit things in
     --overlap=X,Y           how many overlap to allow between subfields in X 
                               and Y direction [0.5, 0.5]
     --border=X,Y            add a border around the subfields"""

help_message['saopt'] = """Saopt options
 -f, --file=FILEPATH         file to store optimized configuration to
     --saifac=FLOAT          intensity drop-off factor considered 'dark' [0.7]
     --rad=RAD               radius of the subaperture pattern, used for 
                               plotting [1024]"""
help_message['saupd'] = """Saupd options
     --mf=FILEPATH           mask to update
     --offsets=FILE          file holding offset vectors for all subapertures"""

help_message['shifts'] = """Shifts options
 -r, --range=INT             shift range to use for cross-correlation [7]
     --safile=FILEPATH       subaperture locations, same format as maskfile
     --sffile=FILEPATH       subfield locations w.r.t. subaperture locations
 -n, --nref=INT              number of references to use [4]"""

help_message['procshifts'] = """Procshifts options
     --safile=FILEPATH       subaperture locations, same format as maskfile
     --sffile=FILEPATH       subfield locations w.r.t. subaperture locations
     --shifts=FILE           shift measurements"""

help_message['tomo'] = """Tomo options
     --shifts=FILE           shift measurements
     --safile=FILE           centroid positions for each subaperture IN 
                               APERTURE SPACE (meter)
     --sffile=FILE           subfield positions on the CCD (pixels)
     --ccdres                angular resolution of the CCD pixels (arcsec/pix)
     --aptr=R                telescope aperture radius (meter)
     --nheights=N1,N2,...    number of heights each layer should be put at
     --layerheights=L1min-L1max,L2min-L2max,...
                             height range of each layer (meters)
     --layercells=W,H        number of cells in a layer"""

help_message['sdimm'] = """SDIMM options
     --shifts=FILE           shift measurements
     --safile=FILE           centroid positions for each subaperture IN 
                               APERTURE SPACE (meter)
     --sffile=FILE           subfield positions on the CCD (pixels)
     --skipsa=SA1,SA2,...    list of subapertures to skip in analysis"""

help_message['examples'] = """Examples
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

 SDIMM+ analysis on shift data:
   pyatk.py sdimm -vv -d sdimm-16x16-dc --shifts 
   subshift-16x16/2009.04.28-run05_wfwfs_test_im28Apr2009_1-999-shifts.npy 
   --safile samask/2009.04.28-run05-samask-ll-centroid.csv --sffile 
   sfmask/2009.04.28-run05-sfmask-16x16.csv"""


### Supported formats
_FORMAT_ANA = 'ana'
_FORMAT_FITS = 'fits'
_FORMAT_PNG = 'png'
_FORMAT_NPY = 'npy'
_INFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_NPY)
_OUTFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_PNG, _FORMAT_NPY)
### Tools available
_TOOLS = ('convert', 'stats', 'shiftoverlay', 'samask', 'sfmask', 'saopt', 'saupd', 'shifts', 'procshifts', 'tomo', 'sdimm')
### Default types to use
_ftype = N.float64
_itype = N.int32

### ==========================================================================
### Startup functions
### ==========================================================================

def main(argv=None):
	#print "Sleeping, attach debuggers now."
	#time.sleep(10)
	beg = time.time()
	# Parse command-line options from sys.argv
	(tool, params, files) = parse_options()
	# Sanity check on parameters
	check_params(tool, params)
	log.prNot(log.NOTICE, "Command: '%s'" % (str(sys.argv)))
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
	elif (tool == 'tomo'): TomoTool(files, params)
	elif (tool == 'sdimm'): SdimmTool(files, params)
	# done
	dur = time.time() - beg
	log.prNot(log.NOTICE, "Completed in %g seconds." % (dur))
	return 0


def parse_options():
	"""
	Process command-line options. See help_message for possible options.
	"""
	
	argv = sys.argv
	# First check whether argv[1] is present (could be -h, --help or a tool)
	try:
		tool = argv[1]
		if tool in ["-h", "--help"]: print_help('common')
		elif tool not in _TOOLS:
			raise Exception
	except:
		print_help('syntaxerr', out=sys.stderr)
		sys.exit(2)
	
	# Parse common options first
	# ==========================
	
	params = get_defaults(tool)
	
	opts, args = getopt.getopt(argv[2:], "vhsi:o:d:f:r:n:l:", ["verbose", "help", "stats", "informat=", "ff=", "fm=", "df=", "dm=", "mf=", "outformat=", "intclip=", "crop=", "file=", "dir=", "scale=", "rad=", "shape=", "pitch=", "xoff=", "origin", "noorigin", "disp=", "plot", "noplot", "norm", "nonorm", "saifac=", "range=", "sffile=", "safile=", "nref=", "sfsize=", "sasize=", "overlap=", "border=", "offsets=", "subap=", "shifts=", "log=", "skip=", "ccdres=", "aptr=", "layerheights=", "nheights=", "layercells=", "skipsa="])
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
		# Tomo / sdimm
		if option in ["--ccdres"]: params['ccdres'] = float(value)
		if option in ["--aptr"]: params['aptr'] = float(value)
		if option in ["--nheights"]: params['nheights'] = value.split(',')
		if option in ["--layerheights"]: params['layerheights'] = value.split(',')
		if option in ["--layercells"]: params['layercells'] = value.split(',')		
		# sdimm options
		if option in ["--skipsa"]: params['skipsa'] = value.split(',')		
	
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
	default['scale'] = 1.0
	if (tool == 'convert'):
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
	elif (tool == 'tomo'):
		default['ccdres'] = 0.45
		default['aptr'] = 0.49
		default['shifts'] = False
		default['nheights'] = ['3','1']
		default['layerheights'] = ['0-1000', '10000-10000']
		default['layercells'] = ['8', '8']
	
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
		try: params['intclip'] = N.array(params['intclip'], dtype=N.float)[[0,1]]
		except: log.prNot(log.ERR, "intclip invalid, should be <float>,<float>.")
	
	if (params.has_key('crop') and params['crop'] is not False):
		try: params['crop'] = N.array(params['crop'], dtype=N.float)[[0,1,2,3]]
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
		try: params['pitch'] = N.array(params['pitch'], dtype=N.float)[[0,1]]
		except: log.prNot(log.ERR, "pitch invalid, should be <float>,<float>.")
	
	if params.has_key('xoff'):
		try: params['xoff'] = N.array(params['xoff'], dtype=N.float)[[0,1]]
		except: log.prNot(log.ERR, "xoff invalid, should be <float>,<float>.")
	
	if params.has_key('disp'):
		try: params['disp'] = N.array(params['disp'], dtype=N.float)[[0,1]]
		except: log.prNot(log.ERR, "disp invalid, should be <float>,<float>.")
	
	if params.has_key('sfsize'):
		try: params['sfsize'] = N.array(params['sfsize'], dtype=N.int)[[0,1]]
		except: log.prNot(log.ERR, "sfsize invalid, should be <int>,<int>.")
	
	if params.has_key('sasize'):
		try: params['sasize'] = N.array(params['sasize'], dtype=N.float)[[0,1]]
		except: log.prNot(log.ERR, "sasize invalid, should be <float>,<float>.")
	
	if params.has_key('overlap'):
		try: params['overlap'] = N.array(params['overlap'], dtype=N.float)[[0,1]]
		except: log.prNot(log.ERR, "overlap invalid, should be <float>,<float>.")
	
	if params.has_key('border'):
		try: params['border'] = N.array(params['border'], dtype=N.int)[[0,1]]
		except: log.prNot(log.ERR, "border invalid, should be <int>,<int>.")
	
	if params.has_key('subap'):
		try: params['subap'] = N.array(params['subap'], dtype=N.int)
		except: log.prNot(log.ERR, "subap invalid, should be list of ints.")
	
	if params.has_key('offsets'):	
		if params['offsets'] and not os.path.exists(params['offsets']):
			log.prNot(log.ERR, "offset file '%s' does not exist." % \
		 		(params['offsets']))
	
	if params.has_key('outformat') and (params['outformat'] not in _OUTFORMATS):
		log.prNot(log.ERR, "Unsupported output '%s'" % (params['outformat']))
	
	if params.has_key('nheights'):
		params['nheights'] = N.array(params['nheights'], dtype=N.int)
	
	if params.has_key('layerheights'):
		try: 
			heights = []
			for hran in params['layerheights']:
				heights.append(hran.split('-'))
				print heights[-1]
			params['layerheights'] = N.array(heights, dtype=N.float)
			print params['layerheights']
		except: 
			log.prNot(log.ERR, "layerheights invalid, should be <float>-<float>,<float>-<float>,....")
	
	if params.has_key('layercells'):
		try: params['layercells'] = \
			 N.array(params['layercells'], dtype=N.int)[[0,1]]
		except: 
			log.prNot(log.ERR, "layercells invalid, should be <int>,<int>.")
		
	if params.has_key('skipsa'):
		params['skipsa'] = N.array(params['skipsa'], dtype=N.int)
	
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
	#elif (tool == 'tomo'):
		# tomo needs ccdres, aptr, shifts
	
	# Done


### ==========================================================================
### Tool classes
### ==========================================================================

class Tool(object):
	"""Generic Tool class with common functions."""
	
	def __init__(self, files, params):
		self.files = files
		self.nfiles = len(files)
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
					#lf.saveOldFile(uri)
					if (os.path.getsize(val) < 1024*1024 and not os.path.exists(uri)):
						import shutil
						shutil.copy2(val, uri)
						#os.link(val, uri)
		# Init list of output files we save
		self.ofiles = {}
		self.ofiles['params'] = lf.saveData(self.mkuri('astooki-params'), \
		 	params, aspickle=True)
		
	
	
	def load(self, filename):
		log.prNot(log.INFO, "Loading '%s'" % (filename))
		if not os.path.isfile(filename):
			log.prNot(log.WARNING, "'%s' is not a regular file." % (filename))
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
		return pyfits.getdata(filename)
	
	
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
		
		# Normalize everything outside the mask
		data[self.mask == False] = N.max(data[self.mask])
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
	"""
	Make a subaperture mask, given a number of geometric input parameters.
	
	This tool calculates a set of coordinates where subapertures or subimages 
	are located and saves these to disk. Optionally, the pattern can be plotted
	as well. The units for the parameters are arbitrary, but should be 
	consistent.
	
	The difference between subapertures and subimages is the plane they are 
	located in. Subapertures are located in the aperture plane (i.e. on the 
	lenslet) while subimages are located on the CCD. For this tool this is 
	obviously irrelevant, but the meaning can be quite different. Unfortunately, 
	some of the nomenclature in astooki is ambigious about this.
	
	@param rad Telescope aperture radius
	@param shape Shape of the telescope aperture (circular or square)
	@param sasize Subapertures size
	@param pitch Pitch of the subaperture positions
	@param xoff x-offset of even and odd rows in units of sasize. Set to [0,0.5]
		to get a brick-pattern grid, useful for hexagonal lenslet arrays.
	@param scale Scale the grid by this factor
	@param disp Displace the grid by this vector
	@param file Base filename to save the pattern to
	"""
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
		(nsa, llpos, cpos, size) = \
			libsh.calcSubaptConf(self.rad, self.sasize, self.pitch, self.shape, \
			self.xoff, self.disp, self.scale)
		# Save to file
		tmp = os.path.splitext(self.file)
		libsh.saveSaSfConf(self.mkuri(tmp[0]+'-origin'+tmp[1]), nsa, [-1,-1], \
		 	size, llpos)
		libsh.saveSaSfConf(self.mkuri(tmp[0]+'-centroid'+tmp[1]), nsa, [-1,-1], \
		 	size, cpos)
		if (self.plot):
			import astooki.libplot as libplot
			plrange = [list(N.array([-self.rad, self.rad]) + self.disp)]*2
			libplot.showSaSfLayout(self.mkuri(self.file+'-plot.eps'), llpos, size, \
				plrange=plrange)
		# Done
		
	


class SubfieldConfTool(Tool):
	"""
	Make a subfield mask, given a number of geometric input parameters.
	
	This is more or less the same as SubaptConfTool(), except that this tool
	generates a subfield mask and takes slightly different parameters. The 
	subfield mask will be a set of (pixel) coordinates relative to the origin of 
	a subaperture that define crops of the subimages on the CCD. Since pixel in 
	a subimage corresponds to a different field of view and view angle, a grid 
	of masks for a subimage correspond to a grid of different fields of view. 
	This grid can be used to calculate the subfield shifts that can consequently 
	be used to analyze the seeing recorded in the data.
	
	The grid generated is relative to the subimage, so the size of such a 
	subimage is relevant: subfields cannot lie *outside* this subimage. Also, 
	when later measuring the shifts for each subfield within a subimage, one 
	will have to move the subfield around to measure the cross-correlation. 
	Therefore one should supply a border which will take care of this and keep a 
	guard range free at the edge of the subimage. The shift range used later 
	should always be less or equal to the border supplied here.
	
	@param sasize The subaperture size to be used
	@param sfsize The subfield size
	@param overlap The overlap in x and y direction. 1 for complete overlap, 0 
	  for no overlap. [0.5, 0.5] gives about 50%% overlap
	@param border The border or guard range to keep clear withins sasize
	"""
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
		else: self.file = 'astooki-subfieldconf.csv'
		
		self.run()
	
	
	def run(self):
		# Generate subfield positions
		effsize = self.sasize - 2*self.border
		pitch = self.sfsize * (1-self.overlap)
		nsf = N.floor((effsize-self.sfsize+pitch) / pitch)
		effpitch = (effsize-self.sfsize)/(nsf-1)
		effpitch[(nsf == 1)] = 0
		
		log.prNot(log.INFO, "Effective size: (%g,%g), effpitch: (%g,%g)" % \
		 	(tuple(effsize) + tuple(effpitch)))
		
		sfpos = self.border + \
			N.indices(nsf, dtype=N.float).reshape(2,-1).T * effpitch
		sfpos = N.floor(sfpos).astype(N.int)
		totnsf = N.product(nsf).astype(N.int)
		
		log.prNot(log.NOTICE, "Found %d x %d subfields." % tuple(nsf))
		log.prNot(log.NOTICE, "Size %d,%d" % tuple(self.sfsize))
		libsh.saveSaSfConf(self.mkuri(self.file), totnsf, [-1,-1], self.sfsize, \
		 	sfpos)
		if (self.plot):
			import astooki.libplot
			libplot.showSaSfLayout(self.mkuri(self.file+'-plot.eps'), sfpos, \
			 	self.sfsize, plrange=[[0, self.sasize[0]], [0, self.sasize[1]]])
		# Done
		
	


class SubaptUpdateTool(Tool):
	"""
	Update subaperture mask with an offset.
	
	Since we want to compare different subfields within each subimage, we need 
	to know the reference direction of each subimage. Because of static 
	aberrations (telescope defocus, instrument issues) we cannot assume that 
	pixel (x,y) in subimage N corresponds to the same field of view as the same 
	pixel in subimage M. Or the other way around: given a granule G on the sun, 
	we want to know at what pixel that granule is located in each of the 
	subimages.
	
	To get these 'static offsets', we take a large field of view in one 
	reference subimage (almost the complete subimage) and cross-correlate that 
	with all subimages. This will give N shift vectors for each frame, N being 
	the number of subapertures. To get better results, it is possible to use 
	multiple subapertures as reference, this should give the same data and gives 
	an indication of the noise or reliability of the shift measurement. The 
	image shifts measured for different reference subapertures will be stored 
	alongside eachother.
	
	Because there is atmospheric seeing which causes tip-tilt of the images, the 
	image shifts will vary strongly from frame to frame. The seeing is 
	statistical, however, which means that the average shift over all frames 
	should be zero. To get the static offsets and identify the relative field of 
	view of all subimages, we average the image shifts over all frames. This 
	gives a offset vector for each subimage which indicates where the subimage 
	is pointing at relative to the other subimages.
	
	If we correct the subimage mask calculated earlier with this list of 
	vectors, each subimage mask will be pointing at the same location on the 
	sun. Once we have established this, we can subdivide the subimages in 
	different subfields knowing exactly where each subfield points and thus 
	providing a reliable method to base subsequent analysis on.
	
	@param maskfile The subaperture mask to be updated
	@param offsets A list of offset vectors
	"""
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
		newpos = (pos + off + maxsh).astype(N.int)
		# Crop the subaperture size by twice the maximum offset
		newsize = (size - maxsh*2).astype(N.int)
		# Store
		libsh.saveSaSfConf(self.mkuri(self.file), nsa, [-1,-1], newsize, newpos)
		if (self.plot):
			import astooki.libplot as libplot
			plfile = self.mkuri(self.file + '-plot.eps')
			# TODO: fix this range
			plran = [[0, N.ceil(1600./512)*512]]*2
			libplot.showSaSfLayout(plfile, newpos, newsize, \
				plrange=plran)
		# Done
	


class SubaptOptTool(Tool):
	"""
	Optimize subaperture configurations on real data.
	
	Although a subaperture mask made with SubaptConfTool() will be fairly 
	accurate if it is supplied with the right parameters, it is still possible 
	that the mask does not match the subimages exactly. To circumvent this 
	problem, SubaptOptTool() can take the statically generated subimage mask and 
	a flatfield image and match the mask onto the flatfield.
	
	The method used here is as follows. Given the mask, at each centroid grid
	position take a vertial slice out of the flatfield that is ~30 pixels wide 
	and twice the subimage high. This slice should cover the whole flatfielded 
	subimage. This slice is then averaged over the width, giving a 1 pixel wide 
	profile vertically across the subimage flatfield. The first pixel to the top 
	and bottom from the center of the slice where the intensity is lower than X 
	times the maximum intensity of the slice is considered to be the edge of the 
	subimage. This will give the height of the subimage. The same routine is 
	also done for horizontally across the subimage to get the width.
	
	The drawback of this routine is that is needs the flatfield and real image 
	to match perfectly. Fortunately, this should and ususally is that case. 
	Furthermore, if there is a speck of dust on the flatfielded image, this can 
	fool this tool. One should therefore always check the output, for example by 
	plotting the grid or testing the grid by using it as a mask on real data.
	
	@param maskfile The subaperture grid to optimize
	@param saifac The intensity dropoff factor to use (X in the above info)
	@param rad The radius of the subaperture grid pattern (used for plotting 
		only)
	
	"""
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
		libsh.saveSaSfConf(self.mkuri(self.file), onsa, [-1,-1], osize, opos)
		if (self.plot):
			import astooki.libplot as libplot
			# TODO: fix this range
			libplot.showSaSfLayout(self.mkuri(self.file+'-plot.eps'), opos, osize, \
				plrange=[[0, 2*self.rad]]*2)
		# Done
	


class ShiftTool(Tool):
	"""
	Calculate shifts between different subfields and subapertures, possibly 
	do some post-processing.
	
	@param shrange Shiftrange to measure in pixels (actual shiftrange will be 
		from -shrange to +shrange)
	@param safile Subaperture CCD position file
	@param sffile Subfield CCD position file
	@param nref Number of reference subapertures to use
	"""
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
		if (len(self.files) < 1):
			log.prNot(log.ERR, "ShiftTool(): Cannot continue without files!")
		
		# Run analysis
		self.run()
	
	
	def run(self):
		### Phase 1: Measure shifts
		### -----------------------
		
		self.calcShifts()
		
		### Phase 2: Process shifts
		### -----------------------
		
		if (self.nsf == 1):
			log.prNot(log.NOTICE, "Starting post-processing.")
			self.postProcess()
		
		### End: store metadata
		metafile = lf.saveData(self.mkuri('astooki-meta-data'), \
			self.ofiles, aspickle=True)
	
	
	def calcShifts(self):
		# Check if shifts already exist, load it if it does
		try: 
			metafile = self.mkuri('astooki-meta-data.pickle')
			# Read metafile, check if entry 'shifts' exists
			meta = lf.loadData(metafile, aspickle=True)
			import pyfits
			hdr = pyfits.getheader(self.mkuri(meta['shifts']['fits']))
			if (hdr.get('NAXIS') == 5 and \
			 	hdr.get('NAXIS1') == 2 and \
			 	hdr.get('NAXIS2') == self.nsf and \
			 	hdr.get('NAXIS3') == self.nsa and \
			 	hdr.get('NAXIS4') == self.nref and \
			 	hdr.get('NAXIS5') == self.nfiles):
				log.prNot(log.NOTICE, "Found previously measured shifts, restoring.")
				data, meta = lf.restoreData(metafile)
				self.ofiles = meta
				self.shifts = data['shifts']
				return
			else:
				raise Exception
		except:
			log.prNot(log.NOTICE, "Could not restore shift data, computing.")
		
		import astooki.clibshifts as ls
		
		allshifts = []
		allfiles = []
		allrefs = []
		for f in self.files:
			base = os.path.basename(f)
			img = self.load(f)
			if (img is None): 
				log.prNot(log.INFO, "Skipping %s, could not read file." % (base))
				continue
			allfiles.append(base)
			dfimg = self.darkflat(img)
			
			log.prNot(log.NOTICE, "Measuring shifts for %s." % (base))
			refaps = []
			imgshifts = ls.calcShifts(dfimg, self.saccdpos, self.saccdsize, \
			 	self.sfccdpos, self.sfccdsize, method=ls.COMPARE_ABSDIFFSQ, \
			 	extremum=ls.EXTREMUM_2D9PTSQ, refmode=ls.REF_BESTRMS, \
			 	refopt=self.nref, shrange=[self.shrange, self.shrange], \
			 	subfields=None, corrmaps=None, refaps=refaps)
			allrefs.append(refaps)
			allshifts.append(imgshifts)
		
		log.prNot(log.NOTICE, "Done, saving results to disk @ '%s'." % (self.file))
		self.shifts = N.array(allshifts)
		
		# Store the list of files where we save data to
		self.ofiles['shifts'] = lf.saveData(self.mkuri('image-shifts'), \
		 	self.shifts, asnpy=True, asfits=True)
		self.ofiles['refaps'] = lf.saveData(self.mkuri('referenace-subaps'), \
			allrefs, asnpy=True, ascsv=True)
		self.ofiles['saccdpos'] = lf.saveData(self.mkuri('subap-ccdpos'), \
		 	self.saccdpos, asnpy=True, asfits=True)
		self.ofiles['sfccdpos'] = lf.saveData(self.mkuri('subfield-ccdpos'), \
		 	self.sfccdpos, asnpy=True, asfits=True)
		self.ofiles['saccdsize'] = lf.saveData(self.mkuri('subap-ccdsize'), \
		 self.saccdsize, asnpy=True, asfits=True)
		self.ofiles['sfccdsize'] = lf.saveData(self.mkuri('subfield-ccdsize'), \
		 	self.sfccdsize, asnpy=True, asfits=True)
		cpos = self.saccdpos.reshape(-1,1,2) + self.sfccdpos.reshape(1,-1,2) + \
			self.sfccdsize.reshape(1,1,2)/2.0
		self.ofiles['sasfpos-c'] = lf.saveData(self.mkuri('sasfpos-c'), cpos, \
		 	asnpy=True, asfits=True)
		self.ofiles['files'] = lf.saveData(self.mkuri('processed-files'), \
		 	allfiles, asnpy=True, ascsv=True, csvfmt='%s')
		# Add meta info
		self.ofiles['path'] = os.path.dirname(os.path.realpath(self.mkuri('tst')))
	
	
	def postProcess(self):
		# Process NaNs and other non-finite numbers
		notfin = N.argwhere(N.isfinite(self.shifts) == False)
		if (notfin.shape[0] > 0):
			log.prNot(log.WARNING, "Found non-finite shifts! Check configuration.")
			notfin_perc = notfin.shape[0]*100./self.shifts.size
			log.prNot(log.WARNING, "%d (%.2g%%) non-finite entries, spread over %d frames, %d subaps, %d subfields." % \
		 	(notfin.shape[0], notfin_perc, N.unique(notfin[:,0]).size, \
		 	N.unique(notfin[:,2]).size, N.unique(notfin[:,3]).size))
		
		# Find the shift variance per subaperture
		log.prNot(log.NOTICE, "Calculating shift statistics.")
		
		# First average over all Nref reference subapertures
		s_ref = N.mean(self.shifts[:,:,:,0,:], axis=1)
		# Now make sure the average *per frame* is zero
		s_avgfr = N.mean(s_ref, axis=1)
		s_norm = s_ref - s_avgfr.reshape(-1,1,2)
		
		# Calculate variance per subaperture
		savar = N.var(s_norm, axis=0)
		log.prNot(log.NOTICE, "Average shift variance: (%g,%g), max: (%g,%g)" % \
			(tuple(N.mean(savar,0)) +  tuple(N.max(savar,0))))
		self.ofiles['shift-var'] = lf.saveData(self.mkuri('shift-variance'), \
			savar, asnpy=True, ascsv=True)
		
		# Now average over all frames to get the offset. Also calculate the error
		log.prNot(log.NOTICE, "Calculating static offsets.")
		soff = N.mean(s_norm, axis=0)
		sofferr = (N.var(s_norm, axis=0))**0.5
		
		self.ofiles['offsets'] = lf.saveData(self.mkuri('static-offsets'), \
			soff, asnpy=True, ascsv=True)
		self.ofiles['offset-err'] = lf.saveData(\
			self.mkuri('static-offset-err'), sofferr, asnpy=True, ascsv=True)
		if (self.plot):
			import astooki.libplot as libplot
			libplot.plotShifts(self.mkuri('static-offset-plot'), self.shifts, \
				self.saccdpos, self.saccdsize, self.sfccdpos, self.sfccdsize, \
				plorigin=(0,0), plrange=(2048, 2048), mag=7.0, allsh=False, \
			 	title='Static offsets, mag=7', legend=True)
			libplot.plotShifts(self.mkuri('static-offset-plot-250'), \
			 	self.shifts[:250], self.saccdpos, self.saccdsize, self.sfccdpos, \
			 	self.sfccdsize, plorigin=(0,0), plrange=(2048, 2048), mag=7.0, \
			 	allsh=False, title='Static offsets, mag=7', legend=True)
	


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
		self.safile = params['safile']
		self.sffile = params['sffile']
		self.shfile = params['shifts']
		# Load safile and sffile
		(self.nsa, self.saccdpos, self.saccdsize) = \
			libsh.loadSaSfConf(self.safile)
		(self.nsf, self.sfccdpos, self.sfccdsize) = \
			libsh.loadSaSfConf(self.sffile)
		# Load shifts
		self.shifts = lf.loadData(self.shfile, asnpy=True)
		self.run()
	
	
	def run(self):
		# Process NaNs and other non-finite numbers
		notfin = N.argwhere(N.isfinite(self.shifts) == False)
		if (notfin.shape[0] > 0):
			log.prNot(log.WARNING, "Found non-finite shifts! Check configuration.")
			notfin_perc = notfin.shape[0]*100./self.shifts.size
			log.prNot(log.WARNING, "%d (%.2g%%) non-finite entries, spread over %d frames, %d subaps, %d subfields." % \
		 	(notfin.shape[0], notfin_perc, N.unique(notfin[:,0]).size, \
		 	N.unique(notfin[:,2]).size, N.unique(notfin[:,3]).size))
		
		# If we have one subfield, treat it as static shift data:
		if (len(self.sfccdpos) == 1):
			log.prNot(log.NOTICE, "Calculating static offsets.")
			(soff, sofferr) = libsh.procStatShift(self.shifts[:,:,:,0,:])
			self.ofiles['offsets'] = lf.saveData(self.mkuri('static-offsets'), \
				soff, asnpy=True, ascsv=True)
			self.ofiles['offset-err'] = lf.saveData(\
				self.mkuri('static-offset-err'), sofferr, asnpy=True, ascsv=True)
			if (self.plot):
				import astooki.libplot as libplot
				libplot.plotShifts(self.mkuri('static-offset-plot'), self.shifts, \
					self.saccdpos, self.saccdsize, self.sfccdpos, self.sfccdsize, \
					plorigin=(0,0), plrange=(2048, 2048), mag=7.0, allsh=False, \
				 	title='Static offsets, mag=7', legend=True)
		# We have subfield data here, process
		else:
			log.prNot(log.WARNING, "Cannot process this data, not in the proper format.")
		
		# Store the new meta file
		metafile = lf.saveData(self.mkuri('astooki-meta-data'), \
			self.ofiles, aspickle=True)
	


class TomoTool(Tool):
	"""Tomographically analyze WFWFS data"""
	def __init__(self, files, params):		
		super(TomoTool, self).__init__(files, params)
		# Load shift data
		self.shifts = lf.loadData(params['shifts'], asnpy=True)
		# Load subaperture centroid positions
		(self.nsa, self.sapos, self.sasize) = \
		 	libsh.loadSaSfConf(params['safile'])
		# Load subfield pixel positions
		(self.nsf, self.sfccdpos, self.sfsize) = \
			libsh.loadSaSfConf(params['sffile'])
		# Store ccd resolution in radians (arcsec -> radian == /60/60 * pi/180)
		self.ccdres = params['ccdres'] * N.pi /60./60./180.
		# Telescope aperture radius
		self.aptr = params['aptr']
		# Geometry (layer heights and number of cells)
		self.nlay = params['nheights']
		lh = params['layerheights']
		self.lcells = params['layercells']
		self.geoms = N.zeros((N.product(self.nlay), len(self.nlay)))
		self.origs = N.zeros((N.product(self.nlay), len(self.nlay)))
		self.sizes = N.zeros((N.product(self.nlay), len(self.nlay)))
		# Calculate effective (subfield) FoV 
		self.fov = (N.max(self.sfccdpos + self.sfsize, axis=0) - \
			N.min(self.sfccdpos, axis=0)) * self.ccdres 
		self.sffov = self.sfsize * self.ccdres
		# Calculate subfield pointing angles, set average to 0
		self.sfang = (self.sfccdpos + self.sfsize/2.0) * self.ccdres
		self.sfang -= N.mean(self.sfang, axis=0)
		# Rotate the coordinate system
		self.sfang *= N.array([-1.0,-1.0])
		
		# Setup all layer configurations
		for l in range(len(self.nlay)):
			thisHeights = N.linspace(lh[l,0], lh[l,1], self.nlay[l])
			# All height permutations such that we can loop over them
			self.geoms[:, l] = N.tile(\
				N.repeat(thisHeights, N.product(self.nlay[l+1:])), \
				 	N.product(self.nlay[:l]))
		
		# Layer origin and sizes
		telang = [0.0, 0.0]
		self.lorigs = self.geoms.reshape(N.product(self.nlay), \
			 	len(self.nlay), 1) * N.tan(telang).reshape(1,1,2)
		self.lsizes = self.aptr + \
				self.geoms.reshape(N.product(self.nlay), len(self.nlay), 1)*\
				N.tan(0.5 * self.fov).reshape(1,1,2)
		
		# Allocate data for reconstruction
		self.recatm = N.zeros((self.shifts.shape[0], N.product(self.nlay), \
		 	len(self.nlay), self.lcells[0], self.lcells[1], 2), dtype=N.float32)
		self.inrms = N.zeros((self.shifts.shape[0], 4))
		self.recrms = N.zeros((self.shifts.shape[0], N.product(self.nlay), 4))
		self.diffrms = N.zeros((self.shifts.shape[0], N.product(self.nlay), 4))
		
		log.prNot(log.NOTICE, "Starting tomographical analysis of WFWFS data stored in '%s' using %d layers with each %dx%d cells." % (params['shifts'], len(self.nlay), self.lcells[0], self.lcells[1]))
		
		self.run()
	
	
	def run(self):
		import astooki.libtomo as lt
		# Setup SVD cache for inverting data
		svdCache = lt.cacheSvd(self.geoms, self.lsizes, self.lorigs, \
			self.lcells, self.sasize, self.sapos, self.sfang, self.sffov, \
			matroot='workdocs/data/matrices/')
		# Check if values are finite
		notfin = N.argwhere(N.isfinite(self.shifts) == False)
		if (notfin.shape[0] > 0):
			log.prNot(log.WARNING, "Some measurements are non-finite, check configuration!")			
			# TODO: setting to zero is a poor solution to NaNs
			for nfidx in notfin:
				self.shifts[tuple(nfidx)] = 0.0
		
		# Setup inversion and forward matrices from SVD components
		modmats = {}
		modmats['inv'] = []
		modmats['fwd'] = []
		log.prNot(log.NOTICE, "Pre-computing inversion and forward matrices...")
		for geom in range(N.product(self.nlay)):
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
		shifts = self.shifts.mean(axis=1)
		# Remove average over all frames
		shifts -= N.mean(shifts, axis=0).reshape(1, shifts.shape[1], shifts.shape[2], 2)
		for it in range(shifts.shape[0]):
			log.prNot(log.NOTICE, "Data frame %d/%d" % (it+1, shifts.shape[0]))
			sh = shifts[it]
			# RMS of input shifts:
			self.inrms[it] = N.r_[\
				self.rms(sh[...,0]), \
				self.rms(sh[...,1]), \
				self.rms(sh[...,0], remdc=True), \
				self.rms(sh[...,1], remdc=True)]
			#log.prNot(log.NOTICE, "Inrms: %g, %g." % tuple(self.inrms[it]))
			# Vectorize shift measurement
			wfwfsX = sh[...,0].flatten()
			wfwfsY = sh[...,1].flatten()
			for geom in range(N.product(self.nlay)):
				#log.prNot(log.NOTICE, "Starting reconstruction for geometry %d/%d." % (geom+1, N.product(self.nlay)))
				# Invert measurements to atmosphere
				# Reconstruct atmospheric x-slopes in current geometry:
				self.recatm[it, geom, :, :, :, 0] = \
				 	(N.dot(modmats['inv'][geom], wfwfsX)).reshape( \
						len(self.nlay),
						self.lcells[0],
						self.lcells[1])
					
				# Reconstruct atmospheric y-slopes in current geometry:
				self.recatm[it, geom, :, :, :, 1] = \
				 	(N.dot(modmats['inv'][geom], wfwfsY)).reshape( \
						len(self.nlay),
						self.lcells[0],
						self.lcells[1])
				
				# Forward calculate all reconstruction to compare model with the 
				# measurements
				#print modmats['inv'][geom].shape, wfwfsX.shape, modmats['fwd'][geom].shape, self.recatm[it, geom, :, :, :, 0].shape
				
				wfwfsXRec = N.dot(modmats['fwd'][geom], \
					self.recatm[it, geom, :, :, :, 0].flatten())
				wfwfsYRec = N.dot(modmats['fwd'][geom], \
					self.recatm[it, geom, :, :, :, 1].flatten())
				
				# Calculate reconstruction RMS
				self.recrms[it, geom] = N.r_[\
					self.rms(wfwfsXRec), \
					self.rms(wfwfsYRec), \
					self.rms(wfwfsXRec, remdc=True), \
					self.rms(wfwfsYRec, remdc=True)]
				
				# Calculate difference RMS
				self.diffrms[it, geom] = N.r_[\
						self.rms(wfwfsXRec - wfwfsX), \
						self.rms(wfwfsYRec - wfwfsY), \
						self.rms(wfwfsXRec - wfwfsX, remdc=True), \
						self.rms(wfwfsYRec - wfwfsY, remdc=True)]
				#log.prNot(log.NOTICE, "Recrms: %g, %g, diffrms: %g, %g." % (tuple(self.recrms[it, geom]) + tuple(self.diffrms[it, geom])))
		
		log.prNot(log.NOTICE, "Done, saving data.")
		self.ofiles['inrms'] = lf.saveData(self.mkuri('tomotool-inrms'), \
		 	self.inrms, asnpy=True, asfits=True)
		self.ofiles['recrms'] = lf.saveData(self.mkuri('tomotool-recrms'), \
		 	self.recrms, asnpy=True, asfits=True)
		self.ofiles['diffrms'] = lf.saveData(self.mkuri('tomotool-diffrms'), \
		 	self.diffrms, asnpy=True, asfits=True)
		self.ofiles['geoms'] = lf.saveData(self.mkuri('tomotool-geoms'), \
		 	self.geoms, asnpy=True, asfits=True)
		# Save metafile
		metafile = lf.saveData(self.mkuri('astooki-meta-data'), \
			self.ofiles, aspickle=True)

	
	
	def rms(self, data, remdc=False):
		if remdc:
			return N.sqrt(N.mean((data-N.mean(data))**2.0))
		else:
			return N.sqrt(N.mean(data**2.0))

	


class SdimmTool(Tool):
	"""
	Perform SDIMM+ analysis on shift data.
	
	- Loop over all rows of subapertures (subapertures that are at the same y 
	  position)
	- For each subap row choose a reference subaperture (i.e. the left-most one)
	- Loop over all rows of subfields (subfields with the same y coordinate)
	- For each subfield row, choose a reference subfield
	- Loop over all other subapertures in this subap row
	- Loop over all other subfields in this subfield row
	- Compare all subaperture-subfield pairs as described in Scharmer & van 
	  Werkhoven
	- Repeat this for all columns
	
	This tool outputs the raw results of the calculations to sdimmraw.<fits|npy> 
	and is a is an N * 9 matrix where each row holds the following information:
    [id, s, a, C_lsa, C_tsa, refsa, sa, refsf, sf]
  with:
	- id=0 for row-wise comparison and 1 for column-wise (as described 
    above),
  - s the scalar distance between the two subapertures in meters,
  - a the scalar angle between the two subfields in pixels (convert with the 
    CCD scale to get a real angle), 
  - C_lsa the longitudinal covariance between the two sequences of 
    differential image shifts (as described in the paper)
  - C_tsa the transversal covariance
  - refsa the index of the reference subaperture used here
  - refsf the ubdex of the reference subfield used here
  - sa the index of the other subaperture used
  - sf the index of the other subfield used
  
  Besides the raw information, this tool also outputs processed information to
  sdimmcol* and sdimmrow* files, where the *col* files have information on the 
  column-wise comparison and the *row* files on the row-wise comparison of the
  data.
	
	sdimm<col|row>.* is a 3 x N x M matrix with N the number of unique
  subaperture distances (s) and M the number of unique angles (a) for this 
  data. For the 85 subaperture lenslet array at the SST, N is 5 for
  column-wise comparison and 9 for row-wise comparison. This is the data that
  should be decomposde in SDIMM basis described in the paper. 
 	
  The first N x M frame holds the longitudinal covariance, the second frame
  holds the transversal covariance and the third frame holds the number of
  covariances each specific cell was averaged over (i.e. given an (s,a)
  coordinate, how many covariances were calculated?
	
  sdimm<col|row>-s.* hold the unique subaperture distances mentioned above, 
  and sdimm<col|row>-a.* hold the unique subfield angles mentioned above.
	
	The following parameters are required as input for this tool:
	@param safile centroid subaperture positions [meter]
	@param sffile centroid subfield positions [pixel]
	@param ccdres angular ccd resolution [arcsec/pix]
	@param aptr aperture radius [meter)]
	"""
	def __init__(self, files, params):
		super(SdimmTool, self).__init__(files, params)
		# Load shift data
		self.shifts = lf.loadData(params['shifts'], asnpy=True)
		# Load subaperture centroid positions
		(self.nsa, self.sapos, self.sasize) = \
		 	libsh.loadSaSfConf(params['safile'])
		# Load subfield pixel positions
		(self.nsf, self.sfccdpos, self.sfsize) = \
			libsh.loadSaSfConf(params['sffile'])
		# Skip these subaps in the analysis
		self.skipsa = N.array(params['skipsa'], dtype=N.int)
		# Store ccd resolution in radians (arcsec -> radian == /60/60 * pi/180)
		#self.ccdres = params['ccdres'] * N.pi /60./60./180.
		# Telescope aperture radius
		#self.aptr = params['aptr']
		# Calculate effective (subfield) FoV 
		# self.fov = (N.max(self.sfccdpos + self.sfsize, axis=0) - \
		# 	N.min(self.sfccdpos, axis=0)) * self.ccdres 
		# self.sffov = self.sfsize * self.ccdres
		# Calculate subfield pointing angles, set average to 0
		# self.sfang = (self.sfccdpos + self.sfsize/2.0) * self.ccdres
		# self.sfang -= N.mean(self.sfang, axis=0)
		# Rotate the coordinate system
		#self.sfang *= N.array([-1.0,-1.0])
		
		log.prNot(log.NOTICE, "Starting SDIMM+ analysis of WFWFS data stored in '%s'." % (params['shifts']))
		
		self.run()
	
	
	def run(self):
		# Check for non-finite values (shouldn't be, just to make sure)
		# notfin = N.argwhere(N.isfinite(self.shifts) == False)
		# if (notfin.shape[0] > 0):
		# 	log.prNot(log.WARNING, 
		# 		"Some measurements are non-finite, check configuration!")
		# 	# TODO: setting to zero is a poor solution to fixing NaNs
		# 	for nfidx in notfin:
		# 		self.shifts[tuple(nfidx)] = 0.0
		
		# Average over number of references
		shifts = self.shifts.mean(axis=1)
		
		# This will hold the sdimm correlation values
		sdimm = []
		
		### Loop over all *rows*
		### ====================
		# Get unique SA row positions
		sarows = N.unique(self.sapos[:,1])
		# Get unique SF row positions
		sfrows = N.unique(self.sfccdpos[:,1])		
		# Loop over all subaperture rows
		for sarowpos in sarows:
			# Get a list of all subapertures at this row (i.e. same y coordinate)
			salist = N.argwhere(self.sapos[:,1] == sarowpos).flatten()
			# Exclude bad subaps
			salist = N.lib.arraysetops.setdiff1d(salist, self.skipsa)
			# Take a reference subaperture in this row (the one on the left)
			refsa = salist[N.argmin(self.sapos[salist][:,0])]
			# Loop over all subapertures in this row
			for rowsa in salist:
				#if (rowsa == refsa): continue
				log.prNot(log.NOTICE, "ROW: Comparing subap %d with subap %d." % \
					(refsa, rowsa))
				# Calculate the distance between these two subaps
				s = self.sapos[rowsa, 0] - self.sapos[refsa, 0]
				# Loop over all subfield rows
				for sfrowpos in sfrows:
					# Get a list of all subfields at this row (i.e. same y coordinate)
					sflist = N.argwhere(self.sfccdpos[:,1] == sfrowpos).flatten()
					# Take a reference subaperture in this row (the one on the left)
					refsf = sflist[N.argmin(self.sfccdpos[sflist][:,0])]
					# Loop over all subfields in this row
					for rowsf in sflist:
						#if (rowsf == refsf): continue
						# Calculate the angle between these subfields (in *pixels*! 
						# multiply with pixel scale to get real angles)
						a = self.sfccdpos[rowsf, 0] - self.sfccdpos[refsf, 0]
						# Compare subap <refsa> with <rowsa> here and subfield <refsf>
						# with <rowsf>. Follow the notation in Scharmer & van Werkhoven:
						# differential image shifts:
						dx_s0 = shifts[:, refsa, refsf, :] - shifts[:, rowsa, refsf, :]
						dx_sa = shifts[:, refsa, rowsf, :] - shifts[:, rowsa, rowsf, :]
						# Unscaled longitudinal and transversal covariance of these shifts
						C_lsa = (N.cov(dx_s0[:,0], dx_sa[:,0]))[0,1]
						C_tsa = (N.cov(dx_s0[:,1], dx_sa[:,1]))[0,1]
						# Add all values to the matrix
						sdimm.append([0, s, a, C_lsa, C_tsa, refsa, rowsa, refsf, rowsf])
		
		### Loop over all *columns*
		### =======================
		sacols = N.unique(self.sapos[:,0])
		sfcols = N.unique(self.sfccdpos[:,0])		
		for sacolpos in sacols:
			salist = N.argwhere(self.sapos[:,0] == sacolpos).flatten()
			refsa = salist[N.argmin(self.sapos[salist][:,1])]
			for colsa in salist:
				#if (colsa == refsa): continue
				log.prNot(log.NOTICE, "COLUMN: Comparing subap %d with subap %d." % \
					(refsa, colsa))
				s = self.sapos[colsa, 1] - self.sapos[refsa, 1]
				for sfcolpos in sfcols:
					sflist = N.argwhere(self.sfccdpos[:,0] == sfcolpos).flatten()
					refsf = sflist[N.argmin(self.sfccdpos[sflist][:,1])]
					for colsf in sflist:
						#if (colsf == refsf): continue
						a = self.sfccdpos[colsf, 1] - self.sfccdpos[refsf, 1]
						dx_s0 = shifts[:, refsa, refsf, :] - shifts[:, colsa, refsf, :]
						dx_sa = shifts[:, refsa, colsf, :] - shifts[:, colsa, colsf, :]
						C_lsa = (N.cov(dx_s0[:,1], dx_sa[:,1]))[0,1]
						C_tsa = (N.cov(dx_s0[:,0], dx_sa[:,0]))[0,1]
						sdimm.append([1, s, a, C_lsa, C_tsa, refsa, colsa, refsf, colsf])
		
		# Convert to numpy array and store to disk
		sdimm = N.array(sdimm)
		self.ofiles['sdimmraw'] = lf.saveData(self.mkuri('sdimmraw'), \
		 	sdimm, asnpy=True, asfits=True)
		
		# Now filter sdimm values so that we have a two matrices of correlations,
		# one for transversal (C_tsa) and one for longitudinal (C_lsa) 
		# correlations
		
		# FIXME: Need to round off 's' values because we get numerical errors
		sdimm[:,1] = N.round(sdimm[:,1], 7)
		sdrow = sdimm[N.argwhere(sdimm[:,0] == 0).flatten()]
		sdcol = sdimm[N.argwhere(sdimm[:,0] == 1).flatten()]
		
		### Process row-wise data here
		### ==========================
		# Unique s and a values:
		uns = N.unique(sdrow[:,1])
		una = N.unique(sdrow[:,2])
		
		# Fill matrix SDimmRowCorrelation
		sdrc = N.zeros(((3,) + uns.shape + una.shape))
		for ns in xrange(len(uns)):
			s = uns[ns]
			for na in xrange(len(una)):
				a = una[na]
				idx = N.argwhere((sdrow[:,1] == s) & (sdrow[:,2] == a))
				if len(idx) == 0:
					log.prNot(log.WARNING, "Warning: found 0 at (%g, %g)!" % (s, a))
				else: 
					sdrc[0, ns, na] = N.mean(sdrow[idx, 3])
					sdrc[1, ns, na] = N.mean(sdrow[idx, 4])
					sdrc[2, ns, na] = len(idx)
		# Save correlation values to disk
		self.ofiles['sdimmrow'] = lf.saveData(self.mkuri('sdimmrow'), \
		 	sdrc, asnpy=True, asfits=True)
		# Save s and a values to disk
		self.ofiles['sdimmrow-s'] = lf.saveData(self.mkuri('sdimmrow-s'), \
			uns, asnpy=True, asfits=True, ascsv=True)
		self.ofiles['sdimmrow-a'] = lf.saveData(self.mkuri('sdimmrow-a'), \
			una, asnpy=True, asfits=True, ascsv=True)
		
		### Process column-wise data here
		### =============================
		# Unique s and a values:
		uns = N.unique(sdcol[:,1])
		una = N.unique(sdcol[:,2])
		
		# Fill SDimmColumnCorrelation
		sdcc = N.zeros(((3,) + uns.shape + una.shape))
		for ns in xrange(len(uns)):
			s = uns[ns]
			for na in xrange(len(una)):
				a = una[na]
				idx = N.argwhere((sdcol[:,1] == s) & (sdcol[:,2] == a))
				if len(idx) == 0: 
					log.prNot(log.WARNING, "Warning: found 0 at (%g, %g)!" % (s, a))
				else: 
					sdcc[0, ns, na] = N.mean(sdcol[idx, 3])
					sdcc[1, ns, na] = N.mean(sdcol[idx, 4])
					sdcc[2, ns, na] = len(idx)
		# Save to disk
		self.ofiles['sdimmcol'] = lf.saveData(self.mkuri('sdimmcol'), \
		 	sdcc, asnpy=True, asfits=True)
		self.ofiles['sdimmcol-s'] = lf.saveData(self.mkuri('sdimmcol-s'), \
			uns, asnpy=True, asfits=True, ascsv=True)
		self.ofiles['sdimmcol-a'] = lf.saveData(self.mkuri('sdimmcol-a'), \
			una, asnpy=True, asfits=True, ascsv=True)
		metafile = lf.saveData(self.mkuri('sdimm-meta-data'), \
			self.ofiles, aspickle=True)
	


class SimulShift(Tool):
	"""Simulate WFWFS measurements using N discrete KL phase screens"""
	def __init__(self, files, params):		
		super(TomoTool, self).__init__(files, params)
		# need: nlayer, ncells, lheights, lstrength, lorigs, lsizes, sapos,
		# sasize, sfpos, sfsize
		
		# Load subaperture centroid positions
		(self.nsa, self.sapos, self.sasize) = \
		 	libsh.loadSaSfConf(params['safile'])
		# Load subfield pixel positions
		(self.nsf, self.sfccdpos, self.sfsize) = \
			libsh.loadSaSfConf(params['sffile'])
		self.aptr = params['aptr']
		
		# Geometry (layer heights and number of cells)
		self.nlay = params['nheights']
		self.lh = params['layerheights']
		self.lcells = params['layercells']
		self.lorigs = N.zeros((N.product(self.nlay), len(self.nlay)))
		self.lsizes = N.zeros((N.product(self.nlay), len(self.nlay)))
		
		# Calculate effective (subfield) FoV 
		self.fov = (N.max(self.sfccdpos + self.sfsize, axis=0) - \
			N.min(self.sfccdpos, axis=0)) * self.ccdres 
		self.sffov = self.sfsize * self.ccdres
		# Calculate subfield pointing angles, set average to 0
		self.sfang = (self.sfccdpos + self.sfsize/2.0) * self.ccdres
		self.sfang -= N.mean(self.sfang, axis=0)
		
		# Layer origin and sizes
		telang = [0.0, 0.0]
		self.lorigs = self.reshape(-1,1) * N.tan(telang).reshape(1,2)
		self.lsizes = self.aptr + \
				self.reshape(-1,1) * N.tan(0.5 * self.fov).reshape(1,1,2)
		
		log.prNot(log.NOTICE, "Starting seeings simulation using %d layers with each %dx%d cells." % (len(self.nlay), self.lcells[0], self.lcells[1]))
		
		self.run()
	
	
	def run(self):
		pass
	
	
	def genCn2(self, maxh=25000, mode=0, n=5000):
		"""
		Generate a C_n^2 profile up to height 'maxh' in 'n' steps. 'mode' 
		determines the type of C_n^2 profile to generate, currently only the H-V 
		5/7 model is supported.
		"""
		# Generate height array and empty C_n^2 profile
		height = linspace(0.0, maxh, n)
		cn2 = zeros(len(height))
		# Generate C_n^2 values
		if mode == 0:
			# Generate H-V 5/7 profile (Tyson p10)
			W = 21.
			A = 1.7e-14
			height /= 1000.0
			cn2 = 5.94e-23 * height**10 * (W/27)**2 * exp(-height) + \
				2.7e-16 * exp(-2 * height / 3) + A * exp(-10 * height)
			height *= 1000.
		else:
			raise ValueError('Unknown mode (should be 0)')
		
		# Return the C_n^2 profile with associated heights
		return array([height, cn2])
	
	
	def readRadialKl(self, filename, nModes=500):
		"""
		Read in the radial KL profiles from disk, limiting the number of modes
		read in to 'nModes'.
		"""
		fd = open(filename, 'r')
		
		n_e = int(fd.readline()) # Number of 'raw' KL modes (unique q's)
		n_r = int(fd.readline()) # Number of radial points per profile
		# Estimate the number of modes present if not set
		n_m = nModes if (nModes) else 2*n_e
		
		# Skip four lines
		for i in range(4):
			fd.next()
		
		# Allocate memory for the various KL quantities
		kl_p = N.zeros(n_m, N.int)
		kl_q = N.zeros(n_m, N.int)
		kl_e = N.zeros(n_m, N.float)
		# Memory for the radial coordinates
		kl_r = N.zeros(n_r, N.float)
		# Memory for the radial KL modes
		kl = N.zeros((n_r, n_m), N.float)
		
		# Read in radial coordinates
		for i in range(n_r):
			kl_r[i] = N.float(fd.next())
		
		# Read in KL modes. File starts with mode 2 (tip), piston is not included.
		jj = 1
		nold = 0
		log.prNot(\
			log.NOTICE, "Starting reading in %d KL modes at %ld" % (n_e, fd.tell()))
		
		for i in xrange(n_e):
			# Use the line listing the KL mode number as consistency check
			n = N.int(fd.next())
			if (n != nold+1):
				raise IOError(-1, ("Reading in KL modes from " + filename + \
					" failed, KL modes do not increment correctly."))
			nold = n
			
			# Read some KL mode properties
			kl_e[jj] = float(fd.next())
			kl_p[jj] = int(fd.next())
			kl_q[jj] = int(fd.next())
			
			# Read in the radial values of the KL mode
			for r in xrange(n_r):
				kl[r,jj] = float(fd.next())
			
			# Increase the number of KL modes read, stop if we have enough
			jj += 1
			if jj >= n_m: break
			
			# If kl_q is not zero, we can re-use this base-mode
			if kl_q[jj-1] != 0:
				kl_e[jj] = kl_e[jj-1]
				kl_p[jj] = kl_p[jj-1]
				kl_q[jj] = -kl_q[jj-1]
				kl[:,jj] = kl[:,jj-1]
				jj += 1
				if jj >= n_m: break
		
		prNot(log.NOTICE, \
			"Read %d modes and %g kilobytes" % (jj, fd.tell()/1000.))
		fd.close()
		
		# Sanity checking here
		if not (N.alltrue(isfinite(kl)) and N.alltrue(kl_e >= 0)):
			raise ValueError("Error reading KL modes, some values are not finite")
		
		return {
			'r': kl_r, 
			'e': kl_e, 
			'p': kl_p, 
			'q': kl_q, 
			'kl': kl
		}
	
	
	

### ==========================================================================
### Helper functions
### ==========================================================================

def print_help(tool, out=sys.stdout):
	if (out == sys.stdout):
		print >> out, help_message['preamble']
		print >> out, help_message['common']
		if (tool != 'common'):
			print >> out, help_message[tool]
		sys.exit(0)
	else:
		print >> out, help_message[tool]
		sys.exit(-1)


if __name__ == "__main__":
	sys.exit(main())
