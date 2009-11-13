#!/usr/bin/env python2.5
# encoding: utf-8

##
# @file pyatk.py
# @brief AsTooki: the astronomical toolkit - imagemagick for astronomers.
# @author Tim van Werkhoven (tim@astro.su.se)
# @date 20090422 14:34
# 
# Created by Tim van Werkhoven on 2009-04-22.
# Copyright (c) 2009 Tim van Werkhoven (tim@astro.su.se)
# 
# This file is licensed under the Creative Commons Attribution-Share Alike
# license versions 3.0 or higher, see
# http://creativecommons.org/licenses/by-sa/3.0/

import sys, os, time
import astooki.libfile as lf
import astooki.liblog as log
import astooki.libsh as libsh
import getopt
import numpy as N
import scipy as S

GITREVISION="v20090626.0-22-gd3205c2"
VERSION = "0.1.0-%s" % (GITREVISION)
AUTHOR = "Tim van Werkhoven (tim@astro.su.se)"
DATE = "20090623"

##
# @brief Run-time command line help for the various tools
help_message = {}
help_message['preamble'] = """astooki version %s (%s) by %s.
Usage: pyatk.py <TOOL> [OPTIONS] [FILES]""" % (VERSION, DATE, AUTHOR)

help_message['syntaxerr'] = """pyatk.py: Syntax incorrect.
Usage: 
   pyatk.py <TOOL> [OPTIONS] [FILES]
or
   pyatk.py [TOOL] --help
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
 -n, --nref=INT              number of references to use [4]
     --comp=METHOD           img shift algorithm to use, 'sqd' for square
                               difference, 'adsq' for abs. difference squared
     --intpl=METHOD          subpixel interpolation method to use, '9pt' for
                               9-point quadratic, '5pt' for 5-point quadratic.
     --mask=MASK             mask to use when correlating (circular or
                               none)"""

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
     --shifts-range=BEG1,END1,BEGN,ENDN  
                             do not use all shift measurements, but only those 
                               from BEG to END
     --shifts-n=N            split up the shifts in parts of N frames, and 
                               process each part equally
     --safile=FILE           centroid positions for each subaperture IN 
                               APERTURE SPACE (meter)
     --sffile=FILE           subfield positions on the CCD (pixels)
 -n, --nref=INT              number of references to use [0 = all]
     --skipsa=SA1,SA2,...    list of subapertures to skip in analysis"""

## @brief ANA data format
_FORMAT_ANA = 'ana'
## @brief FITS data format
_FORMAT_FITS = 'fits'
## @brief PNG image format
_FORMAT_PNG = 'png'
## @brief Binary NumPy data format
_FORMAT_NPY = 'npy'
## @brief Supported input formats
_INFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_NPY)
## @brief Supported output formats
_OUTFORMATS = (_FORMAT_ANA, _FORMAT_FITS, _FORMAT_PNG, _FORMAT_NPY)
## @brief Available tools (used for checking commands)
_TOOLS = ('convert', 'stats', 'shiftoverlay', 'samask', 'sfmask', 'saopt', 'saupd', 'shifts', 'tomo', 'sdimm')
## @brief Default floattype to use
_ftype = N.float64
## @brief Default inttype to use
_itype = N.int32

### ==========================================================================
### Startup functions
### ==========================================================================

## @brief Parse command line arguments and start the right tool class
def main(argv=None):
	print ">> This is astooki %s by %s" % (VERSION, AUTHOR)
	#print "Sleeping, attach debuggers now."
	#time.sleep(10)
	__begin = time.time()
	# Parse command-line options from sys.argv
	(tool, params, files) = parse_options()
	# Sanity check on parameters, init log file
	check_params(tool, params)
	log.prNot(log.NOTICE, ">> This is astooki %s by %s" % (VERSION, AUTHOR))
	cmd = sys.argv[0]
	for arg in sys.argv[1:]:
		if arg in files: break
		cmd += ' ' + arg
	log.prNot(log.NOTICE, "Command: '%s'" % (cmd))
	log.prNot(log.NOTICE, "Parameters: '%s'" % (str(params)))
	log.prNot(log.NOTICE, "Files prefix: '%s'" % (os.path.commonprefix(files)))
	#log.prNot(log.NOTICE, "Files: '%s'" % (str(files)))
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
	elif (tool == 'tomo'): TomoTool(files, params)
	elif (tool == 'sdimm'): SdimmTool(files, params)
	# done
	dur = time.time() - __begin
	log.prNot(log.NOTICE, "Completed in %g seconds." % (dur))
	return 0


## @brief Process command-line options. See help_message for possible options.
def parse_options():
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
	
	opts, args = getopt.getopt(argv[2:], "vhsi:o:d:f:r:n:l:", ["verbose", "help", "stats", "informat=", "ff=", "fm=", "df=", "dm=", "mf=", "outformat=", "intclip=", "crop=", "file=", "dir=", "scale=", "rad=", "shape=", "pitch=", "xoff=", "origin", "noorigin", "disp=", "plot", "noplot", "norm", "nonorm", "saifac=", "range=", "sffile=", "safile=", "nref=", "mask=", "sfsize=", "sasize=", "overlap=", "border=", "offsets=", "subap=", "shifts=", "shifts-n=", "shifts-range=", "log=", "skip=", "ccdres=", "aptr=", "layerheights=", "nheights=", "layercells=", "skipsa=", "intpl=", "comp="])
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
		if option in ["--mask"]: params['mask'] = value
		if option in ["--comp"]: params['comp'] = value
		if option in ["--intpl"]: params['intpl'] = value
		# Tomo
		if option in ["--ccdres"]: params['ccdres'] = float(value)
		if option in ["--aptr"]: params['aptr'] = float(value)
		if option in ["--nheights"]: params['nheights'] = value.split(',')
		if option in ["--layerheights"]: params['layerheights'] = value.split(',')
		if option in ["--layercells"]: params['layercells'] = value.split(',')
		# sdimm options
		if option in ["--skipsa"]: params['skipsa'] = value.split(',')		
		if option in ["--shifts-n"]: params['shifts-n'] = int(value)
		if option in ["--shifts-range"]: params['shifts-range'] = value.split(',')
		
	
	return (tool, params, files)


## @brief Set default parameters for 'tool'
# @param [in] tool The tool class to get the default parameters for
# @return A dict holding the default parameters
def get_defaults(tool):
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
		default['comp'] = 'adsq'
		default['intpl'] = '9pt'
		default['mask'] = 'none'
	elif (tool == 'sdimm'):
		default['nref'] = 0
		#default['shifts-n'] = 1
		#default['shifts-range'] = ['0','1000']
	elif (tool == 'tomo'):
		default['ccdres'] = 0.45
		default['aptr'] = 0.49
		default['shifts'] = False
		default['nheights'] = ['3','1']
		default['layerheights'] = ['0-1000', '10000-10000']
		default['layercells'] = ['8', '8']
	
	return default


## @brief Check if the parsed parameters are valid
def check_params(tool, params):
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
	
	if params.has_key('mask') and \
		(params['mask'] not in ['circular', 'none']):
		log.prNot(log.ERR, "mask invalid, should be 'circular' or 'none'.")
		
	if params.has_key('skipsa'):
		params['skipsa'] = N.array(params['skipsa'], dtype=N.int)
	
	if params.has_key('shifts-range'):
		try: 
			params['shifts-range'] = N.array(params['shifts-range'], \
		 		dtype=N.int)
			tmp = params['shifts-range'].reshape(-1,2)
		except: log.prNot(log.ERR, "shifts-range invalid, should be Nx(<int>,<int>).")
	
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
		#if (not params['shifts']) or (not os.path.exists(params['shifts'])):
		#	log.prNot(log.ERR, "Tool 'shiftoverlay' requires shifts file.")
		if (params['shape'] not in ['box', 'dot']):
			log.prNot(log.ERR, "'shape' should be in ['box', 'dot'].")
	elif (tool == 'shifts'):
		# need safile and sffile
		if (not params['safile']) or (not os.path.exists(params['safile'])):
			log.prNot(log.ERR, "Tool 'shifts' requires safile.")
		if (not params['sffile']) or (not os.path.exists(params['sffile'])):
			log.prNot(log.ERR, "Tool 'shifts' requires sffile.")
		if (params['intpl'] not in ["9pt", "5pt"]):
			log.prNot(log.ERR, "Tool 'shifts' intpl methods are '9pt' or '5pt'.")
		if (params['comp'] not in ["adsq", "sqd"]):
			log.prNot(log.ERR, "Tool 'shifts' intpl methods are 'sqd' or 'adsq'.")
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
	elif (tool == 'sdim'):
		# cannot use shifts-range and shifts-part simultaneously
		if (params['shifts-range']) and (params['shifts-n']):
			log.prNot(log.ERR, "Tool 'sdimm' cannot use 'shifts-range' and 'shifts-n' simultaneously.")
	# Done


### ==========================================================================
### Tool classes
### ==========================================================================

## @brief Generic Tool class with common functions.
# 
# This is the base class for all other 'tools' provided by astooki. This class 
# provides some functions that are used by most/all other tools that subclass 
# this class. Functionality includes loading and saving image data, 
# dark/flatfidling and cropping & masking data.
class Tool(object):
	def __init__(self, files, params):
		# Init starttime
		## @brief Starttime for the tool class
		self.start = time.time()
		## @brief List of files passed on to astooki (can be empty)
		self.files = files
		## @brief Number of files
		self.nfiles = len(files)
		## @brief Raw unparsed parameters dict
		self.params = params
		# Save some options common for all tools
		## @brief Output directory for files
		self.outdir = params['outdir']
		## @brief Input format for any files read
		self.informat = params['informat']
		## @brief Flatfield file to use (if any)
		self.flatfield = params['flatfield']
		## @brief Number of frames summed in flatfield
		self.flatmulti = params['flatmulti']
		self.flatdata = None
		## @brief Darkfield file to use (if any)
		self.darkfield = params['darkfield']
		## @brief Number of frames summed in darkfield
		self.darkmulti = params['darkmulti']
		self.darkdata = None
		self.gaindata = None
		## @brief Subaperture maskfile to use
		self.maskfile = params['maskfile']
		self.mask = None
		## @brief Toggle normalisation of output
		self.norm = params['norm']
		## @brief Toggle cropping of processed files
		self.crop = params['crop']
		## @brief Toggle output of PDF plots of results in certain cases
		self.plot = params['plot']
		## @brief Output filename to use
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
		## @brief Dict of output files created 
		self.ofiles = {}
		self.ofiles['params'] = lf.saveData(self.mkuri('astooki-params'), \
		 	params, aspickle=True)
		
	
	
	## @brief Load a file from disk.
	#  Load a file from disk using self.informat as format, and crop it if 
	#  self.crop is set.
	def load(self, filename):
		log.prNot(log.INFO, "Loading '%s'" % (filename))
		if not os.path.isfile(filename):
			log.prNot(log.WARNING, "'%s' is not a regular file." % (filename))
			return None
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
	
	
	## @brief Load an ANA file
	def __anaload(self, filename):
		import pyana
		return pyana.getdata(filename)
	
	
	## @brief Load a FITS file
	def __fitsload(self, filename):
		import pyfits
		return pyfits.getdata(filename)
	
	
	## @brief Load a NumPy file
	def __npyload(self, filename):
		import numpy
		return numpy.load(filename)
	
	
	## @brief Dark- and flatfield data. 
	#
	#  Initialize dark and gain if necessary.
	def darkflat(self, data):
		if (self.flatfield or self.darkfield):
			self.__initdarkflat()
			tmp = (data-self.darkdata) * self.gaindata
			log.prNot(log.INFO, "Dark-flatfielding data, range: %.4g -- %4g, avg: %.4g" % (N.min(tmp), N.max(tmp), N.mean(tmp)))	
			return tmp
		return data
	
	
	## @brief Initialize dark and flatfield if not already done
	#
	#  Initialize dark- and flatfield if self.darkdata and/or self.flatdata is
	#  not set already. To speed up dark-/flatfielding, create self.gaindata
	#  which is 1/(flat-dark) (with some filtering to prevent division by 0)
	def __initdarkflat(self):
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
	
	
	## @brief Save data as FITS file to disk
	#  @param data Data to save
	#  @param filepath File to save to
	#  @param overwrite Overwrite file (explicitly needed for fitssave)
	def fitssave(self, data, filepath, overwrite=True):
		import pyfits
		log.prNot(log.INFO, "Tool.fitssave(): Saving data to '%s'." % (filepath))
		pyfits.writeto(filepath, data, clobber=overwrite)
	
	
	## @brief Save data as ANA file to disk
	#  @param data Data to save
	#  @param filepath File to save to
	#  @param compressed Toggle compression on or off
	def anasave(self, data, filepath, compressed=1):
		import pyana
		log.prNot(log.INFO, "Tool.anasave(): Saving data to '%s'." % (filepath))
		pyana.fzwrite(filepath, data, compressed)
	
	
	## @brief Save data as NumPy file to disk
	#  @param data Data to save
	#  @param filepath File to save to
	def npysave(self, data, filepath):
		log.prNot(log.INFO, "Tool.npysave(): Saving data to '%s'." % (filepath))
		N.save(data, filepath)
	
	
	## @brief Save data as PNG file to disk
	#  @param data Data to save
	#  @param filepath File to save to
	#  @param scale Scale values from 0 to 255 before saving
	def pngsave(self, data, filepath, scale=True):
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
	
	
	## @brief Make an URI given from 'path'
	def mkuri(self, path):
		_path = os.path.basename(path)
		if (_path != path):
			log.prNot(log.WARNING, "mkuri(): got path instead of filename.")
		
		return os.path.join(self.outdir, _path)
	
	
	## @brief Apply a mask to data
	#
	#  Apply mask to data, setting all values outside the mask to the maximum
	#  value of the pixels within the mask. If self.norm is set, also normalize 
	#  the pixel values within each submask 
	#  @param data Data to mask
	def maskimg(self, data):
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
	

	## @brief Initialize mask for maskimg()
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
		
	


## @brief Make a subaperture mask
# 
# 
# This tool calculates a set of coordinates where subapertures or subimages 
# are located and saves these to disk. Optionally, the pattern can be plotted
# as well. The units for the parameters are arbitrary, but should be 
# consistent. For subimages on a CCD one would use pixels, for subapertures at
# the lenslet, one would use meters.
# 
# The difference between subapertures and subimages is the plane they are 
# located in. Subapertures are located in the aperture plane (i.e. on the 
# lenslet) while subimages are located on the CCD. For this tool this is 
# obviously irrelevant, but the meaning can be quite different.
#
# The following parameters should be defined in the Tool.params dict passed to 
# the initializer of this class:
# 
# @param rad Telescope aperture radius
# @param shape Shape of the telescope aperture (circular or square)
# @param sasize Subapertures size
# @param pitch Pitch of the subaperture positions
# @param xoff x-offset of even and odd rows in units of sasize. [0,0.5] will
# 	give a hexagonal brick-pattern grid while [0, 0] will give a square grid
# @param scale Scale the whole grid by this factor
# @param disp Displace the whole grid by this vector
# @param file Base filename to save the pattern to
class SubaptConfTool(Tool):
	def __init__(self, files, params):
		super(SubaptConfTool, self).__init__(files, params)
		## @brief Subaperture pattern radius
		self.rad = params['rad']
		## @brief Subaperture size
		self.sasize = params['sasize']
		## @brief Subaperture pitch
		self.pitch = params['pitch']
		## @brief Aperture shape, either circular or square
		self.shape = params['shape']
		## @brief Row offset
		self.xoff = params['xoff']
		## @brief Scale factor for complete pattern
		self.scale = params['scale']
		## @brief Displacement vector for complete pattern
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
		
	

## @brief Make a subfield mask
#
# This is more or less the same as SubaptConfTool(), except that this tool
# generates a subfield mask and takes slightly different parameters. The 
# subfield mask will be a set of (pixel) coordinates relative to the origin of 
# a subaperture that define crops of the subimages on the CCD. Since pixels in 
# a subimage corresponds to a different field of view and view angle, a grid 
# of masks for a subimage correspond to a grid of different fields of view. 
# This grid can be used to calculate the subfield shifts that can consequently 
# be used to analyze the seeing recorded in the data.
# 
# The grid generated is relative to the subimage, so the size of such a 
# subimage is relevant: subfields cannot lie *outside* this subimage. Also, 
# when later measuring the shifts for each subfield within a subimage, one 
# will have to move the subfield around to measure the cross-correlation. 
# Therefore one should supply a border which will take care of this and keep a 
# guard range free at the edge of the subimage. The shift range used later 
# should always be less or equal to the border supplied here.
# 
# @param sasize The subaperture size to be used
# @param sfsize The subfield size
# @param overlap The overlap in x and y direction. 1 for complete overlap, 0 
#   for no overlap. [0.5, 0.5] gives about 50% overlap. Note that this might 
#   be changed slightly because of the granularity of pixels
# @param border The border or guard range to keep clear withins sasize
class SubfieldConfTool(Tool):

	def __init__(self, files, params):
		super(SubfieldConfTool, self).__init__(files, params)
		## @brief Subaperture size
		self.sasize = params['sasize']
		## @brief Subfield size
		self.sfsize = params['sfsize']
		## @brief Overlap between different subfields
		self.overlap = params['overlap']
		## @brief Border to keep clear around the subfield pattern
		self.border = params['border']
		
		if params['file']: self.file = params['file']
		else: self.file = 'astooki-subfieldconf.csv'
		
		self.run()
	
	
	def run(self):
		# Generate subfield positions
		effsize = self.sasize - 2*self.border
		pitch = N.round(self.sfsize * (1-self.overlap)).astype(N.int)
		nsf = N.floor((effsize - self.sfsize) / pitch) + 1
		
		log.prNot(log.INFO, "Eff sasize: (%g,%g), pitch: (%g,%g), eff overlap: (%g,%g)" % \
		 	(tuple(effsize) + tuple(pitch) + tuple(1-pitch*1.0/self.sfsize)))
		
		sfpos = self.border + \
			N.indices(nsf, dtype=N.int).reshape(2,-1).T * pitch
		sfpos = N.floor(sfpos).astype(N.int)
		totnsf = N.product(nsf).astype(N.int)
		
		log.prNot(log.NOTICE, "Found %d x %d subfields." % tuple(nsf))
		log.prNot(log.NOTICE, "Size %d,%d" % tuple(self.sfsize))
		libsh.saveSaSfConf(self.mkuri(self.file), totnsf, [-1,-1], self.sfsize, \
		 	sfpos)
		if (self.plot):
			import astooki.libplot as lp
			lp.showSaSfLayout(self.mkuri(self.file+'-plot.eps'), sfpos, \
			 	self.sfsize, plrange=[[0, self.sasize[0]], [0, self.sasize[1]]])
		# Done
		
	

## @brief Optimize subaperture configurations on real data.
#
# Although a subaperture mask made with SubaptConfTool() will be fairly 
# accurate if it is supplied with the right parameters, it is still possible 
# that the mask does not match the subimages exactly. To circumvent this 
# problem, SubaptOptTool() can take the statically generated subimage mask and 
# a flatfield image and tries to match the mask onto the flatfield.
# 
# The method used here is as follows. Given the mask, at each centroid grid
# position take a vertial slice out of the flatfield that is ~30 pixels wide 
# and twice the subimage high. This slice should cover the whole flatfielded 
# subimage. This slice is then averaged over the width, giving a 1 pixel wide 
# profile vertically across the subimage flatfield. The first pixel to the top 
# and bottom from the center of the slice where the intensity is lower than X 
# times the maximum intensity of the slice is considered to be the edge of the 
# subimage. This will give the height of the subimage. The same routine is 
# also done for horizontally across the subimage to get the width.
# 
# The drawback of this routine is that is needs the flatfield and real image 
# to match perfectly. Fortunately, this should and ususally is that case. 
# Furthermore, if there is a speck of dust on the flatfielded image, this can 
# fool this tool. One should therefore always check the output, for example by 
# plotting the grid or testing the grid by using it as a mask on real data.
# 
# @param maskfile The subaperture grid to optimize
# @param saifac The intensity dropoff factor to use (X in the above info)
# @param rad The radius of the subaperture grid pattern (used for plotting 
# 	only)
class SubaptOptTool(Tool):

	def __init__(self, files, params):
		super(SubaptOptTool, self).__init__(files, params)
		# Output file
		tmp = os.path.splitext(self.maskfile)
		if params['file']: self.file = params['file']
		else: self.file = tmp[0]+'-optimized'+tmp[1]
		## @brief The intensity dropoff factor to use
		self.saifac = params['saifac']
		## @brief Radius of the mask pattern (plotting only)
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
	


## @brief Calculate image-shifts between different subfields and subimages
#
# Processing the output of a (Shack-Hartmann) wavefront sensor starts with 
# calculating the image shifts between various subimages and subfields. Given 
# a set of subimage and subfield coordinates, this tool calculats the image 
# shifts for those subimages/subfields.
# 
# The image shifts are calculated per frame. For each frame, the <nref> 
# subimages with the highest RMS value are selected as reference subimages. 
# Each reference subimage is compared with all subimages (also with itself). 
# For each pair of subimages, the image shift is calculated for all subfields. 
# In total this yields a 5-dimensional N_files * N_references * N_subimages 
# * N_subfields * 2 matrix of data. This data is stored in FITS and NumPy 
# binary format in 'image-shifts.<fits|npy>' as 32-bit float values. Note that 
# these files will become rather big, a typical set of 1000 frames with 2 
# references, 85 subimages and about 270 subfields will already result in 
# 350 megabytes of data.
# 
# Calculating the image shifts is done in two steps. First a comparison map is 
# made as function of the shift vector. This is done by comparing a reference 
# subimage with a slightly shifted subimage. The method to compare the image 
# with the reference can be a cross-correlation, the absolute difference 
# squared, or the squared difference of the two images.
#
# Optionally, the two subimages can be masked before comparing them. This can 
# be useful to synchronize the data processing with theory. The only mask
# currently implemented is circular mask.
# 
# After the comparison map is calculated, the maximum value indicates the  
# shift which matches the image with the reference best. To get a subpixel 
# shift vector, the maximum is interpolated from a few pixels around the 
# maximum value. This can be done with a 9-point or a 5-point quadratic 
# interpolation.
# 
# @param shrange Shiftrange to measure in pixels (actual range will be 
# 	from -shrange to +shrange)
# @param safile Subimage mask on the CCD
# @param sffile Subfield mask on the CCD
# @param nref Number of reference subimages to use
# @param comp Comparison method to use. Can be 'adsq' for absolute difference 
#   squared or 'sqd' for square difference.
# @param inpl Subpixel interpolation method to use. Can be '5pt' for 5-point
#   quadratic interpolation or '9pt' for 9-point interpolation.
# @param mask Mask to use when comparing images ('circular' or 'none')
class ShiftTool(Tool):

	def __init__(self, files, params):
		super(ShiftTool, self).__init__(files, params)
		## @brief Shift range to use
		self.shrange = params['shrange']
		## @brief Subimage mask file
		self.safile = params['safile']
		## @brief Subfield mask file
		self.sffile = params['sffile']
		## @brief Number of subimages to use as reference
		self.nref = params['nref']
		## @brief Mask to use before measuring shifts
		self.mask = params['mask']
		# Load safile and sffile
		(self.nsa, self.saccdpos, self.saccdsize) = \
			libsh.loadSaSfConf(self.safile)
		(self.nsf, self.sfccdpos, self.sfccdsize) = \
			libsh.loadSaSfConf(self.sffile)
		if (len(self.files) < 1):
			log.prNot(log.ERR, "ShiftTool(): Cannot continue without files!")
		
		# Parse image comparison and interpolation methods
		import astooki.clibshifts as ls
		## @brief Comparison method to use (abs diff squared, diff squared)
		self.comp = ls.COMPARE_ABSDIFFSQ
		if params['comp'] == 'adsq': self.comp = ls.COMPARE_ABSDIFFSQ
		elif params['comp'] == 'sqd': self.comp = ls.COMPARE_SQDIFF
		else: log.prNot(log.ERR, "ShiftTool(): 'comp' parameter invalid!")
		
		## @brief Subpixel interpolation method to use (9pt or 5pt)
		self.intpl = ls.EXTREMUM_2D9PTSQ
		if params['intpl'] == '9pt': self.intpl = ls.EXTREMUM_2D9PTSQ
		elif params['intpl'] == '5pt': self.intpl = ls.EXTREMUM_2D5PTSQ
		else: log.prNot(log.ERR, "ShiftTool(): 'intpl' parameter invalid!")
		
		# Run analysis
		self.run()
	
	
	def run(self):
		### Phase 1: Measure shifts
		### -----------------------
		
		self.calcShifts()
		
		### Phase 2: Process shifts
		### -----------------------
		
		log.prNot(log.NOTICE, "Starting post-processing.")
		self.postProcess()
		
		### End: store data & metadata
		### FIXME: saving 'shifts' has problems when it is a lot of memory!
		self.ofiles['refaps'] = lf.saveData(self.mkuri('referenace-subaps'), \
			self.allrefs, asnpy=True, ascsv=True)
		self.ofiles['files'] = lf.saveData(self.mkuri('processed-files'), \
		 	self.allfiles, asnpy=True, ascsv=True, csvfmt='%s')
		self.ofiles['shifts'] = lf.saveData(self.mkuri('image-shifts'), \
		 	self.shifts, asnpy=True, asfits=True)
			
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
			pass
		
		try:
			self.shifts = N.load(self.mkuri('image-shifts.npy'))
			if (self.shifts.shape == \
				(self.nfiles, self.nref, self.nsa, self.nsf,2)): 
				log.prNot(log.NOTICE, "Found previously measured shifts, restoring.")
				return
			else: raise Exception
		except:
			pass
		
		import astooki.clibshifts as ls
		
		# Convert mask
		if self.mask == 'circular': self.mask = ls.MASK_CIRC
		else: self.mask = None
		
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
			
			spf = (time.time() - self.start) / len(allfiles)
			togo = (len(self.files) - len(allfiles)) * spf
			eta = time.localtime(time.time() + round(togo))
			log.prNot(log.NOTICE, "Processing frame %d/%d, ETA @ %s, (%g sec, %g spf)" % \
				(len(allfiles), len(self.files), time.strftime("%H:%M:%S", eta), togo, spf))
			refaps = []
			imgshifts = ls.calcShifts(dfimg, self.saccdpos, self.saccdsize, \
			 	self.sfccdpos, self.sfccdsize, method=self.comp, \
			 	extremum=self.intpl, refmode=ls.REF_BESTRMS, \
			 	refopt=self.nref, shrange=[self.shrange, self.shrange], \
				mask=self.mask, refaps=refaps)
			allrefs.append(refaps)
			allshifts.append(imgshifts)
		
		log.prNot(log.NOTICE, "Done, saving results to disk.")
		self.shifts = N.array(allshifts)
		self.allrefs = allrefs
		self.allfiles = allfiles
	
	
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
		
		# Store variance maps to find bad subapertures/subfields
		var = N.var(N.mean(self.shifts,1), 0)
		self.varmap = N.r_[[var[...,0], var[...,1]]]
		self.ofiles['variance'] = lf.saveData(self.mkuri('shifts-variance'), \
			 self.varmap, asnpy=True, asfits=True)
		self.ofiles['saccdpos'] = lf.saveData(self.mkuri('subap-ccdpos'), \
		 	self.saccdpos, asnpy=True, asfits=True)
		self.ofiles['sfccdpos'] = lf.saveData(self.mkuri('subfield-ccdpos'), \
		 	self.sfccdpos, asnpy=True, asfits=True)
		self.ofiles['saccdsize'] = lf.saveData(self.mkuri('subap-ccdsize'), \
		 self.saccdsize, asnpy=True, asfits=True)
		self.ofiles['sfccdsize'] = lf.saveData(self.mkuri('subfield-ccdsize'), \
		 	self.sfccdsize, asnpy=True, asfits=True)
		self.cpos = self.saccdpos.reshape(-1,1,2) + \
		  self.sfccdpos.reshape(1,-1,2) + self.sfccdsize.reshape(1,1,2)/2.0
		self.ofiles['sasfpos-c'] = lf.saveData(self.mkuri('sasfpos-c'), \
			self.cpos, asnpy=True, asfits=True)
		# Add meta info
		self.ofiles['path'] = os.path.dirname(os.path.realpath(self.mkuri('.')))
		
		if (self.nsf == 1):
			# First average over all Nref reference subapertures
			s_ref = N.mean(self.shifts[:,:,:,0,:], axis=1)
			# Now make sure the average *per frame* is zero
			s_avgfr = N.mean(s_ref, axis=1)
			s_norm = s_ref - s_avgfr.reshape(-1,1,2)
			
			# Calculate variance per subaperture
			savar = N.var(s_norm, axis=0)
			log.prNot(log.NOTICE, "Average shift variance: (%g,%g), max: (%g,%g)" %\
				(tuple(N.mean(savar,0)) +  tuple(N.max(savar,0))))
			
			# Now average over all frames to get the offset + error
			log.prNot(log.NOTICE, "Calculating static offsets.")
			soff = N.mean(s_norm, axis=0)
			sofferr = (N.var(s_norm, axis=0))**0.5
			
			self.ofiles['offsets'] = lf.saveData(self.mkuri('static-offsets'), \
				soff, asnpy=True, ascsv=True)
			self.ofiles['offset-err'] = lf.saveData(\
				self.mkuri('static-offset-err'), sofferr, asnpy=True, ascsv=True)
			if (self.plot):
				import astooki.libplot as libplot
				# TODO: fix this range -- make it a parameter?
				plran = [[0, N.ceil(self.saccdpos.max()/512.0)*512]]*2
				libplot.plotShifts(self.mkuri('static-offset-plot'), self.shifts, \
					self.saccdpos, self.saccdsize, self.sfccdpos, self.sfccdsize, \
					mag=7.0, allsh=False, \
				 	title='Static offsets, mag=7', legend=True)
				libplot.plotShifts(self.mkuri('static-offset-plot-250'), \
				 	self.shifts[:250], self.saccdpos, self.saccdsize, self.sfccdpos, \
				 	self.sfccdsize, mag=7.0, \
				 	allsh=False, title='Static offsets, mag=7', legend=True)
	


## @brief Update subimage mask with an offset.
#
# Since we want to compare different subfields within each subimage, we need 
# to know the reference direction of each subimage. Because of static 
# aberrations (telescope defocus, instrument issues) we cannot assume that 
# pixel (x,y) in subimage N corresponds to the same field of view as the same 
# pixel in subimage M. Or the other way around: given a granule on the sun, 
# we want to know at what pixel that granule is located in each of the 
# subimages.
# 
# To get these 'static offsets', we take a large field of view in one 
# reference subimage (almost the complete subimage) and cross-correlate that 
# with all subimages. This will give N shift vectors for each frame, N being 
# the number of subapertures. To get better results, it is possible to use 
# multiple subapertures as reference, this should give the same data and gives 
# an indication of the noise or reliability of the shift measurement. The 
# image shifts measured for different reference subapertures will be stored 
# alongside eachother.
# 
# Because there is atmospheric seeing which causes tip-tilt of the images, the 
# image shifts will vary strongly from frame to frame. The seeing is 
# statistical, however, which means that the average shift over all frames 
# should be zero. To get the static offsets and identify the relative field of 
# view of all subimages, we average the image shifts over all frames. This 
# gives a offset vector for each subimage which indicates where the subimage 
# is pointing at relative to the other subimages.
# 
# If we correct the subimage mask calculated earlier with this list of 
# vectors, each subimage mask will be pointing at the same location on the 
# sun. Once we have established this, we can subdivide the subimages in 
# different subfields knowing exactly where each subfield points and thus 
# providing a reliable method to base subsequent analysis on.
# 
# @param maskfile The subimage mask to be updated
# @param offsets A list of offset vectors
class SubaptUpdateTool(Tool):

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
		newpos = (pos - off + maxsh).astype(N.int)
		# Crop the subaperture size by twice the maximum offset
		newsize = (size - maxsh*2).astype(N.int)
		# Store
		libsh.saveSaSfConf(self.mkuri(self.file), nsa, [-1,-1], newsize, newpos)
		if (self.plot):
			import astooki.libplot as libplot
			plfile = self.mkuri(self.file + '-plot.eps')
			# TODO: fix this range -- make it a parameter?
			plran = [[0, N.ceil(pos.max()/512)*512]]*2
			libplot.showSaSfLayout(plfile, newpos, newsize, plrange=plran)
		# Done
	

## @brief Calculate statistics on files
#
# Calculate some basic statistics on the files passed to this tool, like 
# min, max, average, rms etc.
class StatsTool(Tool):
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
	


## @brief Convert files to different format.
# 
# This tool is similar to the 'convert' program in the Imagemagick suite, 
# except it focuses on astronomical dataformats. It can take all 
# _INFORMATS and save to all _OUTFORMATS (both defined in pyatk.py)
#
# Note that PNG is limited to 2 dimensional data, while the other formats can 
# hold up to 8-dimensional data.
# 
# This tool can also scale, crop and clip images before saving them again,
# making it a good first start for making movies.
#
# @param file If there is only one file to process, save the output here
# @param outformat Use this format to save the converted frames
# @param scale Scale the dimensions of the frame by this factor
# @param intclip Clip the intensity to this range
# @sa Tool.crop Tool.maskfile Tool.darkfield Tool.flatfield
class ConvertTool(Tool):
	def __init__(self, files, params):
		super(ConvertTool, self).__init__(files, params)
		## @brief Save output here, if and only if there is only one input file
		self.file = params['file']
		## @brief Output format to use
		self.outformat = params['outformat']
		## @brief Scale dimensions of the frame before saving
		self.scale = params['scale']
		## @brief Clip the intensity to this range before saving
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
		
	


## @brief Crop out certain subimages, overlay shift vectors, output PNGs
#
# This tool can make a visual representation of the shift vectors measured 
# using ShiftTool. It crops out one or more subimages from the raw frames and 
# overlays a box displaced by the shift vector calculated before. This can be 
# useful to show the image shifts and make a movie out of it.
#
# In addition to the above, this tool can also scale and clip the frames. This  
# is done the same way as in ConvertTool.
#
# @param scale Scale the dimensions of the (cropped) output by this factor
# @param intclip Clip the intensity to this range
# @param shape Determines the shape of the overlay drawn. Can be 'box' for a 
#    box the size of the subfield used to measure the shift, or 'dot' for a 
#    centroid dot a few pixels big.
# @param skip Skip this many frames in the shifts file. By default the first 
#    frame is matched to the first set of shifts in the shifts file. If the 
#    first frame is in fact number 100 in a series though, one should skip 100 
#    indices in the shifts file.
# @param safile Subimage mask on the CCD
# @param sffile Subfield mask on the CCD
# @param shifts Image shifts measured with ShiftTool
class ShiftOverlayTool(Tool):
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
		
		if os.path.exists(self.shifts):
			haveShifts = True
			allshifts = lf.loadData(self.shifts, asnpy=True)
		else:
			haveShifts = False
		
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
			if haveShifts:
				log.prNot(log.NOTICE, "Pre-processing shifts.")
				shifts = allshifts[fidx+self.skip,:,:,:,:]
				notfin = N.argwhere(N.isfinite(shifts) == False)
				for nfidx in notfin: shifts[tuple(nfidx)] = 0.0
				log.prNot(log.NOTICE, "Setting %d non-finite values to 0." % \
					(len(notfin)))
			
			# Average over different references
			if haveShifts:
				shifts = N.mean(shifts, axis=0)
			
			# Loop over subaps to process
			for sa in self.subaps:
				# Crop out subaperture, divide by mean
				data = dfimg[\
					self.saccdpos[sa, 1]: \
					self.saccdpos[sa, 1] + self.saccdsize[1], \
					self.saccdpos[sa, 0]: \
					self.saccdpos[sa, 0] + self.saccdsize[0]]
				data = data / data.mean()
				
				if (self.scale != 1.0):
					# If scale is not unity, calculate the nearest scale that is % 8
					import scipy.ndimage
					sc = self.scale
					orig = data.shape
					nsc = (N.round(orig[1] * sc / 8) * 8)/orig[1]
					data = S.ndimage.zoom(data, nsc, mode='wrap')
					log.prNot(log.NOTICE, "Scaling image by %g (from %d,%d to %d,%d)."%\
					 	(nsc, orig[1], orig[0], data.shape[1], data.shape[0]))
					nsc = (N.array(data.shape)*1.0/N.array(orig))[::-1]
				else: nsc = 1.0
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
						# Visualize shift with moving subfields-sized boxes
						if haveShifts: shpos = (self.sfccdpos[sf] - shifts[sa, sf]) * nsc
						else: shpos = self.sfccdpos[sf] * nsc
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
						if haveShifts: vec = (self.sfccdposc[sf] - shifts[sa, sf]) * nsc
						else: vec = self.sfccdposc[sf] * nsc
						data[vec[1]-1:vec[1]+2, vec[0]-1:vec[0]+2] = dmax
						data[vec[1], vec[0]] = dmin
				
				# Save frame with shifts overlayed
				savefile = f+'-subap%d-shifts.%s' % (sa, self.outformat)
				savefile = self.mkuri(savefile)
				log.prNot(log.NOTICE, "Saving '%s' to '%s'." % \
					(base, os.path.basename(savefile)))
				
				if (self.outformat == _FORMAT_PNG): self.pngsave(data, savefile)
				elif (self.outformat == _FORMAT_FITS): self.fitssave(data, savefile)
				elif (self.outformat == _FORMAT_ANA): self.anasave(data, savefile)
				elif (self.outformat == _FORMAT_NPY): self.npysave(data, savefile)
		
	


## @brief Perform SDIMM+ analysis on shift data.
# 
# - Loop over all rows of subapertures (subapertures that are at the same y 
#   position)
# - For each subap row choose a reference subaperture (i.e. the left-most one)
# - Loop over all rows of subfields (subfields with the same y coordinate)
# - For each subfield row, choose a reference subfield
# - Loop over all other subapertures in this subap row
# - Loop over all other subfields in this subfield row
# - Compare all subaperture-subfield pairs as described in Scharmer & van 
#   Werkhoven
# - Repeat this for all columns
# 
# This tool outputs the raw results of the calculations to sdimm-N.<fits|npy> 
# and is a is an N * 9 matrix where each row holds the following information:
#    [id, s, a, C_lsa, C_tsa, refsa, sa, refsf, sf]
#  with:
#  - id=0 for row-wise comparison and 1 for column-wise (as described 
#    above),
#  - s the scalar distance between the two subapertures in meters,
#  - a the scalar angle between the two subfields in pixels (convert with the 
#    CCD scale to get a real angle), 
#  - C_lsa the longitudinal covariance between the two sequences of 
#    differential image shifts (as described in the paper)
#  - C_tsa the transversal covariance
#  - refsa the index of the reference subaperture used here
#  - refsf the ubdex of the reference subfield used here
#  - sa the index of the other subaperture used
#  - sf the index of the other subfield used
#  - N is the batch being processed, see 'shifts-n' below
#  
#  Besides the raw information, this tool also outputs processed information
#  to sdimmcol* and sdimmrow* files, where the *col* files have information on 
#  the column-wise comparison and the *row* files on the row-wise comparison 
#  of the data.
# 
# sdimm<col|row>.* is a 3 x N x M matrix with N the number of unique
#  subaperture distances (s) and M the number of unique angles (a) for this 
#  data. For the 85 subaperture lenslet array at the SST, N is 5 for
#  column-wise comparison and 9 for row-wise comparison. This is the data that
#  should be decomposde in SDIMM basis described in the paper. 
# 	
#  The first N x M frame holds the longitudinal covariance, the second frame
#  holds the transversal covariance and the third frame holds the number of
#  covariances each specific cell was averaged over (i.e. given an (s,a)
#  coordinate, how many covariances were calculated?
# 
#  sdimm<col|row>-s.* hold the unique subaperture distances mentioned above, 
#  and sdimm<col|row>-a.* hold the unique subfield angles mentioned above.
# 
# The following parameters are required as input for this tool:
# @param shifts the file holding the shift measurements
# @param safile centroid subaperture positions at the lenslet [meter]
# @param sffile centroid subfield positions [pixel]
# @param skipsa Subapertures to skip in analysis (i.e. with high noise)
# @param nref Number of references to use (0 for max)
#
# The following parameters are optional:
# @param shifts-n allows for processing subsets of the shift measurements in 
# batches of 'shifts-n' frames, instead of processing all shift measurements 
# in one go.
class SdimmTool(Tool):
	
	def __init__(self, files, params):
		super(SdimmTool, self).__init__(files, params)
		log.prNot(log.NOTICE, "Starting SDIMM+ analysis of WFWFS data stored in '%s'." % (params['shifts']))
		
		# Load subaperture centroid positions
		(self.nsa, self.sapos, self.sasize) = \
		 	libsh.loadSaSfConf(params['safile'])
		# Load subfield pixel positions
		(self.nsf, self.sfccdpos, self.sfsize) = \
			libsh.loadSaSfConf(params['sffile'])
		## @brief Skip these subaps in the analysis
		self.skipsa = N.array(params['skipsa'], dtype=N.int)
		## @brief Number of references to use for analysis
		self.nref = params['nref']
		## @brief Load shift data here
		self.shifts = lf.loadData(params['shifts'], asnpy=True)
		
		## @brief Number of frames to use per sdimm analysis, allows to split up 
		#  series in smaller subsets
		nframes = self.shifts.shape[0]
		self.shiftsr = []
		if params.has_key('shifts-n'):
			for i in range(nframes/params['shifts-n']):
				self.shiftsr.append([i * params['shifts-n'], \
					(i+1) * params['shifts-n']])
		elif params.has_key('shifts-range'):
			self.shiftsr = params['shifts-range'].reshape(-1,2)
		else:
			self.shiftsr = [[0, nframes]]
		
		self.shiftsr = N.array(self.shiftsr)
		
		log.prNot(log.NOTICE, "Got %d frames, using intervals: %s" % (nframes, self.shiftsr.flatten()))
				
		self.run()
	
	
	def run(self):
		# Calculate the SDIMM+ covariance values
		import astooki.libsdimm as lsdimm
		# Loop over different subsets of the shift measurements
		for r in self.shiftsr:
			log.prNot(log.NOTICE, "Processing frames %d--%d now..." % (r[0], r[1]))
			# Calculate ROW-wise covariance maps
			(slist_r, alist_r, covmap_r) = lsdimm.computeSdimmCovWeave(\
				self.shifts[r[0]:r[1]], self.sapos, self.sfccdpos, refs=self.nref, \
				skipsa=self.skipsa, row=True, col=False)
			
			# Save covariance map to disk
			self.ofiles["sdimmrow-%d--%d" % (r[0], r[1])] = lf.saveData(\
				self.mkuri("sdimmrow-%d--%d" % (r[0], r[1])), covmap_r, asfits=True)
		
			# Calculate COLUMN-wise covariance maps
			(slist_c, alist_c, covmap_c) = lsdimm.computeSdimmCovWeave(\
				self.shifts[r[0]:r[1]], self.sapos, self.sfccdpos, refs=self.nref, \
				skipsa=self.skipsa, row=False, col=True)
			
			# Save covariance map to disk
			self.ofiles["sdimmcol-%d--%d" % (r[0], r[1])] = lf.saveData(\
				self.mkuri("sdimmcol-%d--%d" % (r[0], r[1])), covmap_c, asfits=True)
		
			# Combine ROW and COLUMN covariance maps
			(slist_a, alist_a, covmap_a) = lsdimm.mergeMaps([covmap_r, covmap_c], \
				[slist_r, slist_c], \
				[alist_r, alist_c])
			
			# Save maps to disk
			self.ofiles["sdimm-%d--%d" % (r[0], r[1])] = lf.saveData(\
				self.mkuri("sdimm-%d--%d" % (r[0], r[1])), covmap_a, asfits=True)
		
		# Save ROW s and a values to disk
		self.ofiles['sdimmrow-s'] = lf.saveData(self.mkuri('sdimmrow-s'), \
			slist_r, asfits=True, ascsv=True)
		self.ofiles['sdimmrow-a'] = lf.saveData(self.mkuri('sdimmrow-a'), \
			alist_r, asfits=True, ascsv=True)
		
		# Save COLUMN s and a values to disk
		self.ofiles['sdimmcol-s'] = lf.saveData(self.mkuri('sdimmcol-s'), \
			slist_c, asfits=True, ascsv=True)
		self.ofiles['sdimmcol-a'] = lf.saveData(self.mkuri('sdimmcol-a'), \
			alist_c, asfits=True, ascsv=True)
		
		# Save s and a values to disk
		self.ofiles['sdimm-s'] = lf.saveData(self.mkuri('sdimm-s'), \
			slist_a, asfits=True, ascsv=True)
		self.ofiles['sdimm-a'] = lf.saveData(self.mkuri('sdimm-a'), \
			alist_a, asfits=True, ascsv=True)
				
		# Save meta file to disk
		metafile = lf.saveData(self.mkuri('sdimm-meta-data'), \
			self.ofiles, aspickle=True)
	


## @brief Simulate WFWFS measurements using N discrete KL phase screens
class SimulShift(Tool):
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
		
		# This is not done yet
		log.prNot(log.ERROR, "Not implemented yet")
		
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
	
	
	


## @brief Tomographically analyze WFWFS data
# 
# A different analysis method compared to SDIMM+, this method assumes nothing 
# about the astmophere and tries to invert it using a simple linear model. 
# 
# Currently work in progress.
class TomoTool(Tool):
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
	


### ==========================================================================
### Helper functions
### ==========================================================================

## @brief Show usage information
#
# @param tool Tool to print help about
# @param out File descriptor to print
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


##
# @mainpage Astooki
# @author Tim van Werkhoven (tim@astro.su.se)
# @date 20090624
#
# @section About_sec About
# 
# Astooki, the Astronomical Toolkit is a Python script to process astronomical 
# data. Currently, it focuses strongly on WFWFS data. It can be obtained at 
# http://github.com/tvwerkhoven/astooki-py/
#
# @section Overview_sec Overview
# 
# Astooki provides the following tools to process or generate data with. The 
# name between brackets is the command line option to be used to call this 
# tool:
# - ConvertTool ('convert') to convert data/frames from one format to another
# - StatsTool ('stats') to calculate basic statistics on frames
# - SubaptConfTool ('samask') to generate subimage/aperture masks
# - SubfieldConfTool ('sfmask') to generate subfield masks
# - ShiftOverlayTool ('shiftoverlay') to visualize shift data
# - SubaptOptTool	 ('saopt') to optimize subimage/aperture masks using 
#     flatfield frames
# - SubaptUpdateTool ('saupd') to update subimage/aperture masks with a static 
#     offset
# - ShiftTool ('shifts') to calculate subimage/subfield image shifts 
# - SdimmTool ('sdimm') to analyze the shift data using the SDIMM+ 
#     method
#
# @section Issues_sec Issues
#
# No known issues at the moment
# 
# @section TODO_sec TODO
#
# - Add tomographic inversion analysis method
# - Add seeing simulation routines
# - Implement full SDIMM+ analysis method
# - Implement automatic data-reduction script/tool/framework
#
# @page proc_page Data processing
# @author Tim van Werkhoven (tim@astro.su.se)
# @date 20090624
#
# @section about_sec About
#
# This page gives a short quick & dirty howto on processing WFWFS data using 
# astooki. The datasets from 2009.06.10 are taken as an example, but any other 
# can be used instead.
#
# @section inspect_sec Inspecting the data
# 
# Using the ConvertTool we convert one image \c 
# "../summary/2009.06.10-run00/*best" to a FITS file and store it in \c 
# "img/":
#
# <pre>
# pyatk.py convert -d img -vv \
#  --file 2009.06.10-run00-best.fits \
#  ../summary/2009.06.10-run00/*best
# </pre>
# 
# This will convert one image to a FITS file using a darkfield of 300 summed 
# frames and a flatfield of 300 summed frames.
# 
# <pre>
# pyatk.py convert -vv -d img \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --file 2009.06.10-run00-best-df.fits \
#  ../summary/2009.06.10-run00/*best
# </pre>
# 
# Visually inspecting these files gives the following values for the subimage 
# size and pitch:
# 
# - Subimage size: 177,153
# - Subimage pitch: 193,167
#
# These are necessary to create a subimage mask in the next step which will be 
# used to all further data processing
#
# @section setup_sec Initial setup
# @subsection subimg_sec Generate a subimage mask
#
# With SubaptConfTool we generate a subimage mask that can be used later for 
# data processing. The following parameters are important here:
#
# - 'rad' is the radius of the pattern (1024 in our case with a 2k by 2k pixel
#     CCD), 
# - 'shape' is the shape of the global aperture, 
# - 'sasize' the size of the subimage, 
# - 'pitch' the pitch between two subimage,
# - 'xoff' is the horizontal offset for *even* and *odd* rows on the 
#     subimage pattern. Set these to 0 to get a tile pattern corresponding 
#     to a square lensley layout, setting it to 0.5, 0 gives a brick pattern 
#     corresponding to a hexagonal lenslet layout.
# - 'disp' is the displacement vector used. Since the image coordinates go 
#     from (0, 0) to (2048, 2048) we need to shift the pattern by (1024, 1024) 
#     pixels. We add an offset vector -(7,24) to approximately correct for 
#     alignment errors.
# - 'scale' is a global scale factor that can be used to scale the whole 
#     pattern.
# - 'plot' indicates that this mask should also be plotted to file
#  
# <pre>
# pyatk.py samask -vv -d simask \
#  --file simask.csv \
#  --rad 1024 \
#  --shape circular \
#  --sasize 175,154 \
#  --pitch 194,166 \
#  --xoff 0,0.5 \
#  --disp 1017,1000 \
#  --scale=1 \
#  --plot
# </pre>
#
# @subsection subimg2_sec Using a subimage mask
#
# Using the \c '--mf' option in ConvertTool, one can specify that a subimage 
# mask should be used when processing the image. This will crop out everything
# outside this mask. We use the newly generated mask to process an image to 
# see how well it fits. Since it is a statically generated mask it probably
# will not fit very well, this is solved in the next bit.
# 
# <pre>
# pyatk.py convert -vv -d img \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --mf simask/simask-origin.csv \
#  --file 2009.06.10-run00-best-df-mask.fits \
#  ../summary/2009.06.10-run00/*best
# </pre>
#
# @subsection optmask_sec Optimizing a subimage mask
#
# Using SubaptOptTool and a flatfield frame, we can optimize the mask to make 
# it fit better to a flatfield image. This is done as follows:
# 
# <pre>
# pyatk.py saopt -vv -d simask-fopt \
#  --mf simask/simask-origin.csv \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --file simask-orig-fopt.csv \
#  --saifac 0.8 \
#  --rad 1024 \
#  --plot
# </pre>
#
# Use the optimized mask on the data to see the improvement:
# 
# <pre>
# pyatk.py convert -vv -d img \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --mf simask-fopt/simask-orig-fopt.csv \
#  --file 2009.06.10-run00-best-df-mask-fopt.fits \
#  ../summary/2009.06.10-run00/*best
# </pre>
# 
#
# @section calib_sec Calibration
#
# Since we want to compare different subfields within each subimage, we need 
# to know the reference direction of each subimage. Because of static 
# aberrations (telescope defocus, instrument issues) we cannot assume that 
# pixel (x,y) in subimage N corresponds to the same field of view as the same 
# pixel in subimage M. Or the other way around: given a granule on the sun, we 
# want to know at what pixel that granule is located in each of the subimages.
# 
# See SubaptUpdateTool for more information.
#
# @subsection sfbig_sec Generate subfield mask
# 
# For the static offset correction of any possible telescope aberrations, we
# want to measure image shifts of the complete subimage. To do so, we generate 
# a subfield mask with one big subfield using SubfieldConfTool:
#
# <pre>
# pyatk.py sfmask -vv -d sfmask \
#  --file sfmask-big.csv \
#  --sfsize=113,90 \
#  --sasize=173,150 \
#  --overlap=0,0 \
#  --border=30,30 \
#  --plot
# </pre>
#
# @subsection statoff_sec Measure static offset shifts
#
# We now measure the image shifts for all subimages using the big subfield 
# mask generated above using ShiftTool:
#
# <pre>
# pyatk.py shifts -vv -d statoff-run00 \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --safile simask-fopt/simask-orig-fopt.csv \
#  --sffile sfmask/sfmask-big.csv \
#  --range 7 --nref 1 --mask none --plot \
#  ../2009.06.10-run00/wfwfs_survey_im*00100
# </pre>
#
# @subsection saupd_sec Update subimage mask with offsets
#
# We now use the previously calculated static offsets to update the subimage
# mask. After this update, each pixel in each subimage will point to the same 
# object on the sky. This is done with SubaptUpdateTool as follows:
#
# <pre>
# pyatk.py saupd -vv \
#  -d statoff-run00 \
#  --mf simask-fopt/simask-orig-fopt.csv \
#  --offset statoff-run00/static-offsets.csv \
#  --plot
# </pre>
#
# @section datared_sec Data reduction
#
# Now that we have a fitting subimage mask and determined the static offsets,
# it is time to measure the actual subfield/subimage shifts for the dataset.
#
# @subsection subf_ssec Generate subfield mask
#
# To do so, we first need a subfield mask with smaller subfields. To get 
# unified subfield masks for all datasets, give the *smallest* subimage size 
# for all datasets calculated in the previous step as \c --sasize parameter.
# This is not necessary, but might be desirable.
#
# <pre>
# pyatk.py sfmask -vv -d sfmask \
#  --file sfmask-16x16.csv \
#  --sfsize=16,16 \
#  --sasize=167,140 \
#  --overlap=0.7,0.7 \
#  --border=7,7
# </pre>
#
# @subsection shift_ssec Measure subfield image shifts
#
# Using the subimage- and subfield masks, calculate the image shift for each 
# subfield in each subimage. This data can later be inverted to produce data 
# about the seeing. We use a circular mask to mimic a circular field of view
# as closely as possible:
#
# <pre>
# pyatk.py shifts -vv \
#  -d subshift-16x16-run00 \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --safile statoff-2009.06.10-run00/simask-orig-fopt-updated.csv \
#  --sffile sfmask/sfmask-16x16.csv \
#  --range 7 --nref 2 --mask circular \
#  ../2009.06.10-run00/wfwfs_survey_im*
# </pre>
#
# @section sdimm_sec SDIMM+ analysis
#
# Using the shifts, we can invert the data to produce information about the
# actual seeing, which is what we're after. This can be done with either a 
# statistical SDIMM+ method, or a tomographic method. Here we describe the 
# SDIMM+ method.
#
# First make a subaperture mask of the *lenslets*. This should be a mask of 
# the lenslet coordinates which will later be used for the inversion of the 
# SDIMM+ covariance maps.
#
# <pre>
# pyatk.py samask -vv -d samask\
#  --file samask-ll.csv \
#  --rad 0.52 \
#  --shape circular \
#  --sasize 0.098,0.098 \
#  --pitch 0.098,0.0849 \
#  --xoff 0,0.5 \
#  --disp 0,0 \
#  --scale=1 --plot
# </pre>
#
# Note that radius of 0.52 meter is slightly larger than the real aperture 
# radius, this is because when generating the subaperture mask every 
# subaperture is required to fit 100% within this radius. In the real optical 
# setup, some subapertures are slightly cropped though.
#
# To invert the data using the SDIMM+ method, use:
#
# <pre>
# pyatk.py sdimm -vv \
#  -d sdimm-16x16-run00 \
#  --skipsa -1 \
#  --shifts subshift-16x16-2009.06.10-run00/image-shifts.npy \
#  --safile samask/samask-ll-centroid.csv \
#  --sffile sfmask/sfmask-16x16.csv
# </pre>
#
# The dataset \c run00 for day 2009.06.10 is now reduced to a simple SDIMM+
# covariance map and is ready in \c sdimm-16x16-run00/
#
# @section miscsec Miscellaneous
# 
# These are some miscellaneous recipes.
#
# @subsection moviesec Make a movie 
#
# Crop out a part of the whole frame:
# 
# <pre>pyatk.py convert \
#  -vv -d crop-sa31 \
#  -o png \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --mf simask-fopt/simask-orig-fopt.csv \
#  --scale 1.5 \
#  --intclip 0.8,1.2 \
#  --crop 652,1111,167,141 \
#  ../*_im* </pre>
# 
# Convert the PNG's to a movie:
# 
#<pre> mencoder mf://*.png -mf fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=1000 -o wfwfs_movie.avi</pre>
#
# @subsection movie2sec Make a movie with shift boxes
# Crop out, overlay shifts:
# 
# <pre>pyatk.py shiftoverlay \
#  -vv \
#  -d crop-sa30 \
#  -o png \
#  --ff ../2009.06.10-flats01/*001 --fm 300 \
#  --df ../2009.06.10-darks02/*001 --dm 300 \
#  --safile statoff-2009.06.10-run00/simask-orig-fopt-updated.csv \
#  --sffile sfmask/sfmask-16x16.csv \
#  --shape box \
#  --subap 30 \
#  --skip 0 \
#  --intclip 0.8,1.2 \
#  --shifts subshift-16x16-2009.06.10-run00/image-shifts.npy \
#  --scale 1.5 \
#  ../*_im*</pre>
# 
# Convert to movie:
# <pre>mencoder mf://*.png -mf fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=1000 -o wfwfs_img_shift_movie.avi</pre>