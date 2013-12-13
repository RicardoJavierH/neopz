#! /usr/bin/env python2.7
#***************************************************************************
#*   Copyright (C) 2013 by Edson Borin                                     *
#*   edson@ic.unicamp.br                                                   *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 2 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#**************************************************************************/

# ---------------------------------------------
# Performance test module

# List of rdt files generated by the test
rdtfiles_l = [
	# short_name, option, filename, description
	("sor", "-sor_rdt", "skyline.porder_7.sor.rdt", "SolveSOR: matrix->SolveSOR(...). Solve a linear system using successive over-relaxation method.")
]

def short_description() : return "skyline-solveSOR -- cube1.txt -- polinomial order 7"

def long_description():
	desc =  "Executes the skyline tool collecting statistics for the following steps:"
	for rdtarg in rdtfiles_l :
		desc = desc + '\n\t' + rdtarg[0] + ' (' + rdtarg[1] + ' ' + rdtarg [2] +  ') : ' + rdtarg[3]
	return desc
# ---------------------------------------------

import sys
import os.path
import shlex, subprocess
import resource

# Try to import rdt and stats modules, if available.
import sys

# Variables to be defined by cmake
builddir="@PERFTEST_BASE_DIR@"
datadir="@PERFTEST_SMALL_DATA_DIR@"

def error(message, status):
	sys.stderr.write('ERROR (test.py): '+message+'\n')
        sys.exit(status)

# Setup the command line
def setup_cmd():
	# Check build directory
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 5)
	# Check run directory
	rundir = os.path.join(builddir,'scripts','skyline_tst12')
	if not os.path.isdir(rundir) :
		error(rundir+' is an invalid run directory.', 1)
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 1)
	# Check executable
	executable=os.path.join(builddir,"progs","skyline", "skyline-perf")
	if not os.path.isfile(executable) :
		error(executable+' is an invalid executable file name.', 1)
	# Check input file
	inputfn = os.path.join(datadir,"substruct","inputs","cube1.txt")
	if not os.path.isfile(inputfn) :
		error(inputfn+' is an invalid input file name.', 1)	
	# Put the arguments together
	arguments = ' -if ' + inputfn
        arguments = arguments + ' -porder 7'
	for rdtarg in rdtfiles_l :
		arguments = arguments + ' ' + rdtarg[1] + ' ' + rdtarg[2]
	# TODO: Add arguments to enforce output checking!
	return rundir, executable+arguments

# Limits for this test
# TODO: change cpu limit acording to program execution time
limits = { "cpu"   : (resource.RLIMIT_CPU, 3600, "Max CPU time in seconds"), 
#	   "nofile": (resource.RLIMIT_NOFILE,     7, "The maximum number of open file descriptors for the current process."),
#	   "rss"   : (resource.RLIMIT_RSS,     1024, "The maximum resident set size that should be made available to the process"),
#	   "fsize" : (resource.RLIMIT_FSIZE,      1, "Max size of a file which the process may create"),
#	   "data"  : (resource.RLIMIT_DATA,    1024, "The maximum size (in bytes) of the process's heap"),
#	   "nproc" : (resource.RLIMIT_NPROC,      0, "The maximum number of processes the current process may create")
	 }

# Set the rlimits of the chidren process (see limits above)
# TODO: Improve the handling of sandboxing limits
def setlimits():
	print "Setting resource limit in child"
	for k, v in limits.iteritems() : 
		resource.setrlimit(v[0], (v[1],v[1])) 
		#print k, " : ", v[0], " => ", v[1]

# Sumarizes the RDT (Raw data table) files information
def sumarize_rdt_files(rundir) :
	results = {}
	for f in rdtfiles_l : 
		rdt_id  = f[0]   # Step name
		rdt_fn  = os.path.join(rundir,f[2]) # RDT file name
		rdt_dsc = f[3]   # Description
		results[rdt_id] = (rdt_fn, rdt_dsc)
	return results

# Execute the test.
def run_test(ntimes):
	rundir,cmd=setup_cmd()
        print short_description() 
	print "CMD:",cmd
	args = shlex.split(cmd)
	sout = None
	serr = None
	for i in range(ntimes) : 
		p = subprocess.Popen(args, preexec_fn=setlimits, stdout=sout, stderr=serr, cwd=rundir)
		p.wait()
		if (p.returncode != 0) : 
			return p.returncode, {}
	results = sumarize_rdt_files(rundir)
	return 0, results

