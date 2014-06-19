# This SConstruct launches all Sconscript files inside the library directories 
# The Environment() is created in the SConstruct script
# This dir can be built standalone by executing scons here, or together
# by executing scons in the parent directory
import sys
import os

#Compiler configure
debug = int(ARGUMENTS.get('debug', '0'))
compiler = ARGUMENTS.get('compiler', 'gcc')

#Paths
system_include = '/usr/include'
system_libpath = '/usr/lib'

commons_path = '#/common-libs'
bioinfo_path = '#/bioinfo-libs'

build_tools = ['default']

#Build environment
env = Environment(tools = build_tools,
		  CFLAGS = ' -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE -D_GNU_SOURCE -D_REENTRANT ',
		  CPPPATH = ['.', '#', system_include, '%s/libxml2' % system_include, '%s' % commons_path, '%s' % bioinfo_path], 
		  LIBPATH = [system_libpath],
		  LINKFLAGS = [],
		  LIBS = ['xml2', 'm', 'z', 'curl'])

if os.environ.has_key('CPATH'):
    for dir in os.getenv('CPATH').split(':'):
        env.Append(CPPPATH=[dir])

if os.environ.has_key('LIBRARY_PATH'):
    for dir in os.getenv('LIBRARY_PATH').split(':'):
        env.Append(LIBPATH=[dir])

if compiler == 'icc':
	env['tools'] += ['intelc']
	env['CFLAGS'] += ' -msse4.2 -openmp '
	env['LIBS'] += ['irc']
	env['LINKFLAGS'] += ['-openmp']
else:
	env['CFLAGS'] += ' -fopenmp '
	env['LINKFLAGS'] += ['-fopenmp']

env['objects'] = []
env.Decider('MD5-timestamp')

SConscript('common-libs/SConscript', exports = ['env', 'debug', 'compiler'])
SConscript('bioinfo-libs/SConscript', exports = ['env', 'debug', 'compiler'])
SConscript('math/SConscript', exports = ['debug', 'compiler'])
