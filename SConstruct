# This SConstruct launches all Sconscript files inside the library directories 

import sys
import os

#Compiler configure
debug = int(ARGUMENTS.get('debug', '0'))
compiler = ARGUMENTS.get('compiler', 'gcc')

#Paths
system_include = '/usr/include'
system_libpath = '/usr/lib'
third_party_path = '#/third_party'

build_tools = ['default']
if compiler == 'icc':
    build_tools += ['intelc']


#Build environment
hpg_env = Environment(TOOLS = build_tools,
		  CFLAGS = ' -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE -D_GNU_SOURCE -D_REENTRANT ',
		  CPPPATH = ['.', '#', system_include, '%s/libxml2' % system_include, '%s' % third_party_path], 
		  LIBPATH = [system_libpath],
		  LINKFLAGS = [],
		  LIBS = ['xml2', 'm', 'z', 'curl'])


if os.environ.has_key('CPATH'):
    for dir in os.getenv('CPATH').split(':'):
        hpg_env.Append(CPPPATH=[dir])

if os.environ.has_key('LIBRARY_PATH'):
    for dir in os.getenv('LIBRARY_PATH').split(':'):
        hpg_env.Append(LIBPATH=[dir])

if compiler == 'icc':
	hpg_env['CFLAGS'] += ' -msse4.2 -openmp '
	hpg_env['LIBS'] += ['irc']
	hpg_env['LINKFLAGS'] += ['-openmp']
else:
	hpg_env['CFLAGS'] += ' -fopenmp '
	hpg_env['LINKFLAGS'] += ['-fopenmp']

hpg_env['objects'] = []
hpg_env.Decider('MD5-timestamp')


# Third party
third_party_env = Environment(TOOLS = hpg_env['TOOLS'])
SConscript('third_party/SConscript', exports = ['third_party_env', 'debug', 'compiler'])

hpg_env['CPPPATH'] += third_party_env['CPPPATH']


# our src

SConscript('c/SConscript', exports = ['hpg_env', 'third_party_env', 'debug', 'compiler'])
SConscript('cpp/SConscript', exports = ['hpg_env', 'third_party_env', 'debug', 'compiler'])
