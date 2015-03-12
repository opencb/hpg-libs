# This SConstruct launches all Sconscript files inside the library directories 

import sys
import os

# Compiler configure
debug = int(ARGUMENTS.get('debug', '0'))
compiler = ARGUMENTS.get('compiler', 'gcc')

# Paths
system_include = '/usr/include'
system_libpath = '/usr/lib'
third_party_path = os.getcwd() + '/third_party'

build_tools = ['default']
if compiler == 'intel':
    build_tools += ['intelc']
elif compiler == 'clang':
    build_tools += ['clang','clang++']


# Build C environment
hpg_c_env = Environment(TOOLS = build_tools,
                        CFLAGS = ' -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE -D_GNU_SOURCE -D_REENTRANT ',
                        CPPPATH = ['.', '#', system_include, '%s/libxml2' % system_include, '%s' % third_party_path], 
                        LIBPATH = [system_libpath],
                        LINKFLAGS = [],
                        LIBS = ['xml2', 'm', 'z', 'curl'])


if os.environ.has_key('CPATH'):
    for dir in os.getenv('CPATH').split(':'):
        hpg_c_env.Append(CPPPATH=[dir])

if os.environ.has_key('LIBRARY_PATH'):
    for dir in os.getenv('LIBRARY_PATH').split(':'):
        hpg_c_env.Append(LIBPATH=[dir])

if compiler == 'intel':
    hpg_c_env['CFLAGS'] += ' -msse4.2 -openmp '
    hpg_c_env['LIBS'] += ['irc']
    hpg_c_env['LINKFLAGS'] += ['-openmp']
else:
    hpg_c_env['CFLAGS'] += ' -fopenmp '
    hpg_c_env['LINKFLAGS'] += ['-fopenmp']

if debug == 1:
    hpg_c_env['CFLAGS'] += ' -O0 -g'
else:
    hpg_c_env['CFLAGS'] += ' -O2 '

hpg_c_env['objects'] = []
hpg_c_env.Decider('MD5-timestamp')




# Build C++ environment
hpg_cpp_env = Environment(TOOLS = build_tools,
                          CCFLAGS = ' -Wall -std=c++11 ',
                          CPPPATH = ['.', '#', '#/cpp/src', system_include, '%s/libxml2' % system_include, '%s' % third_party_path], 
                          LIBPATH = [system_libpath],
                          LINKFLAGS = [],
                          LIBS = ['xml2', 'm', 'z', 'curl', 'boost_regex'])

if os.environ.has_key('CPPPATH'):
    for dir in os.getenv('CPPPATH').split(':'):
        hpg_cpp_env.Append(CPPPATH=[dir])

if os.environ.has_key('LIBRARY_PATH'):
    for dir in os.getenv('LIBRARY_PATH').split(':'):
        hpg_cpp_env.Append(LIBPATH=[dir])

if compiler == 'intel':
    hpg_cpp_env['CCFLAGS'] += ' -msse4.2 -openmp '
    hpg_cpp_env['LIBS'] += ['irc']
    hpg_cpp_env['LINKFLAGS'] += ['-openmp']
else:
    hpg_cpp_env['CCFLAGS'] += ' -fopenmp '
    hpg_cpp_env['LINKFLAGS'] += ['-fopenmp']

if debug == 1:
    hpg_cpp_env['CCFLAGS'] += ' -O0 -g'
else:
    hpg_cpp_env['CCFLAGS'] += ' -O3 '

hpg_cpp_env['objects'] = []
hpg_cpp_env.Decider('MD5-timestamp')




# Third party
SConscript('third_party/SConscript', exports = ['hpg_c_env', 'hpg_cpp_env', 'debug', 'compiler'])


# Our source code
SConscript('c/SConscript', exports = ['hpg_c_env', 'debug', 'compiler'])
SConscript('cpp/SConscript', exports = ['hpg_cpp_env', 'debug', 'compiler'])

