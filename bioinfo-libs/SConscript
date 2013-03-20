import os

Import('debug', 'compiler')

# Initialize environment
env = Environment(CC = compiler, 
                  CFLAGS = '-std=c99 -fopenmp -DCP_HAS_PTHREAD_H -DCP_HAS___BUILTIN_CLZ -DCP_HAS_PTHREAD_H -DCP_HAS___BUILTIN_CLZ -DCP_HAS_STRDUP -DCP_HAS_STRNDUP -DCP_HAS_INET_NTOP -DCP_HAS_SYS_TIME_H -DCP_HAS_GETOPT -DCP_HAS_LONG_LONG -DCP_HAS_DLFCN_H -DCP_DBMS_STATIC ',
                  CPPPATH = [os.getcwd(), ARGUMENTS.get('commons-path', os.getcwd() + '/../common-libs/') ],
                  LIBPATH = ['/usr/lib' ])
env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += ' -O0 -g'
else:
    env['CFLAGS'] += ' -O3'

env['objects'] = []

# Targets
SConscript(['bioformats/SConscript',
            ], exports = ['env'])

SConscript(['aligners/SConscript',
            ], exports = ['env'])

env.Library('bioinfo', env['objects'])

# Should traverse the tree and get the *.c files
#env.SharedLibrary('bioformats', env['objects'])

