import os

Import('debug')

env = Environment(CFLAGS = '-std=c99 ')
env.Decider('MD5-timestamp')

# Debug/release mode
if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3'


env['CPPPATH'] = [os.getcwd(), 
                  ARGUMENTS.get('commons-path', os.getcwd() + '/../common-libs/')]
env['objects'] = []

# Targets
SConscript(['bioformats/SConscript',
            #'aligners/SConscript',
            ], exports = 'env')

env.Library('bioinfo', env['objects'])

# Should traverse the tree and get the *.c files
#env.SharedLibrary('bioformats', env['objects'])

