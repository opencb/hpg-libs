import os

Import('debug')

# Initialize environment
vars = Variables('buildvars.py')
vars.Add(PathVariable('CPROPS_INCLUDE_PATH', 'Path to the headers of cprops library', '', PathVariable.PathAccept))
vars.Add(PathVariable('CPROPS_LIBRARY_PATH', 'Path to the compiled cprops library', '', PathVariable.PathAccept))

env = Environment(variables = vars,
                  CFLAGS = '-std=c99 ',
                  CPPPATH = [os.getcwd(), ARGUMENTS.get('commons-path', os.getcwd() + '/../common-libs/'), '$CPROPS_INCLUDE_PATH' ],
                  LIBPATH = ['/usr/lib', '$CPROPS_LIBRARY_PATH' ])
env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3'

env['objects'] = []



# Targets
SConscript(['bioformats/SConscript',
            #'aligners/SConscript',
            ], exports = 'env')

env.Library('bioinfo', env['objects'])

# Should traverse the tree and get the *.c files
#env.SharedLibrary('bioformats', env['objects'])

