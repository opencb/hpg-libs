import os

Import('debug', 'formats', 'aligners', 'compiler')

# Initialize environment
vars = Variables('buildvars.py')
vars.Add(PathVariable('CPROPS_INCLUDE_PATH', 'Path to the headers of cprops library', '', PathVariable.PathAccept))
vars.Add(PathVariable('CPROPS_LIBRARY_PATH', 'Path to the compiled cprops library', '', PathVariable.PathAccept))
vars.Add(PathVariable('SAMTOOLS_INCLUDE_PATH', 'Path to the headers of samtools library', '', PathVariable.PathAccept))
vars.Add(PathVariable('SAMTOOLS_LIBRARY_PATH', 'Path to the compiled samtools library', '', PathVariable.PathAccept))

env = Environment(variables = vars,
      	          CC = compiler, 
                  CFLAGS = '-std=c99 -fopenmp',
                  CPPPATH = [os.getcwd(), ARGUMENTS.get('commons-path', os.getcwd() + '/../common-libs/'), '$CPROPS_INCLUDE_PATH', '$SAMTOOLS_INCLUDE_PATH' ],
                  LIBPATH = ['/usr/lib', '$CPROPS_LIBRARY_PATH', '$SAMTOOLS_LIBRARY_PATH' ])
env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += ' -O0 -g'
else:
    env['CFLAGS'] += ' -O3'

env['objects'] = []

# Targets
SConscript(['bioformats/SConscript',
            ], exports = ['env', 'formats'])

SConscript(['aligners/SConscript',
            ], exports = ['env', 'aligners'])

env.Library('bioinfo', env['objects'])

# Should traverse the tree and get the *.c files
#env.SharedLibrary('bioformats', env['objects'])

