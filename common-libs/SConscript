Import('debug', 'compiler')

# Initialize environment
# Current folder is needed in CPPPATH order to compile containers using commons
vars = Variables('buildvars.py')
vars.Add(PathVariable('CPROPS_INCLUDE_PATH', 'Path to the headers of cprops library', '', PathVariable.PathAccept))
vars.Add(PathVariable('CPROPS_LIBRARY_PATH', 'Path to the compiled cprops library', '', PathVariable.PathAccept))

env = Environment(variables = vars,
                  CC = compiler,
		  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_BSD_SOURCE ',
                  CPPPATH = ['#', '/usr/include/libxml2', '.', '$CPROPS_INCLUDE_PATH' ],
                  LIBPATH = ['/usr/lib', '$CPROPS_LIBRARY_PATH' ])

#env = Environment(CFLAGS='-std=c99 -D_XOPEN_SOURCE=600 -D_BSD_SOURCE ',
#                  CPPPATH = ['#', '/usr/include/libxml2', '.' ])
env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3'



# Targets
commons_obj = env.Object(Glob('commons/*.c'))
containers_obj = env.Object(Glob('containers/*.c'))

env.Library('common', commons_obj + containers_obj)
