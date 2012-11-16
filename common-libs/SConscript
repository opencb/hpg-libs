Import('debug')

# Initialize environment
# Current folder is needed in CPPPATH order to compile containers using commons
env = Environment(CFLAGS='-std=c99 -D_XOPEN_SOURCE=600 -D_BSD_SOURCE ',
                  CPPPATH = ['#', '/usr/include/libxml2', '.' ])
env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3'



# Targets
commons_obj = env.Object(Glob('commons/*.c'))
containers_obj = env.Object(Glob('containers/*.c'))

env.Library('common', commons_obj + containers_obj)
