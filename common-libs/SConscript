Import('debug', 'compiler')

# Initialize environment
env = Environment(CC = compiler,
		  CFLAGS = '-std=c99 -D_GNU_SOURCE -D_XOPEN_SOURCE=600 ',
                  CPPPATH = ['#', '/usr/include/libxml2', '.' ],
                  LIBPATH = ['/usr/lib' ])

env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3'



# Targets
commons_obj = env.Object(Glob('commons/*.c'))
containers_obj = env.Object(Glob('containers/*.c'))


# Initialize environment
cpropsenv = Environment(CC = compiler,
	  		CFLAGS = '-D_REENTRANT -D_GNU_SOURCE -DHAVE_CONFIG_H ',
            		CPPPATH = ['#', '/usr/include/libxml2', '.' ],
			LIBPATH = ['/usr/lib' ])

cpropsenv.Decider('MD5-timestamp')

if debug == 1:
    cpropsenv['CFLAGS'] += '-O0 -g'
else:
    cpropsenv['CFLAGS'] += '-O3'

# Compile cprops objects but ONLY those used for our libraries
cprops_obj = cpropsenv.Object(['containers/cprops/avl.c', 'containers/cprops/hashtable.c', 'containers/cprops/linked_list.c', 'containers/cprops/vector.c', 'containers/cprops/heap.c', 'containers/cprops/mempool.c', 'containers/cprops/rb.c'])


# Objects
env.Library('common', commons_obj + containers_obj + cprops_obj)
