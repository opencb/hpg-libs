Import('debug', 'compiler')

# Initialize environment
env = Environment(CC = compiler,
		  CFLAGS = '-std=c99 -D_GNU_SOURCE -D_XOPEN_SOURCE=600 -DCP_HAS_PTHREAD_H -DCP_HAS___BUILTIN_CLZ -DCP_HAS_PTHREAD_H -DCP_HAS___BUILTIN_CLZ -DCP_HAS_STRDUP -DCP_HAS_STRNDUP -DCP_HAS_INET_NTOP -DCP_HAS_SYS_TIME_H -DCP_HAS_GETOPT -DCP_HAS_LONG_LONG -DCP_HAS_DLFCN_H -DCP_DBMS_STATIC ',
                  CPPPATH = ['#', '/usr/include/libxml2' ],
                  LIBPATH = ['/usr/lib' ])

env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3 -g'



# Targets
commons_obj = env.Object(Glob('commons/*.c'))
containers_obj = env.Object(Glob('containers/*.c'))
cprops_obj = env.Object(['containers/cprops/avl.c', 'containers/cprops/hashtable.c', 'containers/cprops/linked_list.c', 'containers/cprops/vector.c', 'containers/cprops/heap.c'])


env.Library('common', commons_obj + containers_obj + cprops_obj)
