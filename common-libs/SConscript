Import('debug')

env = Environment(CFLAGS='-std=c99 -D_XOPEN_SOURCE=600 -D_BSD_SOURCE ')
env.Decider('MD5-timestamp')


# Debug/release mode
if debug == 1:
    env['CFLAGS'] += '-O0 -g'
else:
    env['CFLAGS'] += '-O3'


# Check dependency libraries: XML2, cURL
conf = Configure(env)

if conf.CheckLib('xml2'):
    # Get include and libs arguments
    env.ParseConfig("pkg-config libxml-2.0 --cflags --libs")
else:
    print 'XML library not found!'
    # TODO: Download XML2 and set its -I, -L, -l flags
    Exit(1)

if conf.CheckLib('curl'):
    # Get include and libs arguments
    env.ParseConfig("pkg-config libcurl --cflags --libs")
else:
    print 'cURL library not found!'
    # TODO: Download cURL and set its -I, -L, -l flags
    Exit(1)

env = conf.Finish()


# Set extra include and libs flags
# Current folder is needed in order to compile containers using commons
env['CPPPATH'] += ['.']


# Targets
commons_obj = env.Object(Glob('commons/*.c'))
containers_obj = env.Object(Glob('containers/*.c'))

env.Library('common', commons_obj + containers_obj)

#env.SharedLibrary('common', Glob('commons/*.c') + Glob('containers/*.c'))

# Also compile tests automatically, can be avoided by expliciting a target
#test = env.Program(target = 'commons/test/http-utils.test', 
                    #source = ['commons/test/test_http_utils.c', 'commons/http_utils.c'], 
                    #LIBS = ['check', 'curl'])
