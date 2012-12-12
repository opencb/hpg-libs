
# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/bioinfo-libs'
commons_path = '#lib/common-libs'
#math_path = '#libs/math'

vars = Variables('buildvars.py')
#vars.Add(PathVariable('ARGTABLE_INCLUDE_PATH', 'Path to the headers of argtable2 library', '', PathVariable.PathAccept))
#vars.Add(PathVariable('ARGTABLE_LIBRARY_PATH', 'Path to the compiled argtable2 library', '', PathVariable.PathAccept))

vars.Add(PathVariable('CPROPS_INCLUDE_PATH', 'Path to the headers of cprops library', '', PathVariable.PathAccept))
vars.Add(PathVariable('CPROPS_LIBRARY_PATH', 'Path to the compiled cprops library', '', PathVariable.PathAccept))

vars.Add(PathVariable('SAMTOOLS_INCLUDE_PATH', 'Path to the headers of samtools library', '', PathVariable.PathAccept))
vars.Add(PathVariable('SAMTOOLS_LIBRARY_PATH', 'Path to the compiled samtools library', '', PathVariable.PathAccept))

#vars.Add(PathVariable('CONFIG_INCLUDE_PATH', 'Path to the headers of config library', '', PathVariable.PathAccept))
#vars.Add(PathVariable('CONFIG_LIBRARY_PATH', 'Path to the compiled config library', '', PathVariable.PathAccept))

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		  CC = compiler,
                  variables = vars,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -g',
                  CPPPATH = ['#', '#src', '#include', '$SAMTOOLS_INCLUDE_PATH', '$CPROPS_INCLUDE_PATH', bioinfo_path, commons_path ],
#                  CPPPATH = ['#', '#src', '#include', '/usr/local/include', '/usr/include/libxml2', '$CPROPS_INCLUDE_PATH', '$ARGTABLE_INCLUDE_PATH', bioinfo_path, commons_path, math_path ],
                  LIBPATH = ['#libs/common-libs/', '$SAMTOOLS_LIBRARY_PATH', '$CPROPS_LIBRARY_PATH', commons_path ],
                  LIBS = ['argtable2', 'common', 'config', 'bam', 'cprops', 'm', 'z'],
                  LINKFLAGS = ['-fopenmp'])
                  
if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

env['objects'] = []


# Targets
formats = ['fastq', 'bam-sam']
aligners = ['bwt', 'sw']

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'formats', 'aligners', 'compiler'])

env.Program('hpg-aligner',
             source = [Glob('src/*.c'),
                       "%s/libcommon.a" % commons_path,
                       "%s/libbioinfo.a" % bioinfo_path
                      ]
           )

env.Install('#bin', 'hpg-aligner')

'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
