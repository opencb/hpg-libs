# This SConstruct launches all Sconscript files inside the library directories 
# The Environment() is created in the SConstruct script
# This dir can be built standalone by executing scons here, or together
# by executing scons in the parent directory
import sys

debug = int(ARGUMENTS.get('debug', '0'))

compiler = ARGUMENTS.get('compiler', 'gcc')

SConscript('common-libs/SConscript', exports = ['debug', 'compiler'])
SConscript('bioinfo-libs/SConscript', exports = ['debug', 'compiler'])
SConscript('math/SConscript', exports = ['debug', 'compiler'])

