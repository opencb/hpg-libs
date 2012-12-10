import os

my_home = os.environ.get('HOME', '')

# Build variables
CPROPS_INCLUDE_PATH = my_home + '/ext/libcprops-0.1.12/include'
CPROPS_LIBRARY_PATH = my_home + '/ext/libcprops-0.1.12/lib'
SAMTOOLS_INCLUDE_PATH = my_home + '/ext/samtools-0.1.18/'
SAMTOOLS_LIBRARY_PATH = my_home + '/ext/samtools-0.1.18/'
