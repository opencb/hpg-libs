import os

my_home = os.environ.get('HOME', '')

# Build variables
CPROPS_INCLUDE_PATH = my_home + '/ext/libcprops-0.1.12/include'
CPROPS_LIBRARY_PATH = my_home + '/ext/libcprops-0.1.12/lib'
