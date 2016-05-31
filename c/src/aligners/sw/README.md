sw
========

HPG-SW executable for the Smith-Waterman algorithm using SSE4.2/AVX2 instrucctions.

BUILDING
Makefile is used as building system to create the 'hpg-sw' executable.
GCC4.4+ are SSE4.2 or AVX2 are the only requirements.
To build the executable 'hpg-sw', download the 'features/avx2' branch of the 'hpg-libs' project:

 $ git clone https://github.com/opencb/hpg-libs.git
 $ git checkout features/avx2
 $ cd c/src/aligners/sw
 $ make clean

For the SSE4.2 version, compile with:
 $ make SIMD=SSE

For the AVX2 version, compile with:
 $ make SIMD=AVX2

The executable 'hpg-sw' will be stored in the directory 'bin'.


DOCUMENTATION
You can find more info about HPG project and bioinfo-libs at:

 http://wiki.opencb.org/projects/hpg/doku.php


CONTACT
You can contact any of the following developers:
 * Joaquín Tárraga (joaquintarraga@gmail.com)
 * Ignacio Medina (igmeca@cipf.es)
