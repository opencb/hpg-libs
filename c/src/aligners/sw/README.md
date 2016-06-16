HPG-SW
========

HPG-SW executable to run the Smith-Waterman algorithm using SSE4.2/AVX2 instrucctions.

BUILDING  
Makefile is used as building system to create the **hpg-sw** executable.  
GCC4.4+ or icc compiler and SSE4.2, AVX2 or AVX312 are the only requirements.
To build the executable 'hpg-sw', download the **features/avx512** branch of the **opencb/hpg-libs** project:  
```
  $ git clone https://github.com/opencb/hpg-libs.git
  $ cd hpg-libs/c/src/aligners/sw
  $ git checkout features/avx512
  $ make clean
```
For the SSE4.2 version, compile with:  
```
  $ make SIMD=SSE
```

For the AVX2 version, compile with:  
```
 $ make SIMD=AVX2
```

For the AVX512 version (only available for Intel icc compiler), compile with:  
```
 $ make SIMD=AVX512 COMPILER=intel
```

Use the input parameter TIMING in order to compute processing times of each thread and stage, e.g.:  
```
 $ make SIMD=AVX512 TIMING=1
```

In addition, you can select the one of the two implementations of the Smith-Waterman algorithm by using the paramenter 'SW_VERSION'.
Some examples:
```
 $ make SIMD=AVX2 SW_VERSION=1
 $ make SIMD=AVX512 SW_VERSION=2 COMPILER=intel
```

The executable **hpg-sw** will be stored in the directory 'bin'.


DOCUMENTATION  
You can find more info about HPG project and bioinfo-libs at:  

 http://wiki.opencb.org/projects/hpg/doku.php


CONTACT  
You can contact any of the following developers:  
 * Joaquín Tárraga (jt645@cam.ac.uk)
 * Ignacio Medina (im411@cam.ac.uk)
