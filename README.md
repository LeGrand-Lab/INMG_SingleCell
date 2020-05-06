# INMG_SingleCell


# Installing FIt-SNE latest version

## Step 1. Fourier transform
FIt-SNE Requires Fourier trasform algorithm, as indicated here:
http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix

* From Terminal, download the ftp `fftw-3.3.8`, extract, go inside and compile: 

```
$ cd programsSingleCell
$ wget http://www.fftw.org/fftw-3.3.8.tar.gz
$ tar xvzf fftw-3.3.8.tar.gz
$ cd fftw-3.3.8
$ sudo su 
# ./configure
# make
# make install
```
If `make` doesn't work, do `gmake`. Typically requires root privileges (unless you specify a different install directory with the `--prefix` flag to `configure`)

## Step 2. FIt-SNE 

We are going to use R implementation, version 1.2.0, from repo https://github.com/KlugerLab/FIt-SNE/releases
```
$ sudo su
# g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member
```
------
References:

George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger, Yuval Kluger. (2019). Fast interpolation-based t-SNE for improved visualization of single-cell RNA-seq data. Nature Methods. (link)
