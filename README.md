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

We are going to use R implementation, from repo https://github.com/KlugerLab/FIt-SNE

```
$ cd programsSingleCell
$ git clone https://github.com/KlugerLab/FIt-SNE.git
$ cd FIt-SNE/
$ sudo su
# g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member
```
Important information about latest releases (version 1.2.0, version 1.2.1): https://github.com/KlugerLab/FIt-SNE/releases

"""
Several changes to default FIt-SNE settings to make it more suitable for embedding large datasets.

Major changes to default values:
-Learning rate increased from the fixed value of 200 to max(200, N/early_exag_coeff).
-Iteration number decreased from 1000 to 750.
-Initialization is set to PCA (computed via fast SVD implementations like ARPACK).

Minor changes:
-Late exaggeration start is set to the end of early exaggeration (if late exaggeration coefficient is provided).
-Limiting max step size to 5 (solves problem encountered when learning rate set too high and attractive forces cause a small number of points to overshoot)

""" *Dr. Linderman*

----------
References:

George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger, Yuval Kluger. (2019). Fast interpolation-based t-SNE for improved visualization of single-cell RNA-seq data. Nature Methods. 

Kobak, D. & Berens, P. The art of using t-SNE for single-cell transcriptomics. Nat. Commun. 10, (2019).


