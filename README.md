# INMG_SingleCell

!! under construction !!

## Technical notes

**UMAP vs FIt-SNE**

Among cluster visualization algorithms, UMAP and FIt-SNE have demonstrated their superiority, with an slight advantage of UMAP in terms of time and consistence. However some authors highlight the capacity of tSNE to distinguish clusters embedded into bigger clusters (case of "containment")<https://pair-code.github.io/understanding-umap/>, so the relatively recent FIt-SNE which largely outperforms traditional tSNE is a very good choice, not inferior to UMAP. A new version of FIt-SNE (v.1.2) exists but as no r wrapper is available to date,  previous version (v1.1.0) has been used for homeostasis (D0) analysis of published papers in this repo.

#### Installing FIt-SNE 

#### Step 1. Fourier transform
FIt-SNE Requires Fourier trasform algorithm, as indicated here:
http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix. For MacOS and Linux the installation follows exactly the same steps.

* From Terminal, download the ftp `fftw-3.3.8`, extract, go inside and compile: 

```
$ cd programsSC
$ wget http://www.fftw.org/fftw-3.3.8.tar.gz
$ tar xvzf fftw-3.3.8.tar.gz
$ cd fftw-3.3.8
$ sudo su 
# ./configure
# make
# make install
```
If `make` doesn't work, do `gmake`. Typically requires root privileges (unless you specify a different install directory with the `--prefix` flag to `configure`)

#### Step 2. FIt-SNE 

We are going to use R implementation, from repo https://github.com/KlugerLab/FIt-SNE

```
$ cd programsSC
$ git clone https://github.com/KlugerLab/FIt-SNE.git
$ cd FIt-SNE/
$ sudo su
# g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member
```
Important information about latest releases (version 1.2.0, version 1.2.1): https://github.com/KlugerLab/FIt-SNE/releases

Quoting Dr.Linderman (team KlugerLab):
"""
Several changes to default FIt-SNE settings to make it more suitable for embedding large datasets.

Major changes to default values:
-Learning rate increased from the fixed value of 200 to max(200, N/early_exag_coeff).
-Iteration number decreased from 1000 to 750.
-Initialization is set to PCA (computed via fast SVD implementations like ARPACK).

Minor changes:
-Late exaggeration start is set to the end of early exaggeration (if late exaggeration coefficient is provided).
-Limiting max step size to 5 (solves problem encountered when learning rate set too high and attractive forces cause a small number of points to overshoot)

""" 

----------
References:

George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger, Yuval Kluger. (2019). Fast interpolation-based t-SNE for improved visualization of single-cell RNA-seq data. Nature Methods. 

Kobak, D. & Berens, P. The art of using t-SNE for single-cell transcriptomics. Nat. Commun. 10, (2019).
----------------
----------------
Johanna GL

