On the Implementation of Collaborative Total Variation Regularization

Joan Duran, joan.duran@uib.es, University of Balearic Islands (Spain)
Michael Moeller, michael.moeller@in.tum.de, Technical University of Munich (Germany)
Catalina Sbert, catalina.sbert@uib.es, University of Balearic Islands (Spain)
Daniel Cremers, cremers@tum.de, Technical University of Munich (Germany)

# OVERVIEW

This C source code accompanies with Image Processing On Line (IPOL) article
“On the Implementation of Collaborative Total Variation Regularization“ at 

    http://www.ipol.im/pub/pre/

This code is used by the online IPOL demo:

    http://dev.ipol.im/

This program reads and writes PNG images, but can be easily adapted to any
other file format. Only 8bit PNG images are handled. 

Three programs are provided:

* ‘denoisingPDHG_ipol' reads a noisy png image and then denoises it using the
Primal-Dual Hybrid Gradient Method for vectorial TV regularization with collaborative norms.

* ‘addnoise_ipol’ adds zero-mean Gaussian noise to a noisy-free png image.

* 'imdiff_ipol' visualizes the difference between two images in such a way that
the error range is linearly transformed from [-4*sigma, 4*sigma] to [0, 255].
Errors outside this range are saturated to 0 and 255, respectively. It also
computes the Root Mean Squared Error and the Peak Signal-to-Noise Ratio:    
    - RMSE = (1/N sum |A[i] - B[i]|^2)^1/2,
    - PNSR = 10 log(255^2 / RMSE^2).


# USAGE

Usage: denoisingPDHG_ipol noisy.png denoised.png lambda p q r

noisy.png    :: input noisy image.
denoised.png :: denoised image.
lambda       :: balancing parameter.
p q r        :: collaborative norm indices. For the supremum norm use -1 (for
                example, use (p,q,r)=(-1,1,1) for the \ell^{\infty,1,1} norm),
                for the (S^1,\ell^1) norm use (p,q,r)=(1,1,0), and for the
                (S^{\infty},\ell^1) norm use (p,q,r)=(-1,1,0).

Usage: addnoise_ipol input.png noisy.png std

input.png  :: input noise-free image.
noisy.png  :: noisy image.
std        :: noise standard deviation.

Usage: imdiff_ipol image1.png image2.png imdiff.png sigma 

image1.png : first image.
image2.png : second image.
imdiff.png : difference image.
sigma      : noise standard deviation.

This program also provides on screen the RMSE and the PSNR values.


#LICENSE

Files mt199937.ar.c and mt19937.ar.h are copyright Makoto Matsumoto and 
Takuji Nishimura. Files io_png.c and io_png.h are copyright Nicolas Limare.
These files are distributed under the BSD license conditions described in
the corresponding headers files.

All the other files are distributed under the terms of the GNU General
Public License Version 3 (the "GPL").


# REQUIREMENTS

The code is written in ANSI C and C++, and should compile on any system 
with an ANSI C/C++ compiler.

The libpng header and libraries are required on the system for compilation
and execution. 

The implementation uses OPENMP which not supported by old versions of
gcc (older than the 4.2). 


# COMPILATION

Simply use the provided makefile, with the command 'make'.


# EXAMPLE
Add noise to the input image and denoises it using the \ell^{\infty,1,1} norm
  
    ./addnoise_ipol traffic.png noisy.png 25
    ./denoisingPDHG_ipol noisy.png denoised.png 0.05 -1 1 1
    ./imdiff_ipol traffic.png denoised.png difference.png 25
