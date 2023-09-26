/*
 * Copyright 2009-2015 IPOL Image Processing On Line http://www.ipol.im/
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @mainpage On the Implementation of Collaborative Total Variation
 *           Regularization
 *
 * README.txt:
 * @verbinclude README.txt
 */


/**
 * @file   denoisingPDHG_ipol.cpp
 * @brief  Main executable file
 *
 * @author Joan Duran <joan.duran@uib.es>
 */

#include "libPDHG.h"
#include "io_png.h"

// Usage: denoisingPDHG_ipol noisy.png denoised.png lambda p q r

int main(int argc, char **argv)
{
    if(argc < 7)
	{
		printf("usage: denoisingPDHG_ipol noisy.png denoised.png lambda p q "
               "r\n");
        printf("noisy.png    :: input noisy image.\n");
        printf("denoised.png :: denoised image.\n");
		printf("lambda       :: balancing parameter.\n");
        printf("p q r        :: collaborative norm indices. For the supremum "
               "norm use -1 (for example, use (p,q,r)=(-1,1,1) for the "
               "l^{infty,1,1} norm),\n"
               "                for the (S^1,l^1) norm use (p,q,r)=(1,1,0), "
               " and for the (S^{infty},l^1) norm use (p,q,r)=(-1,1,0).\n");
	    
        return EXIT_FAILURE;
	}
    
    // Read input noisy image
    size_t nx, ny, nc;
    float *d_v = NULL;

    d_v = io_png_read_f32(argv[1], &nx, &ny, &nc);

    if(!d_v) 
	{
	    fprintf(stderr, "Error - %s not found  or not a correct png image.\n",
                argv[1]);
       	return EXIT_FAILURE;
    }

    // Image variables
	int d_w = (int) nx;
	int d_h = (int) ny;
	int d_c = (int) nc;
	int d_wh = d_w * d_h;
	int d_whc = d_c * d_w * d_h;
    
    if(d_c != 3) // We only used color images with three channels
    {
        fprintf(stderr, "Error - %s not a three-color image.\n", argv[1]);
        return EXIT_FAILURE;
   	}

    // Parameters
    float lambda = atof(argv[3]);

	if(lambda < 0.0f)
    {
	    fprintf(stderr, "Error - %s must be nonnegative.\n", argv[3]);
       	return EXIT_FAILURE;
    }
    
    int p = atoi(argv[4]);
    int q = atoi(argv[5]);
    int r = atoi(argv[6]);

    // PDHG parameters
    float tau = 0.5f;
    float sigma = 0.5f;
    float theta = 1.0f;
    float tol = 0.0001f;
    int maxIter = 250;
    
    // TV-denoising using PDHG algorithm
    float **noisy = new float*[d_c];
    float **denoised = new float*[d_c];

    for(int c = 0; c < d_c; c++) 
	{
        noisy[c] = &d_v[c*d_wh];
        denoised[c] = new float[d_wh];
    }
    
    if(PDHG(noisy, denoised, tau, sigma, theta, lambda, p, q, r, tol, maxIter,
            d_c, d_w, d_h) != 1)
	    return EXIT_FAILURE;

    // Save denoised png image
    float *denoised_png = new float[d_whc];
    int k = 0;
	for(int c = 0; c < d_c; c++)
        for(int i = 0; i < d_wh; i++)
        {
            denoised_png[k] = denoised[c][i];
            k++;
        }
    
    if(io_png_write_f32(argv[2], denoised_png, (size_t) d_w, (size_t) d_h,
                        (size_t) d_c) != 0) 
        fprintf(stderr, "Error - Failed to save png image %s.\n", argv[5]);

	// Free memory
    free(d_v);
    delete[] noisy;
    
    for(int c = 0; c < d_c; c++)
        delete[] denoised[c];

	return EXIT_SUCCESS;
}
