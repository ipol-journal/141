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
 * Add zero-mean Gaussian noise to an image.
 *
 * README.txt:
 * @verbinclude README.txt
 */


/**
 * @file   addnoise_ipol.cpp
 * @brief  Main executable file
 *
 * @author Joan Duran <joan.duran@uib.es>
 */

#include "libauxiliar.h"
#include "io_png.h"

// Usage: addnoise_ipol input.png noisy.png std

int main(int argc, char **argv)
{
    if(argc < 4)
	{
		printf("usage: addnoise_ipol option input.png noisy.png sigma\n");
        printf("input.png  :: input noise-free image.\n");
        printf("noisy.png  :: noisy image.\n");
        printf("std        :: noise standard deviation.\n");
	    
        return EXIT_FAILURE;
	}
    
    // Read input image
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

    // Noise standard deviation
    float std = atof(argv[3]);

    if(std < 0.0f)
    {
	    fprintf(stderr, "Error - %s must be nonnegative.\n", argv[2]);
      	return EXIT_FAILURE;
   	}

    // Add noise to input image
   	float *noisy = new float[d_whc];
        
    for(int i = 0; i < d_c; i++)
        fiAddNoise(&d_v[i*d_wh], &noisy[i*d_wh], std, i, d_wh);

    // Save noisy png image
	if(io_png_write_f32(argv[2], noisy, (size_t) d_w, (size_t) d_h,
                        (size_t) d_c) != 0)
		    fprintf(stderr, "Error - Failed to save png image %s.\n", argv[2]);
    
	// Free memory
    free(d_v);
    delete[] noisy;

	return EXIT_SUCCESS;
}
