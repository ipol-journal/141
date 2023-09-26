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
 * @file libauxiliar.cpp
 * @brief Auxiliar functions.
 * @author Joan Duran <joan.duran@uib.es>
 */

#include "libauxiliar.h"

/**
 * \brief  Compute discrete gradient operator via forward differences.
 *
 * @param[in]  u  input vector : the first pointer accounts for the channel and
 *             the second one for the pixel position.
 * @param[out] ux, uy  horizontal and vertical forward differences : the first
 *             pointer accounts for the channel and the second one for the
 *             pixel position.
 * @param[in]  num_channels  number of channels of the image.
 * @param[in]  width, height  image size.
 *
 */

void gradient(float **u, float **ux, float **uy, int num_channels, int width,
              int height)
{
    // Compute discrete gradient using forward differences
    for(int k = 0; k < num_channels; k++)
    {
#pragma omp parallel for
        for(int j = 0; j < height; j++)
        {
            int l = j * width;
            
            for(int i = 0; i < width; i++)
            {
                // Derivatives in the x-direction
                if(i != width-1)
                    ux[k][i+l] = u[k][i+1+l] - u[k][i+l];
                else
                    ux[k][i+l] = 0.0f;
                
                // Derivatives in the y-direction
                if(j != height-1)
                    uy[k][i+l] = u[k][i+width+l] - u[k][i+l];
                else
                    uy[k][i+l] = 0.0f;
            }
        }
    }
}

/**
 * \brief  Compute divergence operator as @f$ \langle -\mbox{div} p, u \rangle =
 *         \langle p, \nabla  u\rangle @f$.
 *
 * @param[in]  px, py  dual variables : the first pointer accounts for the
 *             channel and the second one for the pixel position.
 * @param[out] div  divergence operator : the first pointer accounts for the
 *             channel and the second one for the pixel position.
 * @param[in]  num_channels  number of channels of the image.
 * @param[in]  width, height  image size.
 *
 */

void divergence(float **px, float **py, float **div, int num_channels,
                int width, int height)
{
    // Image size
    int dim = width * height;
    
    // Compute divergence as <-divp, u> = <p, grad u>
#pragma omp parallel for
    for(int k = 0; k < num_channels; k++)
    {
        // Fill the output vector with zeros
        fpClear(div[k], 0.0f, dim);
        
        for(int j = 0; j < height; j++)
        {
            int l = j * width;
            
            for(int i = 0; i < width; i++)
            {
                if(i != width-1)
                {
                    // ux[k][i+l] = u[k][i+1+l] - u[k][i+l]
                    div[k][i+1+l] -= px[k][i+l];
                    div[k][i+l] += px[k][i+l];
                }
                
                if(j != height-1)
                {
                    // uy[k][i+l] = u[k][i+width+l] - u[k][i+l]
                    div[k][i+width+l] -= py[k][i+l];
                    div[k][i+l] += py[k][i+l];
                }
            }
        }
    }
}

/**
 * \brief  Compute the sign of a float value.
 *
 * @param[in]  value  input float value.
 * @return sign of value.
 */

int SIGN(float value)
{
    int sgn;
    
    if(value > 0.0f)  sgn = 1;
    else if(value < 0.0f)  sgn = -1;
    else  sgn = 0;
    
    return sgn;
}

/**
 * \brief  Initialize a float vector.
 *
 * @param[in]  u  vector input.
 * @param[out] u  vector output.
 * @param[in]  value  value inserted.
 * @param[in]  dim  size of the vector.
 *
 */

void fpClear(float *u, float value, int dim)
{
    for(int i = 0; i < dim; i++)
        u[i] = value;
}

/**
 * \brief  Copy the values of a float vector into another.
 *
 * @param[in]  input  vector input.
 * @param[out] output  vector output.
 * @param[in]  dim  size of vectors.
 *
 */

void fpCopy(float *input, float *output, int dim)
{
    if(input != output)
        memcpy((void *) output, (const void *) input, dim * sizeof(float));
}

/**
 * \brief  Add white Gaussian noise to an image.
 *
 * @param[in]  u  original image.
 * @param[out] v  noised image.
 * @param[in]  std  noise standard deviation.
 * @param[in]  randinit  random parameter.
 * @param[in]  dim  image size.
 *
 */

void fiAddNoise(float *u, float *v, float std, long int randinit, int dim)
{
    mt_init_genrand((unsigned long int) time (NULL) +
                    (unsigned long int) getpid() +
                    (unsigned long int) randinit);
    
    for(int i = 0; i < dim; i++)
    {
        double a = mt_genrand_res53();
        double b = mt_genrand_res53();
        double z = (double)(std) * sqrt(-2.0 * log(a)) * cos(2.0 * M_PI * b);
        
        v[i] =  u[i] + (float) z;
    }
}
