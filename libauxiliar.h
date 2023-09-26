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
 * @file libauxiliar.h
 * @brief Auxiliar functions.
 * @author Joan Duran <joan.duran@uib.es>
 */

#ifndef _LIBAUXILIAR_H_
#define _LIBAUXILIAR_H_

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <time.h>
#include <unistd.h>
#include "mt19937ar.h"

#define fTiny 0.00000001f
#define fLarge 100000000.0f

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

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
              int height);

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
                int width, int height);

/**
 * \brief  Compute the sign of a float value.
 *
 * @param[in]  value  input float value.
 * @return sign of value.
 */

int SIGN(float value);

/**
 * \brief  Initialize a float vector. 
 *
 * @param[in]  u  vector input.
 * @param[out] u  vector output.
 * @param[in]  value  value inserted.
 * @param[in]  dim  size of the vector.
 *
 */

void fpClear(float *u, float value, int dim);

/**
 * \brief  Copy the values of a float vector into another.
 *
 * @param[in]  input  vector input.
 * @param[out] output  vector output.
 * @param[in]  dim  size of vectors.
 *
 */

void fpCopy(float *input, float *output, int dim);

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

void fiAddNoise(float *u, float *v, float std, long int randinit, int dim);

#endif
