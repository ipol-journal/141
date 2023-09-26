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
 * @file libProximal.h
 * @brief Implementation of proximal mappings.
 * @author Joan Duran <joan.duran@uib.es>
 */

#ifndef _LIBPROXIMAL_H_
#define _LIBPROXIMAL_H_

#include "libauxiliar.h"

#define INFNORM -1

/**
 * \brief  Compute the solution to the proximal mapping of the fidelity term.
 *
 * @param[out] u  solution to the proximal problem : the first pointer accounts
 *             for the channel and the second one for the pixel position.
 * @param[in]  v  argument of the proximal problem : the first pointer accounts
 *             for the channel and the second one for the pixel position.
 * @param[in]  f  noisy image : the first pointer accounts for the channel
 *             and the second one for the pixel position.
 * @param[in]  taulambda  step-size parameter multiplied by the trade-off
 *             parameter that controls the balancing between the fidelity term
 *             and the TV-regularization term.
 * @param[in]  num_channels  number of channels of the image.
 * @param[in]  dim  image size.
 *
 */

void proxG(float **u, float **v, float **f, float taulambda, int num_channels,
           int dim);

/**
 * \brief  Compute the solution to the proximal mapping of the TV-regularization
 *         term using collaborative norms. Unless stated otherwise, p denotes
 *         the channels, q denotes the spatial derivatives and r denotes the
 *         pixels of the images.
 *
 * @param[out] gx, gy  solution to the proximal problem : the first pointer
 *             accounts for the channel and the second one for the pixel
 *             position.
 * @param[in]  vx, vy  argument of the proximal problem : the first pointer
 *             accounts for the channel and the second one for the pixel
 *             position.
 * @param[in]  sigma  step-size parameter.
 * @param[in]  p, q, r  collaborative norm indices used for vectorial TV.
 * @param[in]  num_channels  number of channels of the image.
 * @param[in]  dim  image size.
 *
 */

void proxF(float **gx, float **gy, float **vx, float **vy, float sigma,
           int p, int q, int r, int num_channels, int dim);

/**
 * \brief  Compute the solution to the proximal mapping of the TV-regularization
 *         term using Schatten matrix norms. Only
 *         @f$ \left(S^1(col,der),\ell^1(pix)\right) @f$ and
 *         @f$ \left(S^{\infty}(col,der),\ell^1(pix)\right) @f$ are considered.
 *
 * @param[out] gx, gy  solution to the proximal problem : the first pointer
 *             accounts for the channel and the second one for the pixel
 *             position.
 * @param[in]  vx, vy  argument of the proximal problem : the first pointer
 *             accounts for the channel and the second one for the pixel
 *             position.
 * @param[in]  sigma  step-size parameter.
 * @param[in]  p  index of the Schatten @f$ S^p @f$ norm. Only
 *             @f$ p\in\{1,\infty\} @f$ considered.
 * @param[in]  num_channels  number of channels of the image.
 * @param[in]  dim  image size.
 *
 */

void Spl1(float **gx, float **gy, float **vx, float **vy, float sigma, int p,
          int num_channels, int dim);

/**
 * \brief  Compute the projection of a one-dimensional vector onto the unit
 *         @f$ \ell^1 @f$ norm ball.
 *
 * @param[in]   u  input vector.
 * @param[out]  u  projected vector.
 * @param[in]   dim  vector size.
 *
 */

void l1projection(float *u, int dim);

/**
 * \brief  Compute the projection of a two-dimensional matrix onto the unit
 *         @f$ \ell^{1,2} @f$ norm ball.
 *
 * @param[in]   projx, projy  input matrix.
 * @param[out]  projx, projy  projected matrix.
 * @param[in]   sigma  step-size parameter of related proximal mapping.
 * @param[in]   dim  vector size.
 *
 */

void l12projection(float *projx, float *projy, float sigma, int dim);

#endif
