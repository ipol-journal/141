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
 * @file libPDHG.h
 * @brief Implementation of the primal-dual hybrid gradient algorithm for
 *        minimizing the functional composed of an @f$ \ell^2 @f$ fidelity term
 *        and a vectorial total variation regularization term that penalizes
 *        the discrete gradient of the underlying image by means of
 *        collaborative norms: @f$ \ell^{p,q,r} @f$ or @f$ (S^{p},\ell^q) @f$.
 * @author Joan Duran <joan.duran@uib.es>
 */

#ifndef _LIBPDHG_H_
#define _LIBPDHG_H_

#include "libProximal.h"
#include "libauxiliar.h"

/**
 * \brief   Apply the primal-dual hybrid gradient method to find the solution
 *          of the vectorial TV problem using collaborative norms.
 *
 * @param[in]   f  noisy image : the first pointer accounts for the number of
 *              channels and the second one for the pixel position.
 * @param[out]  u  denoised image : the first pointer accounts for the number
 *              of channels and the second one for the pixel position.
 * @param[in]   tau  initial step-size parameter for the fidelity term.
 * @param[in]   sigma  initial step-size parameter for the TV-regularization
 *              term.
 * @param[in]   lambda  balancing parameter that controls the trade-off
 *              between the fidelity term and the TV-regularization term.
 * @param[in]   p, q, r  collaborative norm indices used for vectorial TV.
 * @param[in]   tol  tolerance of the PDHG algorithm.
 * @param[in]   maxIter  maximum number of iterations of the PDHG algorithm.
 * @param[in]   num_channels  number of channels of the image.
 * @param[in]   width, height  image size.
 * @return 1 if exit success.
 *
 */

int PDHG(float **f, float **u, float tau, float sigma, float theta,
         float lambda, int p, int q, int r, float tol, int maxIter,
         int num_channels, int width, int height);

#endif
