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
 * @file libPDHG.cpp
 * @brief Implementation of the primal-dual hybrid gradient algorithm for
 *        minimizing the functional composed of an @f$ \ell^2 @f$ fidelity term
 *        and a vectorial total variation regularization term that penalizes
 *        the discrete gradient of the underlying image by means of
 *        collaborative norms: @f$ \ell^{p,q,r} @f$ or @f$ (S^{p},\ell^q) @f$.
 * @author Joan Duran <joan.duran@uib.es>
 */

#include "libPDHG.h"

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
         int num_channels, int width, int height)
{
    // Image size
    int dim = width * height;
    
    // Auxiliar vectors
    float **u_upd = new float*[num_channels];
    float **gx = new float*[num_channels];
    float **gy = new float*[num_channels];
    float **gx_upd = new float*[num_channels];
    float **gy_upd = new float*[num_channels];
    float **qx = new float*[num_channels];
    float **qy = new float*[num_channels];
    float **qx_upd = new float*[num_channels];
    float **qy_upd = new float*[num_channels];
    float **v = new float*[num_channels];
    float **vx = new float*[num_channels];
    float **vy = new float*[num_channels];
    float **gradx = new float*[num_channels];
    float **grady = new float*[num_channels];
    float **gradx_upd = new float*[num_channels];
    float **grady_upd = new float*[num_channels];
    float **gradx_ubar = new float*[num_channels];
    float **grady_ubar = new float*[num_channels];
    float **div = new float*[num_channels];
    float **div_upd = new float*[num_channels];
    
    for(int k = 0; k < num_channels; k++)
    {
        u_upd[k] = new float[dim];
        gx[k] = new float[dim];
        gy[k] = new float[dim];
        gx_upd[k] = new float[dim];
        gy_upd[k] = new float[dim];
        qx[k] = new float[dim];
        qy[k] = new float[dim];
        qx_upd[k] = new float[dim];
        qy_upd[k] = new float[dim];
        v[k] = new float[dim];
        vx[k] = new float[dim];
        vy[k] = new float[dim];
        gradx[k] = new float[dim];
        grady[k] = new float[dim];
        gradx_upd[k] = new float[dim];
        grady_upd[k] = new float[dim];
        gradx_ubar[k] = new float[dim];
        grady_ubar[k] = new float[dim];
        div[k] = new float[dim];
        div_upd[k] = new float[dim];
    }
    
    // Initializations
    for(int k = 0; k < num_channels; k++)
    {
        fpClear(u[k], 0.0f, dim);
        fpClear(gx[k], 0.0f, dim);
        fpClear(gy[k], 0.0f, dim);
        fpClear(qx[k], 0.0f, dim);
        fpClear(qy[k], 0.0f, dim);
        fpClear(gradx[k], 0.0f, dim);
        fpClear(grady[k], 0.0f, dim);
        fpClear(div[k], 0.0f, dim);
        }
    
    // Backtracking parameters
    float s = 1.0f;
    float gamma = 0.75f;
    float beta = 0.95f;
    float alpha0 = 0.2f;
    float alpha = alpha0;
    float delta = 1.5f;
    float eta = 0.95f;
    
    // PDHG algorithm parameters
    float taulambda = tau * lambda;
    float divtau = 1.0f / tau;
    float divsigma = 1.0f / sigma;
    float theta1 = 1.0f + theta;
    
    // Apply Primal-Dual Hybrid Gradient scheme
    int iter = 0;
    float residual = fLarge;
    
    while((iter < maxIter) && (residual >= tol))
    {
        // Argument of proximal mapping of fidelity term
        for(int k = 0; k < num_channels; k++)
            for(int i = 0; i < dim; i++)
                v[k][i] = u[k][i] + tau * div[k][i];
        
        // Proximal solution of fidelity term
        proxG(u_upd, v, f, taulambda, num_channels, dim);
        
        // Gradient of updated primal variable
        gradient(u_upd, gradx_upd, grady_upd, num_channels, width, height);
        
        // Argument of proximal mapping of regularization term
        for(int k = 0; k < num_channels; k++)
            for(int i = 0; i < dim; i++)
            {
                float ubarx = theta1 * gradx_upd[k][i] - theta * gradx[k][i];
                float ubary = theta1 * grady_upd[k][i] - theta * grady[k][i];
                
                vx[k][i] = ubarx + divsigma * qx[k][i];
                vy[k][i] = ubary + divsigma * qy[k][i];
                
                gradx_ubar[k][i] = ubarx;
                grady_ubar[k][i] = ubary;
            }
        
        proxF(gx_upd, gy_upd, vx, vy, sigma, p, q, r, num_channels, dim);
        
        // Update dual variable
        for(int k = 0; k < num_channels; k++)
            for(int i = 0; i < dim; i++)
            {
                qx_upd[k][i] = qx[k][i] + sigma * (gradx_ubar[k][i]
                                                   - gx_upd[k][i]);
                qy_upd[k][i] = qy[k][i] + sigma * (grady_ubar[k][i]
                                                   - gy_upd[k][i]);
            }
        
        // Divergence of updated dual variable
        divergence(qx_upd, qy_upd, div_upd, num_channels, width, height);
        
        // Compute primal residual, dual residual, and backtracking condition
        float resprimal = 0.0f;
        float resdual = 0.0f;
        
        float product = 0.0f;
        float unorm = 0.0f;
        float qnorm = 0.0f;
        
        for(int k = 0; k < num_channels; k++)
            for(int i = 0; i < dim; i++)
            {
                float udiff = u[k][i] - u_upd[k][i];
                float qxdiff = qx[k][i] - qx_upd[k][i];
                float qydiff = qy[k][i] - qy_upd[k][i];
                float divdiff = div[k][i] - div_upd[k][i];
                float gradxdiff = gradx[k][i] - gradx_upd[k][i];
                float gradydiff = grady[k][i] - grady_upd[k][i];
                
                resprimal += fabs(divtau * udiff + divdiff);
                resdual += fabs(divsigma * qxdiff - gradxdiff);
                resdual += fabs(divsigma * qydiff - gradydiff);
                
                unorm += (udiff * udiff);
                qnorm += (qxdiff * qxdiff + qydiff * qydiff);
                product += (gradxdiff * qxdiff + gradydiff * qydiff);
            }
        
        float b = (2.0f * tau * sigma * product) / (gamma * sigma * unorm +
                                                    gamma * tau * qnorm);
        
        // Adapt step-size parameters
        float dual_dot_delta = resdual * s * delta;
        float dual_div_delta = (resdual * s) / delta;
        
        if(b > 1)
        {
            // Decrease step-sizes to fit balancing principle
            tau = (beta * tau) / b;
            sigma = (beta * sigma) / b;
            alpha = alpha0;
            
            for(int k = 0; k < num_channels; k++)
            {
                fpCopy(u[k], u_upd[k], dim);
                fpCopy(gx[k], gx_upd[k], dim);
                fpCopy(gy[k], gy_upd[k], dim);
                fpCopy(qx[k], qx_upd[k], dim);
                fpCopy(qy[k], qy_upd[k], dim);
                fpCopy(gradx[k], gradx_upd[k], dim);
                fpCopy(grady[k], grady_upd[k], dim);
                fpCopy(div[k], div_upd[k], dim);
            }
            
        } else if(resprimal > dual_dot_delta)
        {
            // Increase primal step-size and decrease dual step-size
            tau = tau / (1.0f - alpha);
            sigma = sigma * (1.0f - alpha);
            alpha = alpha * eta;
            
        } else if(resprimal < dual_div_delta)
        {
            // Decrease primal step-size and increase dual step-size
            tau = tau * (1.0f - alpha);
            sigma = sigma / (1.0f - alpha);
            alpha = alpha * eta;
        }
        
        // Update variables
        taulambda = tau * lambda;
        divsigma = 1.0f / sigma;
        divtau = 1.0f / tau;
        iter++;
        
        for(int k = 0; k < num_channels; k++)
        {
            fpCopy(u_upd[k], u[k], dim);
            fpCopy(gx_upd[k], gx[k], dim);
            fpCopy(gy_upd[k], gy[k], dim);
            fpCopy(qx_upd[k], qx[k], dim);
            fpCopy(qy_upd[k], qy[k], dim);
            fpCopy(gradx_upd[k], gradx[k], dim);
            fpCopy(grady_upd[k], grady[k], dim);
            fpCopy(div_upd[k], div[k], dim);
        }
        
        // Compute residual at current iteration
        residual = (resprimal + resdual) / ((float) (dim * num_channels));
    }
    
    // Delete allocated memory
    for(int k = 0; k < num_channels; k++)
    {
        delete[] u_upd[k];
        delete[] gx[k];
        delete[] gy[k];
        delete[] gx_upd[k];
        delete[] gy_upd[k];
        delete[] qx[k];
        delete[] qy[k];
        delete[] qx_upd[k];
        delete[] qy_upd[k];
        delete[] v[k];
        delete[] vx[k];
        delete[] vy[k];
        delete[] gradx[k];
        delete[] grady[k];
        delete[] gradx_upd[k];
        delete[] grady_upd[k];
        delete[] gradx_ubar[k];
        delete[] grady_ubar[k];
        delete[] div[k];
        delete[] div_upd[k];
    }
    
    delete[] u_upd;
    delete[] gx;
    delete[] gy;
    delete[] gx_upd;
    delete[] gy_upd;
    delete[] qx;
    delete[] qy;
    delete[] qx_upd;
    delete[] qy_upd;
    delete[] v;
    delete[] vx;
    delete[] vy;
    delete[] gradx;
    delete[] grady;
    delete[] gradx_upd;
    delete[] grady_upd;
    delete[] gradx_ubar;
    delete[] grady_ubar;
    delete[] div;
    delete[] div_upd;
    
	return 1;
}
