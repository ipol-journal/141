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
 * @file libProximal.cpp
 * @brief Implementation of proximal mappings.
 * @author Joan Duran <joan.duran@uib.es>
 */

#include "libProximal.h"

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
           int dim)
{
    float constant = 1.0f + taulambda;
    
    for(int k = 0; k < num_channels; k++)
        for(int i = 0; i < dim; i++)
            u[k][i] = (v[k][i] + taulambda * f[k][i]) / constant;
}

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
           int p, int q, int r, int num_channels, int dim)
{
    // Constant
    float divsigma = 1.0f / sigma;
    
    if((p == 1) && (q == 1) && (r == 1))
    {
        // $\ell^{1,1,1}$-TV regularization
        for(int k = 0; k < num_channels; k++)
            for(int i = 0; i < dim; i++)
            {
                gx[k][i] = SIGN(vx[k][i]) * MAX(fabs(vx[k][i]) - divsigma,
                                                0.0f);
                gy[k][i] = SIGN(vy[k][i]) * MAX(fabs(vy[k][i]) - divsigma,
                                                0.0f);
            }
        
    } else if((p == 2) && (q == 1) && (r == 1))
    {
        // $\ell^{2,1,1}$-TV regularization
        for(int i = 0; i < dim; i++)
        {
            float normx = 0.0f;
            float normy = 0.0f;
            
            for(int k = 0; k < num_channels; k++)
            {
                float xvalue = vx[k][i];
                float yvalue = vy[k][i];
                
                normx += xvalue * xvalue;
                normy += yvalue * yvalue;
            }
            
            normx = sqrtf(normx);
            normy = sqrtf(normy);
            
            float weightx, weighty;
            
            if(normx > fTiny)
                weightx = MAX(normx - divsigma, 0.0f) / normx;
            else
                weightx = 0.0f;
            
            if(normy > fTiny)
                weighty = MAX(normy - divsigma, 0.0f) / normy;
            else
                weighty = 0.0f;
            
            for(int k = 0; k < num_channels; k++)
            {
                gx[k][i] = vx[k][i] * weightx;
                gy[k][i] = vy[k][i] * weighty;
            }
        }
        
    } else if((p == 2) && (q == 2) && (r == 1))
    {
        // $\ell^{2,2,1}$-TV regularization
        for(int i = 0; i < dim; i++)
        {
            float norm = 0.0f;
            
            for(int k = 0; k < num_channels; k++)
            {
                float xvalue = vx[k][i];
                float yvalue = vy[k][i];
                
                norm += xvalue * xvalue;
                norm += yvalue * yvalue;
            }
            
            norm = sqrtf(norm);
            
            float weight;
            
            if(norm > fTiny)
                weight = MAX(norm - divsigma, 0.0f) / norm;
            else
                weight = 0.0f;
            
            for(int k = 0; k < num_channels; k++)
            {
                gx[k][i] = vx[k][i] * weight;
                gy[k][i] = vy[k][i] * weight;
            }
        }
        
    } else if((p == INFNORM) && (q == 1) && (r == 1))
    {
        // $\ell^{\infty,1,1}$-TV regularization
        float *projx, *projy;
        
#pragma omp parallel for private(projx, projy)
        for(int i = 0; i < dim; i++)
        {
            projx = new float[num_channels];
            projy = new float[num_channels];
            
            for(int k = 0; k < num_channels; k++)
            {
                projx[k] = sigma * fabs(vx[k][i]);
                projy[k] = sigma * fabs(vy[k][i]);
            }
            
            l1projection(projx, num_channels);
            l1projection(projy, num_channels);
            
            for(int k = 0; k < num_channels; k++)
            {
                gx[k][i] = vx[k][i] - divsigma * SIGN(vx[k][i]) * projx[k];
                gy[k][i] = vy[k][i] - divsigma * SIGN(vy[k][i]) * projy[k];
            }
        }
        
        delete[] projx;
        delete[] projy;
        
    } else if((p == INFNORM) && (q == 2) && (r == 1))
    {
        // $\ell^{\infty,2,1}$-TV regularization
        float *projx, *projy;
        
#pragma omp parallel for private(projx, projy)
        for(int i = 0; i < dim; i++)
        {
            projx = new float[num_channels];
            projy = new float[num_channels];
            
            for(int k = 0; k < num_channels; k++)
            {
                projx[k] = fabs(vx[k][i]);
                projy[k] = fabs(vy[k][i]);
            }
            
            l12projection(projx, projy, sigma, num_channels);
            
            for(int k = 0; k < num_channels; k++)
            {
                gx[k][i] = vx[k][i] - divsigma * SIGN(vx[k][i]) * projx[k];
                gy[k][i] = vy[k][i] - divsigma * SIGN(vy[k][i]) * projy[k];
            }
        }
        
        delete[] projx;
        delete[] projy;
        
    } else if((p == INFNORM) && (q == INFNORM) && (r == 1))
    {
        // $\ell^{\infty,\infty,1}$-TV regularization
        int dimproj = 2 * num_channels;
        float *proj;
        
#pragma omp parallel for private(proj)
        for(int i = 0; i < dim; i++)
        {
            proj = new float[dimproj];
            
            for(int k = 0, l = 0; k < num_channels; k++, l+=2)
            {
                proj[l] = sigma * fabs(vx[k][i]);
                proj[l+1] = sigma * fabs(vy[k][i]);
            }
            
            l1projection(proj, dimproj);
            
            for(int k = 0, l = 0; k < num_channels; k++, l+=2)
            {
                gx[k][i] = vx[k][i] - divsigma * SIGN(vx[k][i]) * proj[l];
                gy[k][i] = vy[k][i] - divsigma * SIGN(vy[k][i]) * proj[l+1];
            }
        }
        
        delete[] proj;
        
    } else if((p == 2) && (q == INFNORM) && (r == 1))
    {
        // $\ell^{2,\infty,1}$-TV regularization
        // The order of the dimensions in this case is given by (der, col, pix)
        float *norm, *proj;
        
#pragma omp parallel for private(proj, norm)
        for(int i = 0; i < dim; i++)
        {
            norm = new float[num_channels];
            proj = new float[num_channels];
            
            for(int k = 0; k < num_channels; k++)
            {
                float xvalue = vx[k][i];
                float yvalue = vy[k][i];
                float normvalue = sqrtf(xvalue * xvalue + yvalue * yvalue);
                
                norm[k] = normvalue;
                proj[k] = sigma * normvalue;
            }
            
            l1projection(proj, num_channels);
            
            for(int k = 0; k < num_channels; k++)
            {
                float weight;
                
                if(norm[k] > fTiny)
                    weight = MAX(norm[k] - divsigma * proj[k], 0.0f) / norm[k];
                else
                    weight = 0.0f;
                
                gx[k][i] = vx[k][i] * weight;
                gy[k][i] = vy[k][i] * weight;
            }
        }
        
        delete[] norm;
        delete[] proj;
        
    } else if((q == 1) && (r == 0))
    {
        // $(S^p,\ell^1)$-TV regularization
        Spl1(gx, gy, vx, vy, sigma, p, num_channels, dim);
    }
}

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
          int num_channels, int dim)
{
    // Constants
    float divsigma = 1.0f / sigma;
    
    // Auxiliar vector
    float *proj;
    
    // (S^p, \ell^1) norm decouples at each pixel
#pragma omp parallel for private(proj)
    for(int i = 0; i < dim; i++)
    {
        proj= new float[2];
        
        // Compute matrix $M\in\R^{2\times 2}$
        float M1 = 0.0f;
        float M2 = 0.0f;
        float M3 = 0.0f;
        
        for(int k = 0; k < num_channels; k++)
        {
            float valuex = vx[k][i];
            float valuey = vy[k][i];
            
            M1 += (valuex * valuex);
            M2 += (valuex * valuey);
            M3 += (valuey * valuey);
        }
        
        // Compute eigenvalues of M
        float T = M1 + M3;
        float D = M1 * M3 - M2 * M2;
        float det = sqrtf(MAX((T * T / 4.0f) - D, 0.0f));
        float eig1 = MAX((T / 2.0f) + det, 0.0f);
        float eig2 = MAX((T / 2.0f) - det, 0.0f);
        float sig1 = sqrtf(eig1);
        float sig2 = sqrtf(eig2);
        
        // Compute normalized eigenvectors
        float V1, V2, V3, V4;
        V1 = V2 = V3 = V4 = 0.0f;
        
        if(M2 != 0.0f)
        {
            float v0 = M2;
            float v1 = eig1 - M3;
            float v2 = eig2 - M3;
            
            float mu1 = sqrtf(v0 * v0 + v1 * v1);
            float mu2 = sqrtf(v0 * v0 + v2 * v2);
            
            if(mu1 > fTiny)
            {
                V1 = v1 / mu1;
                V3 = v0 / mu1;
            }
            
            if(mu2 > fTiny)
            {
                V2 = v2 / mu2;
                V4 = v0 / mu2;
            }
            
        } else
        {
            if(M1 > M3)
            {
                V1 = V4 = 1.0f;
                V2 = V3 = 0.0f;
                
            } else
            {
                V1 = V4 = 0.0f;
                V2 = V3 = 1.0f;
            }
        }
        
        // Compute prox_p of the diagonal entries
        float sig1_upd, sig2_upd;
        sig1_upd = sig2_upd = 0.0f;
        
        if(p == 1)
        {
            sig1_upd = MAX(sig1 - divsigma, 0.0f);
            sig2_upd = MAX(sig2 - divsigma, 0.0f);
            
        } else if(p == INFNORM)
        {
            proj[0] = sigma * fabs(sig1);
            proj[1] = sigma * fabs(sig2);
            l1projection(proj, 2);
            
            sig1_upd = sig1 - divsigma * proj[0];
            sig2_upd = sig2 - divsigma * proj[1];
        }
        
        // Compute the diagonal entries of $\widehat{\Sigma}\Sigma^{\dagger}_0$
        if(sig1 > fTiny)
            sig1_upd /= sig1;
        
        if(sig2 > fTiny)
            sig2_upd /= sig2;
        
        // Compute solution
        float t1, t2, t3;
        
        t1 = sig1_upd * V1 * V1 + sig2_upd * V2 * V2;
        t2 = sig1_upd * V1 * V3 + sig2_upd * V2 * V4;
        t3 = sig1_upd * V3 * V3 + sig2_upd * V4 * V4;
        
        for(int k = 0; k < num_channels; k++)
        {
            gx[k][i] = vx[k][i] * t1 + vy[k][i] * t2;
            gy[k][i] = vx[k][i] * t2 + vy[k][i] * t3;
        }
    }
    
    // Delete allocated memory
    delete[] proj;
}

/**
 * \brief  Compute the projection of a one-dimensional vector onto the unit
 *         @f$ \ell^1 @f$ norm ball.
 *
 * @param[in]   u  input vector.
 * @param[out]  u  projected vector.
 * @param[in]   dim  vector size.
 *
 */

void l1projection(float *u, int dim)
{
    float sum = fLarge;
    int num = 0;
    float shrinkfactor = 0.0f;
    
    while(sum > 1.0f)
    {
        sum = 0.0f;
        num = 0;
        
        for(int i = 0; i < dim; i++)
        {
            u[i] = MAX(u[i] - shrinkfactor, 0.0f);
            
            sum += fabs(u[i]);
            
            if(u[i]!= 0.0f)
                num++;
        }
        
        if(num > 0)
            shrinkfactor = (sum - 1.0f) / num;
        else
            break;
    }
}

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

void l12projection(float *projx, float *projy, float sigma, int dim)
{
    float *vx = new float[dim];
    float *vy = new float[dim];
    
    fpCopy(projx, vx, dim);
    fpCopy(projy, vy, dim);
    
    float Ax, Ay, b1x, b1y;
    Ax = Ay = b1x = b1y = 0.0f;
    
    for(int i = 0; i < dim; i++)
    {
        Ax += vx[i];
        Ay += vy[i];
    }
    
    float *b2x = new float[dim];
    float *b2y = new float[dim];
    
    fpClear(b2x, 0.0f, dim);
    fpClear(b2y, 0.0f, dim);
    
    float *zx = new float[dim];
    float *zy = new float[dim];
    
    for(int iter = 0; iter < 100; iter++)
    {
        float dx = Ax + b1x;
        float dy = Ay + b1y;
        
        float norm = sqrtf(dx * dx + dy * dy);
        
        if(norm > 1.0f)
        {
            dx /= norm;
            dy /= norm;
        }
        
        for(int i = 0; i < dim; i++)
        {
            zx[i] = MAX(projx[i] + b2x[i], 0.0f);
            zy[i] = MAX(projy[i] + b2y[i], 0.0f);
        }
        
        for(int i = 0; i < dim; i++)
        {
            float sumx = 0.0f;
            float sumy = 0.0f;
            
            for(int j = 0; j < dim; j++)
            {
                float alpha;
                
                if(i == j)  alpha = 0.4f;
                else  alpha = -0.1f;
                
                sumx += alpha * (sigma * vx[j] + dx - b1x + zx[j] - b2x[j]);
                sumy += alpha * (sigma * vy[j] + dy - b1y + zy[j] - b2y[j]);
            }
            
            projx[i] = sumx;
            projy[i] = sumy;
        }
        
        Ax = 0.0f;
        Ay = 0.0f;
        
        for(int i = 0; i < dim; i++)
        {
            Ax += projx[i];
            Ay += projy[i];
        }
        
        b1x = b1x + Ax - dx;
        b1y = b1y + Ay - dy;
        
        for(int i = 0; i < dim; i++)
        {
            b2x[i] = b2x[i] + projx[i] - zx[i];
            b2y[i] = b2y[i] + projy[i] - zy[i];
        }
    }
    
    delete[] vx;
    delete[] vy;
    delete[] b2x;
    delete[] b2y;
    delete[] zx;
    delete[] zy;
}


