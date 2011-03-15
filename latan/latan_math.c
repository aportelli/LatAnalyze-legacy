/* latan_math.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <latan/latan_math.h>
#include <latan/latan_includes.h>

unsigned int binomial(const unsigned int n, const unsigned int p)
{
    unsigned int *b;
    unsigned int i,j;
    unsigned int res;
    
    MALLOC_NOERRET(b,unsigned int*,n+1);
    
    b[0]=1;
    for (i=1;i<=n;i++)
    {
        b[i] = 1;
        for (j=i-1;j>0;j--)
        {
            b[j] += b[j-1];
        }
    }
    res = b[p];
    
    FREE(b);
    
    return res;
}

latan_errno finite_diff(mat *ddat, const mat *dat)
{
    size_t i,j;
    double ddat_ij;
    
    if (nrow(ddat) != nrow(dat) - 2)
    {
        LATAN_ERROR("derivative matrix have wrong dimensions",LATAN_EBADLEN);
    }
    
    FOR_VAL(ddat,i,j)
    {
        ddat_ij = 0.5*(mat_get(dat,i+2,j) - mat_get(dat,i,j));
        mat_set(ddat,i,j,ddat_ij);
    }
    
    return LATAN_SUCCESS;
}

size_t coord_to_rowmaj(const size_t *x, const size_t *dim, const size_t ndim)
{
    size_t d,ind;
    
    ind = dim[1]*x[0];
    
    for (d=1;d<ndim-1;d++)
    {
        ind = dim[d+1]*(x[d] + ind);
    }
    ind += x[ndim-1];
    
    return ind;
}

void rowmaj_to_coord(size_t *x, const size_t *dim, const size_t ndim,\
                     const size_t ind)
{
    size_t j,dimprod;
    int d;
    
    j = ind;
    dimprod = 1;
    
    for (d=(int)(ndim)-1;d>=0;d--)
    {
        x[d] = (j/dimprod)%dim[d];
        j -= dimprod*x[d];
        dimprod *= dim[d];
    }
}
