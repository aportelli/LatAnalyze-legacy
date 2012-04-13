/* ex_stat.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
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

#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define SEQ_LENGTH 1000
#define NDIM 4
#define NBOOT 2000

int main(void)
{
    mat **gvec;
    mat *mean,*var;
    rs_sample *s_mean;
    size_t i,j;
    double sigma;
    
    mean    = mat_create(NDIM,1);
    s_mean  = rs_sample_create(NDIM,NBOOT);
    var     = mat_create(NDIM,NDIM);
    gvec    = mat_ar_create_from_dim(SEQ_LENGTH,mean);
    randgen_init_from_time();
    
    printf("-- generating %d gaussian %d-vectors...\n",SEQ_LENGTH,NDIM);
    for (j=0;j<NDIM;j++)
    {
        sigma = DRATIO(j+1,5);
        printf("dimension %d variance\t: %f\n",(int)j,SQ(sigma));
        for (i=0;i<SEQ_LENGTH;i++)
        {
            mat_set(gvec[i],j,0,rand_n(0.0,sigma));
        }
    }
    
    printf("-- computing mean...\n");
    mat_mean(mean,gvec,SEQ_LENGTH);
    printf("mean =\n");
    mat_print(mean,"%f");
    printf("\n");
    
    printf("-- computing variance...\n");
    mat_var(var,gvec,SEQ_LENGTH);
    printf("variance =\n");
    mat_print(var,"%f");
    printf("\n");
    
    printf("-- resampling mean...\n");
    resample(s_mean,gvec,SEQ_LENGTH,&rs_mean,BOOT,NULL);
    
    printf("-- computing variance from resampled sample...\n");
    rs_sample_var(var,s_mean);
    mat_eqmuls(var,SEQ_LENGTH);
    printf("variance =\n");
    mat_print(var,"%f");
    printf("\n");
    
    mat_destroy(mean);
    rs_sample_destroy(s_mean);
    mat_destroy(var);
    mat_ar_destroy(gvec,SEQ_LENGTH);
    
    return EXIT_SUCCESS;
}
