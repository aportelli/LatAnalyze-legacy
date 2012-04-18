/* ex_fit.c, part of LatAnalyze library
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <latan/latan_fit.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_models.h>
#include <latan/latan_plot.h>
#include <latan/latan_rand.h>

#define ERR 0.05
#define NDATA 20
#define XMAX 5
#define NPAR 2
#define LIB MINUIT

int main(void)
{
    plot *p;
    mat *var,*x,*y,*real_param,*fit_param;
    size_t i;
    minalg_no alg;
    double step,buf;
    fit_data *d;
    strbuf plotcmd;
    
    step = DRATIO(XMAX,NDATA);
    
    var        = mat_create(NDATA,1);
    x          = mat_create(NDATA,1);
    y          = mat_create(NDATA,1);
    real_param = mat_create(NPAR,1);
    fit_param  = mat_create(NPAR,1);
    d          = fit_data_create(NDATA,fm_expdec.nxdim,fm_expdec.nydim);
    p          = plot_create();
    
    randgen_init_from_time();
    mat_set(real_param,0,0,0.5);
    mat_set(real_param,1,0,log(5.0));
    mat_cst(var,SQ(ERR));
    fit_data_set_model(d,&fm_expdec,NULL);
    fit_data_set_y_covar(d,0,0,var);
    printf("-- generating exponential decay data with gaussian errors...\n");
    for (i=0;i<NDATA;i++)
    {
        fit_data_set_x(d,i,0,i*step);
        buf = fit_data_model_eval(d,0,i,real_param)+rand_n(0.0,ERR);
        fit_data_set_y(d,i,0,buf);
    }
    fit_data_fit_all_points(d,true);
    latan_set_verb(DEBUG1);
    for (alg=0;alg<NMINALG;alg++)
    {
        mat_set(fit_param,0,0,0.3);
        mat_set(fit_param,1,0,7.0);
        minimizer_set_alg(alg);
        printf("-- fitting datas...\n");
        data_fit(fit_param,d);
        printf("exact parameters=\n");
        mat_print(real_param,"%f");
        printf("\nfit parameters=\n");
        mat_print(fit_param,"%f");
        printf("\nchi2/dof= %e\n",fit_data_get_chi2pdof(d));
        sprintf(plotcmd,"%e*exp(-%e*x)",exp(mat_get(fit_param,1,0)),\
                mat_get(fit_param,0,0));
        plot_add_plot(p,plotcmd);
    }
    mat_eqsqrt(var);
    fit_data_get_x_k(x,d,0);
    fit_data_get_y_k(y,d,0);
    plot_add_dat(p,x,y,NULL,var,"","");
    plot_disp(p);
    plot_destroy(p);
    
    fit_data_destroy(d);
    mat_destroy(var);
    mat_destroy(x);
    mat_destroy(y);
    mat_destroy(real_param);
    mat_destroy(fit_param);
    
    return EXIT_SUCCESS;
}
