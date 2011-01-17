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
    mat *var,*real_param,*fit_param;
    size_t i;
    minalg_no alg;
    double step;
    fit_data *d;
    strbuf plotcmd;
    
    step = DRATIO(XMAX,NDATA);
    
    var = mat_create(NDATA,1);
    real_param = mat_create(NPAR,1);
    fit_param = mat_create(NPAR,1);
    d = fit_data_create(NDATA,fm_expdec.ndim);
    p = plot_create();
    
    randgen_init_from_time();
    mat_set(real_param,0,0,0.5);
    mat_set(real_param,1,0,5.0);
    mat_cst(var,SQ(ERR));
    fit_data_set_model(d,&fm_expdec);
    fit_data_set_data_var(d,var);
    printf("-- generating exponential decay data with gaussian errors...\n");
    for (i=0;i<NDATA;i++)
    {
        fit_data_set_x(d,i,0,i*step);
        fit_data_set_data(d,i,fit_data_model_eval(d,i,real_param)\
                          +rand_n(0.0,ERR));
    }
    fit_data_fit_all_points(d,true);
    latan_set_verb(DEBUG);
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
        sprintf(plotcmd,"%e*exp(-%e*x)",mat_get(fit_param,1,0),\
                mat_get(fit_param,0,0));
        plot_add_plot(p,plotcmd);
    }
    mat_eqsqrt(var);
    plot_add_dat_yerr(p,fit_data_pt_x(d),fit_data_pt_data(d),var,"","");
    plot_disp(p);
    plot_destroy(p);
    
    fit_data_destroy(d);
    mat_destroy(var);
    mat_destroy(real_param);
    mat_destroy(fit_param);
    
    return EXIT_SUCCESS;
}