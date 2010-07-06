#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
	plot p;
	mat x, data, var, real_param, fit_param;
	size_t i;
	double step;
	fit_data d;
	stringbuf plotfmt,plotcmd;
	
	step = DRATIO(XMAX,NDATA);
	
	var = mat_create(NDATA,1);
	real_param = mat_create(NPAR,1);
	fit_param = mat_create(NPAR,1);
	d = fit_data_create(NDATA,fm_expdec.ndim);
	
	randgen_init_from_time();
	mat_set(real_param,0,0,0.5);
	mat_set(real_param,1,0,5.0);
	mat_set(fit_param,0,0,0.32);
	mat_set(fit_param,1,0,6.1);
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
	printf("-- fitting datas...\n");
	data_fit(fit_param,d);
	printf("exact parameters=\n");
	mat_print(real_param);
	printf("\nfit parameters=\n");
	mat_print(fit_param);
	printf("\nchi2/dof= %e\n",fit_data_get_chi2pdof(d));
	p = plot_create();
	fit_model_get_plot_fmt(plotfmt,fit_data_pt_model(d));
	sprintf(plotcmd,plotfmt,mat_get(fit_param,0,0),\
			mat_get(fit_param,1,0));
	mat_eqsqrt(var);
	plot_add_plot(p,plotcmd);
	plot_add_dat_yerr(p,fit_data_pt_x(d),fit_data_pt_data(d),var,"","");
	plot_disp(p);
	plot_destroy(p);
	
	fit_data_destroy(d);
	mat_destroy(var);
	mat_destroy(real_param);
	mat_destroy(fit_param);
	
	return EXIT_SUCCESS;
}