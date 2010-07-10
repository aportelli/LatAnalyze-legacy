#include <stdlib.h>
#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_plot.h>

int main(void)
{
	plot *p;
	
	p = plot_create();
	
	plot_add_plot(p,"x**2");
	plot_disp(p);
	plot_print(p);
	
	plot_destroy(p);
	
	return EXIT_SUCCESS;
}