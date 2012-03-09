#include <stdlib.h>
#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_plot.h>
#include <latan/latan_tabfunc.h>

const double sq_table[] =
{
    1.0,\
    4.0,\
    9.0,\
    16.0
};
tabfunc sq={1.0,4.0,sq_table,4};

int main(void)
{
    plot *p;
    double xerr[2];
    
    p = plot_create();

    xerr[0] = 0.2;
    xerr[1] = 0.5;
    
    plot_add_vlineaerr(p,3,xerr,"rgb 'green'");
    plot_add_plot(p,"x**2 lc rgb 'red'");
    plot_add_func(p,&tabfunc_lineval,&sq,sq.xmin,sq.xmax,1000,"x**2 table",\
                  "rgb 'blue'");
    plot_set_scale_ymanual(p,0.0,16.0);
    plot_disp(p);
    plot_print(p);
    
    plot_destroy(p);
    
    return EXIT_SUCCESS;
}
