#include <latan/latan_minimizer.h>

double polynom(const mat *var, void *param);

double polynom(const mat *var, void *param)
{
    double x;
    
    param = NULL;
    x = mat_get(var,0,0);
    
    return (x-1.0)*(x-2.0)*(x-3.0)*(x-4.1);
}

int main(void)
{
    mat *res;
    double fval;
    
    res = mat_create(1,1);
    
    mat_set(res,0,0,2.5);
    
    minimize(res,&fval,&polynom,NULL);
    
    printf("%f\n",mat_get(res,0,0));
    
    return 0;
}
