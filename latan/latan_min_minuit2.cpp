#include <latan/latan_includes.h>
#ifdef HAVE_MINUIT2
#include <latan/latan_min_minuit2.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/VariableMetricMinimizer.h>
#include <Minuit2/SimplexMinimizer.h>
#include <Minuit2/ScanMinimizer.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnPlot.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnSimplex.h>

#ifndef INIT_RERROR
#define INIT_RERROR 0.5
#endif
#ifndef STRATEGY
#define STRATEGY 2
#endif
#ifndef FIT_TOL
#define FIT_TOL 1.0e-2
#endif 

class Minuit2MinFunc: public ROOT::Minuit2::FCNBase
{
public:
    Minuit2MinFunc(min_func *init_f, void *init_param);
    ~Minuit2MinFunc(void);
    
    virtual double operator()(const std::vector<double>& v_var) const;
    virtual double Up(void) const;
    
private:
    min_func *f;
    void *param;
};

Minuit2MinFunc::Minuit2MinFunc(min_func *init_f, void *init_param)
{
    f     = init_f;
    param = init_param;
}

Minuit2MinFunc::~Minuit2MinFunc(void)
{
}

double Minuit2MinFunc::operator()(const std::vector<double>& v_x) const
{
    size_t i;
    double res;
    mat *x_buf;
   
    x_buf = mat_create(v_x.size(),1);
    
    for (i=0;i<nrow(x_buf);i++)
    {
        mat_set(x_buf,i,0,v_x[i]);
    }
    res = f(x_buf,param);

    mat_destroy(x_buf);
    
    return res;
}

double Minuit2MinFunc::Up(void) const
{
    return 1.0;
}

latan_errno minimize_minuit2(mat *x, double *f_min, min_func *f, void *param)
{
    latan_errno status;
    strbuf buf;
    std::string name;
    size_t ndim;
    size_t i;
    double x_i;
    ROOT::Minuit2::MnUserParameters Init_x;
    ROOT::Minuit2::MnApplication *Minimizer;
  
    status        = LATAN_SUCCESS;
    ndim          = nrow(x);
    
    for (i=0;i<ndim;i++)
    {
        sprintf(buf,"p%lu",(long unsigned)i);
        name = buf;
        Init_x.Add(name,mat_get(x,i,0),fabs(mat_get(x,i,0))*INIT_RERROR);
    }

    Minuit2MinFunc F(f,param);
    ROOT::Minuit2::MnMigrad  Migrad(F,Init_x,STRATEGY);
    ROOT::Minuit2::MnSimplex Simplex(F,Init_x,STRATEGY);
    switch (minimizer_get_alg())
    {
        case MIN_MIGRAD:
            Minimizer = &Migrad;
            break;
        case MIN_SIMPLEX:
            Minimizer = &Simplex;
            break;
        default:
            LATAN_ERROR("invalid MINUIT minimization algorithm flag",
                        LATAN_EINVAL);
            break;
    }
    latan_printf(DEBUG,"(MINUIT) Minimizing...\n");
    ROOT::Minuit2::FunctionMinimum Min = (*Minimizer)();
    if (!Min.IsValid())
    {
        LATAN_WARNING("MINUIT library reported that minimization result is not valid",\
                      LATAN_FAILURE);
        status = LATAN_FAILURE;
    }
    for (i=0;i<ndim;i++)
    {
        x_i = Min.UserParameters().Parameter((unsigned int)i).Value();
        mat_set(x,i,0,x_i);
    }
    *f_min = Min.Fval();
          
    latan_printf(DEBUG,"(MINUIT) Scan around last position :\n");
    if (latan_get_verb() == DEBUG)
    {
        std::vector<std::pair<double, double> > ScanRes;
        ROOT::Minuit2::MnPlot Plot;

        ROOT::Minuit2::MnScan DScanner(F,Min.UserParameters(),STRATEGY);
        std::cout << "--------------------------------------------------------";
        std::cout << std::endl;
        for (i=0;i<ndim;i++)
        {
            std::cout << "Parameter p" << (int)i << std::endl;
            ScanRes = DScanner.Scan((unsigned int)i);
            Plot(ScanRes);
        }
        std::cout << "--------------------------------------------------------";
        std::cout << std::endl;
    }
    latan_printf(DEBUG,"(MINUIT) Minimizer call :\n");
    if (latan_get_verb() == DEBUG)
    {
        std::cout << "--------------------------------------------------------";
        std::cout << Min;
        std::cout << "--------------------------------------------------------";
        std::cout << std::endl;
    }
    
    return status;
}

#endif
