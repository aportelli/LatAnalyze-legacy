/* latan_min_minuit2.cpp, part of LatAnalyze library
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

using namespace std;
using namespace ROOT;
using namespace Minuit2;

class Minuit2MinFunc: public FCNBase
{
public:
    Minuit2MinFunc(min_func *init_f, void *init_param);
    ~Minuit2MinFunc(void);
    
    virtual double operator()(const vector<double>& v_var) const;
    virtual double Up(void) const;
    
private:
    min_func *f;
    void *param;
};

Minuit2MinFunc::Minuit2MinFunc(min_func *init_f, void *init_param)
{
    f            = init_f;
    param        = init_param;
}

Minuit2MinFunc::~Minuit2MinFunc(void)
{
}

double Minuit2MinFunc::operator()(const vector<double>& v_x) const 
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

latan_errno minimize_minuit2(mat *x, const mat *x_limit, double *f_min,\
                             min_func *f, void *param)
{
    latan_errno status;
    strbuf buf;
    string name;
    size_t ndim;
    size_t i;
    double x_i,err_i;
    bool is_xl_l_nan,is_xl_u_nan;
    MnUserParameters Init_x;
    MnApplication *Minimizer;
  
    status        = LATAN_SUCCESS;
    ndim          = nrow(x);
    
    for (i=0;i<ndim;i++)
    {
        sprintf(buf,"p%lu",(long unsigned)i);
        name = buf;
        Init_x.Add(name,mat_get(x,i,0),fabs(mat_get(x,i,0))*INIT_RERROR);
        if (x_limit)
        {
            is_xl_l_nan = latan_isnan(mat_get(x_limit,i,0));
            is_xl_u_nan = latan_isnan(mat_get(x_limit,i,1));
            if ((!is_xl_l_nan)&&(is_xl_u_nan))
            {
                Init_x.SetLowerLimit((unsigned int)i,mat_get(x_limit,i,0));
            }
            else if ((is_xl_l_nan)&&(!is_xl_u_nan))
            {
                Init_x.SetUpperLimit((unsigned int)i,mat_get(x_limit,i,1));
            }
            else if ((!is_xl_l_nan)&&(!is_xl_u_nan))
            {
                Init_x.SetLimits((unsigned int)i,mat_get(x_limit,i,0),\
                                 mat_get(x_limit,i,1));
            }
        }
    }

    Minuit2MinFunc F(f,param);
    MnSimplex Simplex1(F,Init_x,0);
    Minimizer = &Simplex1;
    latan_printf(DEBUG2,"(MINUIT) Minimizing...\n");
    FunctionMinimum Min = (*Minimizer)();
    latan_printf(DEBUG2,"(MINUIT) Pre-minimizer call :\n");
    if (latan_get_verb() == DEBUG2)
    {
        cout << "--------------------------------------------------------";
        cout << Min;
        cout << "--------------------------------------------------------";
        cout << endl;
    }
    for (i=0;i<ndim;i++)
    {
        x_i   = Min.UserParameters().Value((unsigned int)i);
        err_i = Min.UserParameters().Error((unsigned int)i);
        Init_x.SetValue((unsigned int)i,x_i);
        Init_x.SetError((unsigned int)i,err_i);
    }
    MnMigrad Migrad2(F,Init_x,2);
    MnSimplex Simplex2(F,Init_x,2);
    switch (minimizer_get_alg())
    {
        case MIN_MIGRAD:
            Minimizer = &Migrad2;
            break;
        case MIN_SIMPLEX:
            Minimizer = &Simplex2;
            break;
        default:
            LATAN_ERROR("invalid MINUIT minimization algorithm flag",
                        LATAN_EINVAL);
            break;
    }
    Min = (*Minimizer)();
    for (i=0;i<ndim;i++)
    {
        x_i = Min.UserParameters().Parameter((unsigned int)i).Value();
        mat_set(x,i,0,x_i);
    }
    if (!Min.IsValid())
    {
        LATAN_WARNING("MINUIT library reported that minimization result is not valid",\
                      LATAN_FAILURE);
        status = LATAN_FAILURE;
    }
    *f_min = Min.Fval();
          
    latan_printf(DEBUG2,"(MINUIT) Minimizer call :\n");
    if (latan_get_verb() == DEBUG2)
    {
        cout << "--------------------------------------------------------";
        cout << Min;
        cout << "--------------------------------------------------------";
        cout << endl;
    }
    latan_printf(DEBUG2,"(MINUIT) Scan around last position :\n");
    if (latan_get_verb() == DEBUG2)
    {
        vector<pair<double, double> > ScanRes;
        MnPlot Plot;

        MnScan Scanner(F,Min.UserParameters(),2);
        cout << "--------------------------------------------------------";
        cout << endl;
        for (i=0;i<ndim;i++)
        {
            cout << "Parameter p" << (int)i << endl;
            ScanRes = Scanner.Scan((unsigned int)i);
            Plot(ScanRes);
        }
        cout << "--------------------------------------------------------";
        cout << endl;
    }
    
    
    return status;
}

#endif
