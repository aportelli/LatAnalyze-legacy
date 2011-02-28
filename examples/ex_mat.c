#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_rand.h>
#include <latan/latan_io.h>

int main(void)
{
    mat *a,*b,*c,*d,*e;
    const double b_init[] =
    {
        2.5,\
        1.0,\
        5.0
    };
    const double c_init[] =
    {
        1.0,10.0,6.7
    };
    
    a = mat_create(3,3);
    b = mat_create_from_ar(b_init,3,1);
    c = mat_create_from_ar(c_init,1,3);
    d = mat_create(2,2);
    e = mat_create(3,3);
    
    randgen_init_from_time();
    printf("b =\n");
    mat_print(b,"%f");
    printf("c =\n");
    mat_print(c,"%f");
    printf("a <- b*c\n");
    mat_mul(a,b,'n',c,'n');
    printf("a =\n");
    mat_print(a,"%f");
    printf("d <- a(1:2,0:1)\n");
    mat_get_subm(d,a,1,0,2,1);
    printf("d =\n");
    mat_print(d,"%f");
    printf("a(1,0:1) <- c\n");
    mat_set_subm(a,c,1,0,1,2);
    printf("a =\n");
    mat_print(a,"%f");
    printf("a <- random(-6,6)\n");
    mat_rand_u(a,-6.0,6.0);
    printf("a =\n");
    mat_print(a,"%f");
    printf("e <- a^(-1)\n");
    mat_inv(e,a);
    printf("e =\n");
    mat_print(e,"%f");
    printf("e <- e*a\n");
    mat_mul(e,e,'n',a,'n');
    printf("e =\n");
    mat_print(e,"%f");
    
    printf("a : ncol= %d tda= %d\n",(int)ncol(a),(int)a->data_cpu->tda);
    mat_print(a,"%f");
    mat_rand_u(e,-6.0,6.0);
    printf("e: ncol= %d tda= %d\n",(int)ncol(e),(int)e->data_cpu->tda);
    mat_print(e,"%f");
    mat_eqadd(e,a);
    printf("e <- e + a\n");
    mat_print(e,"%f");

    mat_destroy(a);
    mat_destroy(b);
    mat_destroy(c);
    mat_destroy(d);
    mat_destroy(e);
    
    return EXIT_SUCCESS;
}
