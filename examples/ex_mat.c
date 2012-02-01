#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_rand.h>
#include <latan/latan_io.h>

int main(void)
{
    mat *a,*b,*c,*d,*e,*f,*g,*h,*i,*j;
    const double b_init[] =
    {
        2.5,\
        2.0,\
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
    f = mat_create(3,3);
    g = mat_create(5,3);
    h = mat_create(3,5);
    i = mat_create(5,3);
    j = mat_create(3,1);
    
    randgen_init(12);

    printf("---------- %-30s ----------\n","matrix product");
    printf("b =\n");
    mat_print(b,"%8.2f");
    printf("\n");
    printf("c =\n");
    mat_print(c,"%8.2f");
    printf("\n");
    printf("a <- b*c\n");
    mat_mul(a,b,'n',c,'n');
    printf("a =\n");
    mat_print(a,"%8.2f");
    printf("\n");

    printf("\n---------- %-30s ----------\n","sub-matrices");
    printf("d <- a(1:2,0:1)\n");
    mat_get_subm(d,a,1,0,2,1);
    printf("d =\n");
    mat_print(d,"%9.2f");
    printf("\n");
    printf("a(1,0:2) <- c\n");
    mat_set_subm(a,c,1,0,1,2);
    printf("a =\n");
    mat_print(a,"%8.2f");
    printf("\n");

    printf("\n---------- %-30s ----------\n","symmetric matrix product");
    printf("e <- a*t(a)\n");
    mat_mul(e,a,'n',a,'t');
    mat_print(e,"%8.2f");
    printf("\n");
    printf("e is assumed symmetric\n\n");
    mat_assume_sym(e,true);
    printf("j <- e*b\n");
    mat_mul(j,e,'n',b,'t');
    mat_print(j,"%8.2f");
    printf("\n");
    printf("f <- a*e\n");
    mat_mul(f,a,'n',e,'n');
    mat_print(f,"%8.2f");
    printf("\n");
    printf("f <- e*t(a)\n");
    mat_mul(f,e,'n',a,'t');
    mat_print(f,"%8.2f");
    printf("\n");
    printf("g <- random(-1.0,1.0)\n");
    mat_rand_u(g,-1.0,1.0);
    mat_print(g,"%8.2f");
    printf("\n");
    printf("h <- e*t(g) (buffer should be created)\n");
    mat_mul(h,e,'n',g,'t');
    mat_print(h,"%8.2f");
    printf("\n");
    printf("i <- g*e\n");
    mat_mul(i,g,'n',e,'n');
    mat_print(i,"%8.2f");
    printf("\n");
    printf("e is not assumed symmetric anymore\n");
    mat_assume_sym(e,false);
    
    printf("\n---------- %-30s ----------\n","matrix inversion");
    printf("*** LU decomposition\n");
    printf("a <- random(-6.0,6.0)\n");
    mat_rand_u(a,-6.0,6.0);
    printf("a =\n");
    mat_print(a,"%8.2f");
    printf("\n");
    printf("e <- a^(-1)\n");
    mat_inv_LU(e,a);
    printf("e =\n");
    mat_print(e,"%8.2f");
    printf("\n");
    printf("f <- e*a\n");
    mat_mul(f,e,'n',a,'n');
    printf("f =\n");
    mat_print(f,"%8.2f");
    printf("\n");
    printf("*** Cholesky decomposition\n");
    printf("a <- a*t(a)\n");
    mat_mul(a,a,'n',a,'t');
    mat_print(a,"%8.2f");
    printf("\n");
    printf("a is assumed symmetric\n\n");
    mat_assume_sym(a,true);
    printf("e <- a^(-1)\n");
    mat_inv_symChol(e,a);
    printf("e =\n");
    mat_print(e,"%8.2f");
    printf("\n");
    printf("f <- e*a\n");
    mat_mul(f,e,'n',a,'n');
    printf("f =\n");
    mat_print(f,"%8.2f");
    printf("\n");
    
    mat_destroy(a);
    mat_destroy(b);
    mat_destroy(c);
    mat_destroy(d);
    mat_destroy(e);
    mat_destroy(f);
    mat_destroy(g);
    mat_destroy(h);
    mat_destroy(i);
    mat_destroy(j);
    
    return EXIT_SUCCESS;
}

