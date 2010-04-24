#include <stdio.h>
#include <latan/mat.h>
#include <latan/rand.h>
#include <latan/io.h>

int main(int argc, char *argv[])
{
	mat a,b,c,d,e,f;
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
	mat_print(b);
	printf("c =\n");
	mat_print(c);
	printf("a <- b*c\n");
	mat_mul(a,b,c);
	printf("a =\n");
	mat_print(a);
	printf("d <- a(1:2,0:1)\n");
	mat_cp_subm(d,a,1,0,2,1);
	printf("d =\n");
	mat_print(d);
	printf("a(1,0:1) <- c\n");
	mat_set_subm(a,c,1,0,1,2);
	printf("a =\n");
	mat_print(a);
	printf("a <- random(-6,6)\n");
	mat_rand_u(a,-6.0,6.0);
	printf("a =\n");
	mat_print(a);
	printf("e <- a^(-1)\n");
	mat_inv(e,a);
	printf("e =\n");
	mat_print(e);
	printf("e <- e*a\n");
	mat_eqmul_r(e,a);
	printf("e =\n");
	mat_print(e);
	
	mat_destroy(a);
	mat_destroy(b);
	mat_destroy(c);
	mat_destroy(d);
	mat_destroy(e);
	
	return EXIT_SUCCESS;
}
