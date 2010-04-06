#include <stdio.h>
#include <latan/mat.h>
#include <latan/rand.h>
#include <latan/io.h>

int main(int argc, char *argv[])
{
	mat a,b,c,d,e,f;
	
	mat_create(&a,3,3);
	mat_create(&b,3,1);
	mat_create(&c,1,3);
	mat_create(&d,2,2);
	mat_create(&e,3,3);
	
	rand_timeinit();
	
	mat_set(b,0,0,2.0);
	mat_set(b,1,0,1.5);
	mat_set(b,2,0,5.0);
	mat_set(c,0,0,1.0);
	mat_set(c,0,1,10.0);
	mat_set(c,0,2,6.7);
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
	
	mat_destroy(&a);
	mat_destroy(&b);
	mat_destroy(&c);
	mat_destroy(&d);
	mat_destroy(&e);
	
	return EXIT_SUCCESS;
}
