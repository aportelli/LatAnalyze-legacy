#include <stdio.h>
#include <latan/mat.h>

int main(int argc, char *argv[])
{
	mat a,b,c,d;
	
	mat_create(&a,3,2);
	mat_create(&b,3,1);
	mat_create(&c,1,2);
	mat_create(&d,2,2);
	
	mat_set(b,0,0,2.0);
	mat_set(b,2,0,5.0);
	mat_set(c,0,0,1.0);
	mat_set(c,0,1,3.0);
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
	mat_set_subm(a,c,1,0,1,1);
	printf("a =\n");
	mat_print(a);
	
	mat_destroy(&a);
	mat_destroy(&b);
	mat_destroy(&c);
	mat_destroy(&d);
	
	return EXIT_SUCCESS;
}
