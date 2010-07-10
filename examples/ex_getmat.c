#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_io.h>

int main(int argc, char *argv[])
{
    mat *a,*b;
	int dat_nrow;
	
	dat_nrow = mat_load_nrow("FOO","id1","ex_getmat_f1");
	a = mat_create((size_t)(dat_nrow),2);
	dat_nrow = mat_load_nrow("BAR","testid","ex_getmat_f1");
	b = mat_create((size_t)(dat_nrow),2);
	printf("Getting a from FOO id1 in file ex_getmat_f1...\n");
	mat_load(a,"FOO","id1","ex_getmat_f1");
	printf("Getting b from BAR testid in file ex_getmat_f1...\n\n");
	mat_load(b,"BAR","testid","ex_getmat_f1");
	printf("a =\n");
	mat_print(a);
	printf("\nb =\n");
	mat_print(b);
	printf("\n");
	
	mat_destroy(a);
	mat_destroy(b);
	
	mat **c;
	stringbuf ffname;
	int nfile, i;
	
	nfile = get_nfile("ex_getmat_man");
	
	get_firstfname(ffname,"ex_getmat_man");
	dat_nrow = mat_load_nrow("FOO","id1",ffname);
	c = mat_ar_create((size_t)(nfile),(size_t)(dat_nrow),2);
	printf("Filling array c from FOO id1 matrices in files referenced in file ex_getmat_man...\n\n");
	mat_load_ar(c,"FOO","id1","ex_getmat_man");
	for (i=0;i<nfile;i++)
	{
		printf("c[%i] =\n",i);
		mat_print(c[i]);
		printf("\n");
	}
	
	mat_ar_destroy(c,nfile);
	
	return EXIT_SUCCESS;
}