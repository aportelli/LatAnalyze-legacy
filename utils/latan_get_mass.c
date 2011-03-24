#include <stdlib.h>
#include <stdio.h>
#include <latan/latan_mass.h>

int main(int argc, char *argv[])
{
    double mass[2];
    
    if (argc != 2)
    {
        fprintf(stderr,"usage: %s <particle name>\n",argv[0]);
        return EXIT_FAILURE;
    }
    
    get_exp_mass(mass,argv[1]);
    printf("M_%s = %f +/- %f MeV\n",argv[1],mass[0],mass[1]);
    
    return EXIT_SUCCESS;
}
