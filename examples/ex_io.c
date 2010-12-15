#include <stdlib.h>
#include <stdio.h>
#include <latan/latan_io.h>

int main(void)
{
    mat *prop;

    prop = mat_create(32,1);
    
    prop_load(prop,ch_PP,0,0,ss_G,ss_G,"ex_io.xml");
    mat_print(prop);
    
    mat_destroy(prop);

    return EXIT_SUCCESS;
}