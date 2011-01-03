#include <stdlib.h>
#include <latan/latan_rand.h>
#include <latan/latan_io.h>

#define OUTF_NAME argv[1]

int main(int argc, char *argv[])
{
    rg_state state;
    
    if (argc != 2)
    {
        fprintf(stderr,"usage: %s <output_file>\n",argv[0]);
        return EXIT_FAILURE;
    }
    
    randgen_init_from_time();
    randgen_get_state(state);
    randgen_save_state(OUTF_NAME,state);
    
    return EXIT_SUCCESS;
}