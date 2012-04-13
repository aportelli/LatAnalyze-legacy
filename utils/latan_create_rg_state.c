#include <stdlib.h>
#include <latan/latan_rand.h>
#include <latan/latan_io.h>

#define OUTF_NAME argv[1]

int main(int argc, char *argv[])
{
    rg_state state;
    int itime;
    strbuf latan_path;
    
    if (argc != 2)
    {
        fprintf(stderr,"usage: %s <output_file/name>\n",argv[0]);
        return EXIT_FAILURE;
    }

    io_init();
    itime = randgen_init_from_time();
    randgen_get_state(state);
    sprintf(latan_path,"%s%c%d",OUTF_NAME,LATAN_PATH_SEP,itime);
    randgen_save_state(latan_path,'w',state);
    io_finish();

    return EXIT_SUCCESS;
}
