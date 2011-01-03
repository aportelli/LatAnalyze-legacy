#include <stdio.h>
#include <math.h>
#include <latan/latan_mat.h>
#include <latan/latan_math.h>
#include <latan/latan_statistics.h>
#include <latan/latan_rand.h>
#include <latan/latan_plot.h>

#define GENTEST_SEQ_LENGTH 25
#define GENTEST_SAVE_STEP 9
#define DIS_SEQ_LENGTH 1000000
#define DICE_NFACE 6
#define GAUSS_MU 0.0
#define GAUSS_SIG 1.0
#define GAUSS_HIST_MAX 5.0
#define HIST_CONT_NINT 500

int main(void)
{
    int dice,hist_dice[DICE_NFACE];
    int i,f;
    double randd;
    rg_state state;
    mat *rseq,*hist_cont,*x;
    plot *dist_plot;
    strbuf plotcmd;
    
    rseq = mat_create(DIS_SEQ_LENGTH,1);
    hist_cont = mat_create(HIST_CONT_NINT,1);
    x = mat_create(HIST_CONT_NINT,1);
    
    printf("- GENERATOR STATE I/O TESTS\n");
    randgen_init_from_time();
    printf("-- generating a %d steps random sequence...\n",GENTEST_SEQ_LENGTH);
    for (i=0;i<GENTEST_SEQ_LENGTH;i++) 
    {
        if (i == GENTEST_SAVE_STEP)
        {
            randgen_get_state(state);
            randgen_save_state("ex_rand",state);
            printf("generator state after step %d saved in ex_rand.rand\n",\
                   GENTEST_SAVE_STEP-1);
        }
        randd = rand_u(0.0,1.0);
        printf("step %d\t: %e\n",i,randd);
    }
    printf("-- messing up the generator...\n");
    randgen_init_from_time();
    randgen_get_state(state);
    printf("-- reloading state from ex_rand.rand...\n");
    randgen_load_state(state,"ex_rand.rand");
    randgen_set_state(state);
    printf("-- generating a %d steps random sequence...\n",GENTEST_SEQ_LENGTH);
    for (i=0;i<GENTEST_SEQ_LENGTH;i++) 
    {
        randd = rand_u(0.0,1.0);
        printf("step %d\t: %e\n",i,randd);
    }
    
    printf("- DISTRIBUTIONS TESTS\n");
    printf("-- DISCRET UNIFORM DISTRIBUTION\n");
    randgen_init_from_time();
    for (f=0;f<DICE_NFACE;f++)
    {
        hist_dice[f] = 0;
    }
    printf("-- Throwing a %d faces dice %d times...\n",DICE_NFACE,\
           DIS_SEQ_LENGTH);
    for (i=0;i<DIS_SEQ_LENGTH;i++)
    {
        dice = rand_ud(DICE_NFACE);
        for (f=0;f<DICE_NFACE;f++)
        {
            if (dice == f)
            {
                hist_dice[f]++;
            }
        }
    }
    printf("distribution :\n");
    for (f=0;f<DICE_NFACE;f++)
    {
        printf("face %d\t: %f%%\n",f+1,DRATIO(hist_dice[f],DIS_SEQ_LENGTH)*100.0);
    }
    printf("-- NORMAL DISTRIBUTION\n");
    printf("-- Generating %d gaussian numbers with mean %.2f and width %.2f...\n",\
           DIS_SEQ_LENGTH,GAUSS_MU,GAUSS_SIG);
    for (i=0;i<DIS_SEQ_LENGTH;i++)
    {
        mat_set(rseq,i,0,rand_n(GAUSS_MU,GAUSS_SIG));
    }
    printf("-- Making sequence histogram...\n");
    histogram(hist_cont,rseq,-GAUSS_HIST_MAX,GAUSS_HIST_MAX,HIST_CONT_NINT);
    mat_eqmuls(hist_cont,\
               1.0/DIS_SEQ_LENGTH*HIST_CONT_NINT/(2.0*GAUSS_HIST_MAX));
    dist_plot = plot_create();
    mat_set_step(x,-GAUSS_HIST_MAX,2.0*GAUSS_HIST_MAX/HIST_CONT_NINT);
    plot_add_dat(dist_plot,x,hist_cont,"rand_n distribution","");
    sprintf(plotcmd,"exp(-(x-%e)**2/(2.0*%e))/%e title \"theoretical distribution\"",\
            GAUSS_MU,SQ(GAUSS_SIG),sqrt(2*C_PI)*GAUSS_SIG);
    plot_add_plot(dist_plot,plotcmd);
    plot_disp(dist_plot);
    plot_destroy(dist_plot);
    
    mat_destroy(rseq);
    mat_destroy(hist_cont);
    mat_destroy(x);
    
    return EXIT_SUCCESS;
}