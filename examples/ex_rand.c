/* ex_rand.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <math.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_math.h>
#include <latan/latan_statistics.h>
#include <latan/latan_rand.h>
#include <latan/latan_plot.h>

#define GENTEST_SEQ_LENGTH 25
#define GENTEST_SAVE_STEP 9
#define DIS_SEQ_LENGTH 10000
#define DICE_NFACE 6
#define GAUSS_MU 0.0
#define GAUSS_SIG 1.0
#define GAUSS_HIST_MAX 5.0
#define HIST_CONT_NINT 100

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

    io_init();
    printf("- GENERATOR STATE I/O TESTS\n");
    randgen_init_from_time();
    printf("-- generating a %d steps random sequence...\n",GENTEST_SEQ_LENGTH);
    for (i=0;i<GENTEST_SEQ_LENGTH;i++) 
    {
        if (i == GENTEST_SAVE_STEP)
        {
            randgen_get_state(state);
            randgen_save_state("ex_rand.seed:ex_rand",'w',state);
            printf("generator state after step %d saved in ex_rand.seed\n",\
                   GENTEST_SAVE_STEP-1);
        }
        randd = rand_u(0.0,1.0);
        printf("step %d\t: %e\n",i,randd);
    }
    printf("-- messing up the generator...\n");
    randgen_init_from_time();
    randgen_get_state(state);
    printf("-- reloading state from ex_rand.seed...\n");
    randgen_load_state(state,"ex_rand.seed:ex_rand");
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
        printf("face %d\t: %f%%\n",f+1,DRATIO(hist_dice[f],\
               DIS_SEQ_LENGTH)*100.0);
    }
    printf("-- NORMAL DISTRIBUTION\n");
    printf("-- Generating %d gaussian numbers with mean %.2f and width %.2f...\n",\
           DIS_SEQ_LENGTH,GAUSS_MU,GAUSS_SIG);
    for (i=0;i<DIS_SEQ_LENGTH;i++)
    {
        mat_set(rseq,i,0,rand_n(GAUSS_MU,GAUSS_SIG));
    }
    printf("-- Making sequence histogram...\n");
    histogram(hist_cont,rseq,NULL,-GAUSS_HIST_MAX,GAUSS_HIST_MAX,\
              HIST_CONT_NINT);
    dist_plot = plot_create();
    plot_add_histogram(dist_plot,hist_cont,-GAUSS_HIST_MAX,GAUSS_HIST_MAX,\
                       DIS_SEQ_LENGTH,true,"rand_n distribution","rgb 'red'");
    sprintf(plotcmd,"exp(-(x-%e)**2/(2.0*%e))/%e t 'theoretical distribution'",\
            GAUSS_MU,SQ(GAUSS_SIG),sqrt(2*C_PI)*GAUSS_SIG);
    plot_add_plot(dist_plot,plotcmd,"");
    plot_disp(dist_plot);
    plot_destroy(dist_plot);
    io_finish();

    mat_destroy(rseq);
    mat_destroy(hist_cont);
    mat_destroy(x);
    
    return EXIT_SUCCESS;
}
