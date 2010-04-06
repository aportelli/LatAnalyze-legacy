#include <stdio.h>
#include <latan/rand.h>

#define GENTEST_SEQ_LENGTH 25
#define GENTEST_SAVE_STEP 9
#define DIS_SEQ_LENGTH 1000000
#define DICE_NFACE 6

int main(void)
{
	int dice,hist_dice[DICE_NFACE];
	int i,f;
	double randd;
	rand_gen_state state;
	
	printf("- GENERATOR STATE I/O TESTS\n");
	rand_timeinit();
	printf("-- generating a %d steps random sequence...\n",GENTEST_SEQ_LENGTH);
	for (i=0;i<GENTEST_SEQ_LENGTH;i++) 
	{
		if (i == GENTEST_SAVE_STEP)
		{
			rand_get_gen_state(state);
			rand_save_gen_state("ex_rand",state);
			printf("generator state after step %d saved in ex_rand.rand\n",\
				   GENTEST_SAVE_STEP-1);
		}
		randd = rand_u(0.0,1.0);
		printf("step %d\t: %e\n",i,randd);
	}
	printf("-- messing up the generator...\n");
	rand_timeinit();
	rand_get_gen_state(state);
	printf("-- reloading state from ex_rand.rand...\n");
	rand_load_gen_state(state,"ex_rand");
	rand_set_gen_state(state);
	printf("-- generating a %d steps random sequence...\n",GENTEST_SEQ_LENGTH);
	for (i=0;i<GENTEST_SEQ_LENGTH;i++) 
	{
		randd = rand_u(0.0,1.0);
		printf("step %d\t: %e\n",i,randd);
	}
	printf("- DISTRIBUTIONS TESTS\n");
	rand_timeinit();
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
	
	return EXIT_SUCCESS;
}