#include <stdio.h>
#include <latan/rand.h>

#define SEQ_LENGTH 1000000
#define DICE_NFACE 6

int main(void)
{
	int dice,hist_dice[DICE_NFACE];
	int i,f;
	
	rand_timeinit();
	for (f=0;f<DICE_NFACE;f++)
	{
		hist_dice[f] = 0;
	}
	
	printf("-- Throwing a %d face(s) dice %d time(s)...\n",DICE_NFACE,SEQ_LENGTH);
	for (i=0;i<SEQ_LENGTH;i++)
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
		printf("face %d\t: %f%%\n",f+1,DRATIO(hist_dice[f],SEQ_LENGTH)*100.0);
	}
	
	return EXIT_SUCCESS;
}