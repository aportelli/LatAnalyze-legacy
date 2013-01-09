/* ex_endian.c, part of LatAnalyze library
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
#include <stdlib.h>
#include <latan/latan_globals.h>

int main(void)
{
    FILE *f;
    double dtab[4];
    int itab[4];
    int i;
    
    f = fopen("ex_endian_be.bin","r");
    fread(dtab,sizeof(double),4,f);
    fread(itab,sizeof(int),4,f);
    fclose(f);
    for (i=0;i<4;i++)
    {
        dtab[i] = latan_conv_endianness_d(dtab[i],BE);
        printf("dtab[%d] = %e\n",i,dtab[i]);
    }
    for (i=0;i<4;i++)
    {
        itab[i] = latan_conv_endianness_i(itab[i],BE);
        printf("itab[%d] = %d\n",i,itab[i]);
    }
    f = fopen("ex_endian_le.bin","w");
    fwrite(dtab,sizeof(double),4,f);
    fwrite(itab,sizeof(int),4,f);
    fclose(f);
    
    return EXIT_SUCCESS;
}
