/* ex_io.c, part of LatAnalyze library
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
#include <latan/latan_mat.h>
#include <latan/latan_rand.h>
#include <latan/latan_io.h>

#define FNAME "ex_io.dat"

int main(void)
{
    mat *m;
    size_t dim[2];
    
    m = mat_create(6,6);
    
    io_init();
    
    io_set_fmt(IO_ASCII);
    mat_load(m,dim,"ex_io.dat:m");
    printf("%dx%d matrix loaded from %s:\n",(int)(dim[0]),(int)(dim[1]),FNAME);
    mat_print(m,"% .15e");
    io_set_fmt(IO_XML);
    mat_save("ex_io.xml:m",'w',m);

    io_finish();
    
    mat_destroy(m);
    
    return EXIT_SUCCESS;
}
