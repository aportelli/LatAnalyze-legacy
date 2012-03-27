#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <latan/latan_math.h>

#define MAX_MAT_SIZE 200
#define MAT_SIZE_STEP 10
#define TEST_TIME 2

static latan_errno mat_mul_nn(mat *a, const mat *b, const mat *c)
{
    return mat_mul(a,b,'n',c,'n');
}

#define NFLOP_MAT_EQADD(row_b,col_b)\
(2.0*(double)(row_b)*(double)(col_b))
#define NFLOP_MAT_MUL_NN(row_b,col_c,col_b)\
(2.0*(double)(row_b)*(double)(col_c)*(double)(col_b))

#define FUNC_BENCH_2ARG(perf,nops,func,func_flop,a,b)\
{\
    clock_t end_clk,timer;\
    \
    nops    = 0;\
    \
    end_clk = clock() + TEST_TIME*CLOCKS_PER_SEC;\
    timer   = clock();\
    while (clock() < end_clk)\
    {\
        (func)(a,b);\
        nops++;\
    }\
    nops  /= (double)(timer)/(double)(CLOCKS_PER_SEC);\
    perf   = func_flop*nops/(1.0e9);\
}
#define FUNC_BENCH_3ARG(perf,nops,func,func_flop,a,b,c)\
{\
    clock_t end_clk,timer;\
    \
    nops    = 0;\
    \
    end_clk = clock() + TEST_TIME*CLOCKS_PER_SEC;\
    timer   = clock();\
    while (clock() < end_clk)\
    {\
        (func)(a,b,c);\
        nops++;\
    }\
    timer  = clock() - timer;\
    nops  /= (double)(timer)/(double)(CLOCKS_PER_SEC);\
    perf   = func_flop*nops/(1.0e9);\
}

int main(void)
{
    strbuf name,version;
    size_t mat_size;
    double func_flop,perf,nops;
    mat *a,*b,*c;
    
    latan_get_name(name);
    latan_get_version(version);
    printf("\n");
    printf("%s v%s matrix operations benchmark\n",name,version);
    printf("\n");
    printf("test time : %d sec\n",TEST_TIME);
    printf("\n");
    
    printf("mat_mul_nn test :\n");
    printf("-----------------\n");
    printf("general matrix\n");
    printf("%6s %10s %10s\n","size","Gflop/s","Nop/s");
    for (mat_size=MAT_SIZE_STEP;mat_size<=MAX_MAT_SIZE;mat_size+=MAT_SIZE_STEP)
    {
        func_flop = NFLOP_MAT_MUL_NN(mat_size,mat_size,mat_size);
        
        a = mat_create(mat_size,mat_size);
        b = mat_create_from_dim(a);
        c = mat_create_from_dim(a);
        
        mat_rand_u(b,0.0,1.0);
        mat_rand_u(c,0.0,1.0);
        
        FUNC_BENCH_3ARG(perf,nops,mat_mul_nn,func_flop,a,b,c);
        printf("%6d %10f %10.1f\n",(int)mat_size,perf,nops);
        
        mat_destroy(a);
        mat_destroy(b);
        mat_destroy(c);
    }
    printf("\n");
    
    printf("mat_eqadd test :   \n");
    printf("-----------------\n");
    printf("%6s %10s\n","size","Gflop/s");
    for (mat_size=MAT_SIZE_STEP;mat_size<=MAX_MAT_SIZE;mat_size+=MAT_SIZE_STEP)
    {
        func_flop = NFLOP_MAT_EQADD(mat_size,mat_size);
        
        a = mat_create(mat_size,mat_size);
        b = mat_create_from_dim(a);

        mat_rand_u(a,0.0,1.0);
        mat_rand_u(b,0.0,1.0);
        
        FUNC_BENCH_2ARG(perf,nops,mat_eqadd,func_flop,a,b);
        printf("%6d %10f %10.1f\n",(int)mat_size,perf,nops);
        
        mat_destroy(a);
        mat_destroy(b);
    }
    printf("\n");
    
    return EXIT_SUCCESS;
}
