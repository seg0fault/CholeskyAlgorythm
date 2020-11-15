#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fenv.h>
#include "matrix.h"
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/sysinfo.h>
#include <sched.h>
#include <unistd.h>

//n m r s file
int get_closest_divisor(int a, int b);

int get_closest_divisor(int a, int b)
{
    int guess1 = b;
    int guess2 = b;
    if(!(a % b))
        return b;
    for(; a % guess1; guess1--);
    for(; a % guess2; guess2++);
    return guess2 - b < b - guess1 ? guess2 : guess1;
}

int main(int argc, char *argv[])
{
    FILE *fp = nullptr;
    matrix A;
    matrix R;
    matrix V;
    matrix A_reverse;
    double *D;
    double t;
    double t1;
    int n;
    int m;
    int r;
    int s;
    int res;
    cpu_set_t cpu;
    int nprocs = get_nprocs();
    CPU_ZERO(&cpu);
    CPU_SET(nprocs - 1, &cpu);
    sched_setaffinity(getpid(), sizeof(cpu), &cpu);
    feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    if(argc < 5)
    {
        printf("Usage n m r s [filename]\n");
        return 1;
    }
    if(sscanf(argv[1], "%d", &n) != 1)
    {
        printf("Incorrect n\n");
        return 2;
    }
    if(sscanf(argv[2], "%d", &m) != 1)
    {
        printf("Incorrect m\n");
        return 2;
    }
    if(sscanf(argv[3], "%d", &r) != 1)
    {
        printf("Incorrect r\n");
        return 2;
    }
    if(sscanf(argv[4], "%d", &s) != 1)
    {
        printf("Incorrect s\n");
        return 2;
    }
    if(n % m)
    {
        printf("n mod m != 0 not supported yet :(\n");
        m = get_closest_divisor(n, m);
        printf("Warning! m redefined as %d\n", m);
    }
    if(r > n)
        r = n;
    if(!s)
    {
        if(argc >= 6)
        {
            fp = fopen(argv[5], "r");
            if(!fp)
            {
                printf("Cannot open file %s\n", argv[5]);
                return 4;
            }
        }
        else
        {
            printf("s = 0, but file not given\n");
            return 3;
        }
    }

    if((res = A.matrix_allocate(n, m)) != 0)
    {
        if(res < 0)
        {
            printf("Allocation A call error: %d\n", res);
            if(fp)
                fclose(fp);
            return 5;
        }
        else
        {
            printf("Not enough memory to allocate matrix A\n");
            if(fp)
                fclose(fp);
            return 6;
        }
    }
    if((res = R.matrix_allocate(n, m)) != 0)
    {
        if(res < 0)
        {
            printf("Allocation R call error: %d\n", res);
            if(fp)
                fclose(fp);
            return 5;
        }
        else
        {
            printf("Not enough memory to allocate matrix R\n");
            if(fp)
                fclose(fp);
            return 6;
        }
    }
    if((res = V.matrix_allocate(n, m)) != 0)
    {
        if(res < 0)
        {
            printf("Allocation R call error: %d\n", res);
            if(fp)
                fclose(fp);
            return 5;
        }
        else
        {
            printf("Not enough memory to allocate matrix R\n");
            if(fp)
                fclose(fp);
            return 6;
        }
    }
    if((res = A_reverse.matrix_allocate(n, m)) != 0)
    {
        if(res < 0)
        {
            printf("Allocation R call error: %d\n", res);
            if(fp)
                fclose(fp);
            return 5;
        }
        else
        {
            printf("Not enough memory to allocate matrix R\n");
            if(fp)
                fclose(fp);
            return 6;
        }
    }
    if(fp)
    {
        if((res = A.init_matrix_file(fp)))
        {
            printf("Read error code %d\n", res);
            if(fp)
                fclose(fp);
            return 6;
        }
    }
    else
    {
        switch(s)
        {
        case 1: A.init_matrix_formula(*function_1); break;
        case 2: A.init_matrix_formula(*function_2); break;
        case 3: A.init_matrix_formula(*function_3); break;
        case 4: A.init_matrix_formula(*function_4); break;
        default: printf("Unknown s\n"); return 7;
        }
    }
    D = (double*) calloc (n, sizeof (double));
    printf("Matrix initialized OK\n");
    printf("A:\n");
    A.print_matrix(r);
//    A.print_debug();
    printf("Matrix norm = %10.3e \n", A.get_norm());
    t = clock();
    if(A.get_cholesky(R, D, V))
    {
        printf("Bad matrix\n");
        if(fp)
            fclose(fp);
        free(D);
        return 1;
    }
    R.get_reversed_r(V);
    V.get_reversed_a(A_reverse, D);
    t = clock() - t;

    printf("R: \n");
    R.print_upper(r);
//    R.print_debug();
    printf("D: \n");
    for(res = 0; res < r; res++)
        printf("%10.3e ", D[res]);
    printf("\n\n");

    printf("R_reversed: \n");
    V.print_upper(r);

    printf("A_reversed: \n");
    A_reverse.print_matrix(r);
    t1 = clock();
    printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], R.residual(A, A_reverse), t/CLOCKS_PER_SEC, s, n, m);
    t1 = clock() - t;
    printf("residual time = %f\n", t1/CLOCKS_PER_SEC);
    if(fp)
        fclose(fp);
    free(D);
    return 0;
}






