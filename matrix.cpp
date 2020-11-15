#include "matrix.h"

double function_1(int n, int i, int j){return n - (i < j ? j : i);}
double function_2(int n, int i, int j){(void)n; return (i < j ? j : i);}
double function_3(int n, int i, int j){(void)n; return fabs(i - j);}
double function_4(int n, int i, int j){(void)n; return 1.0/(double)(i + j + 1);}

int
matrix::matrix_allocate(int n, int m)
{
    if (data)
        return -4;
    if (n <= 0 || m <= 0)
        return -1;
    if (m > n)
        return -2;
    matrix_size = n;
    block_size = m;
    r_block_size = n % m;
    R_ii.size_it = block_size;
    R_ij.size_it = block_size;
    if (r_block_size)
        matrix_size_blocks = n / m + 1;
    else
      {
        matrix_size_blocks = n / m;
        r_block_size = block_size;
      }
    blocks_count = matrix_size_blocks * (matrix_size_blocks + 1) / 2;
    elements_count = blocks_count * block_size * block_size;
    data = (double*) calloc (elements_count, sizeof (double));
    block_vector = (double*) malloc (block_size * sizeof (double));
    if(!data)
        return 1;
    return 0;
}

void
matrix::init_matrix_formula(double (*f)(int, int, int))
{
    int i;
    int j;
    int offset_i = 0;
    int offset_j = 0;
    int block_i;
    int block_j;
    int index = 0;
    for(block_i = 0; block_i < matrix_size_blocks - 1; block_i++)
    {
        offset_j = offset_i;
        for(block_j = block_i; block_j < matrix_size_blocks - 1; block_j++)
        {
            for(i = 0; i < block_size; i++)
            {
                for(j = 0; j < block_size; j++)
                {
                    data[index] = f(matrix_size, offset_i + i, offset_j + j);
                    index++;
                }
            }
            offset_j += block_size;
        }
        for(i = 0; i < block_size; i++)
        {
            for(j = 0; j < r_block_size; j++)
            {
                data[index] = f(matrix_size, offset_i + i, offset_j + j);
                index++;
            }
            index += (block_size - r_block_size);
        }
//        index += (block_size - r_block_size) * block_size; //???
        offset_i += block_size;
    }
    for(i = 0; i < r_block_size; i++)
    {
        for(j = 0; j < r_block_size; j++)
        {
            data[index] = f(matrix_size, offset_i + i, offset_j + j);
            index++;
        }
        index += (block_size - r_block_size);
    }
}

int
matrix::init_matrix_file(FILE *fp)
{
    int i;
    int j;
    int k;
    double tmp;
    double *block;
    for (i = 0; i < matrix_size; i++)
    {
        for (j = 0; j < i; j++)
        {
            if(fscanf(fp, "%lf", &tmp) != 1)
                return 1;
        }
        for (; j < matrix_size; j++)
        {
            if(fscanf(fp, "%lf", &tmp) != 1)
                return 1;
            *(get_element_address(i, j)) = tmp;
        }
    }
    //cure
    for(i = 0; i < matrix_size_blocks; i++)
    {
        block = get_block_address(i, i);
        for(j = 0; j < block_size; j++)
        {
            for(k = j + 1; k < block_size; k++)
            {
                block[k * block_size + j] = block[j * block_size + k];
            }
        }
    }
    return 0;
}

double matrix::get_element(int i, int j)
{
    int block_i;
    int block_j;
    int inside_i;
    int inside_j;
    int index;
    block_i = i / block_size;
    block_j = j / block_size;
    inside_i = i % block_size;
    inside_j = j % block_size;
    index = (matrix_size_blocks * block_i - block_i * (block_i + 1) / 2 + block_j) * block_size * block_size + inside_i * block_size + inside_j;
    return data[index];
}

double* matrix::get_element_address(int i, int j)
{
    int block_i;
    int block_j;
    int inside_i;
    int inside_j;
    int index;
    block_i = i / block_size;
    block_j = j / block_size;
    inside_i = i % block_size;
    inside_j = j % block_size;
    index = (matrix_size_blocks * block_i - block_i * (block_i + 1) / 2 + block_j) * block_size * block_size + inside_i * block_size + inside_j;
    return data + index;
}

double* matrix::get_block_address(int block_i, int block_j)
{
    int tmp;
    if(block_j < block_i)
    {
        tmp = block_i;
        block_i = block_j;
        block_j = tmp;
    }
    return data + (matrix_size_blocks * block_i - block_i * (block_i + 1) / 2 + block_j) * block_size * block_size;
}

void
matrix::print_matrix(int max_print)
{
    int i;
    int j;
    if(max_print > matrix_size)
        max_print = matrix_size;
    for(i = 0; i < max_print; i++)
    {
        for(j = 0; j < i; j++)
            printf("%10.3e ", get_element(j, i));
        for(; j < max_print; j++)
            printf("%10.3e ", get_element(i, j));
        printf("\n");
    }
    printf("\n");
}

void
matrix::print_upper(int max_print)
{
    int i;
    int j;
    if(max_print > matrix_size)
        max_print = matrix_size;
    for(i = 0; i < max_print; i++)
    {
        for(j = 0; j < i; j++)
            printf("%10.3e ", 0.);
        for(; j < max_print; j++)
            printf("%10.3e ", get_element(i, j));
        printf("\n");
    }
    printf("\n");
}

void
matrix::print_debug()
{
    for(int i = 0; i < elements_count; i++)
    {
        printf("%10.3e ", data[i]);
    }
    printf("\n");
}
























