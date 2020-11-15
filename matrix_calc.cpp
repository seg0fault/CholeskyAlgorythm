#include "matrix.h"

double
matrix::get_norm()
{
    int block_i;
    int block_j;
    int i;
    int j;
    double r0;
    double max = 0;
    double *it;
    double *ptr;
    int vert_step;
    for(block_i = 0; block_i < matrix_size_blocks; block_i++)
    {
        memset(block_vector, 0, block_size * sizeof (double));
        it = data + block_i * block_size * block_size;
        vert_step = (matrix_size_blocks - block_i - 1) * block_size * block_size;
        for (block_j = 0; block_j < block_i; block_j++) //идем вертикально
        {
            for(i = 0; i < block_size; i++, it++)
            {
                ptr = it;
                r0 = 0;
                for(j = 0; j < block_size; j++, ptr +=  block_size)
                {
                    r0 += fabs(*ptr);
//                    printf("[v]%10.3e \n", *ptr);
                }
                block_vector[i] += r0;
            }
            it = ptr - block_size + 1;
            it += vert_step + (block_i - block_j - 1) * block_size * block_size;
        }
        for (; block_j < matrix_size_blocks - 1; block_j++) //идем горизонтально до маленького блока
        {
            for(i = 0; i < block_size; i++)
            {
                r0 = 0;
                for(j = 0; j < block_size; j++, it++)
                {
                    r0 += fabs(*it);
//                    printf("[h]%10.3e \n", *it);
                }
                block_vector[i] += r0;
            }
        }
        for(i = 0; i < block_size; i++, it += block_size - r_block_size) //обрабатываем маленький блок
        {
            r0 = 0;
            for(j = 0; j < r_block_size; j++, it++)
            {
                r0 += fabs(*it);
//                printf("[l]%10.3e \n", *it);
            }
            block_vector[i] += r0;
            if(block_vector[i] > max)
                max = block_vector[i];
        }
    }
    epsilon = max * 1e-16;
    return max;
}

int matrix::get_cholesky(matrix& R, double *D, matrix& V)
{
    int i;
    int j;
    int k;
    double *R_ki;
    R_ii.allocate_block(block_size);
    for (i = 0; i < matrix_size_blocks; i++)
    {
        R_ii.get_block(get_block_address(i,i));
        for (k = 0; k < i; k++)
        {
            R_ki = R.get_block_address(k,i);
            R_ii.A_RtDR(R_ki, R_ki, D + block_size * k);
        }
        R.R_ii.link_block(R.get_block_address(i,i), block_size);
        if(R_ii.cholesky_block(R.R_ii, D + i * block_size, epsilon))
        {
            printf("ERROR! BLOCK (%d, %d) looks bad\n", i, i);
            return -1;
        }
        R.R_ii.reverse_block_r(R_ii);
        R_ii.put_block(V.get_block_address(i,i));
        for(j = i + 1; j < matrix_size_blocks; j++)
        {
            R.R_ij.link_block(R.get_block_address(i, j), block_size);
            R.R_ij.get_block(get_block_address(i, j));
            for (k = 0; k < i; k++)
            {
                R.R_ij.A_RtDR(R.get_block_address(k,i), R.get_block_address(k,j), D + block_size * k);
            }
            R.R_ij.DRA(R_ii, D + i * block_size);
        }
    }
    R_ii.free_block();
    R.R_ii.unlink_block();
    return 0;
}

void
matrix::get_reversed_r(matrix& V)
{
    int i;
    int j;
    int k;
    R_ii.allocate_block(block_size);
    for(i = 0; i < matrix_size_blocks; i++)
    {
        for(j = i + 1; j < matrix_size_blocks; j++)
        {
            R_ii.clear_block();
            R_ij.link_block(V.get_block_address(i, j), block_size);
            for(k = i; k < j; k++)
            {
                R_ii.minusVR(V.get_block_address(i,k), get_block_address(k, j));
            }
            R_ij.multiplyAB(R_ii.get_start(), V.get_block_address(j, j));
        }
    }
    R_ii.free_block();
}

void
matrix::get_reversed_a(matrix& Ar, double *D)
{
    int i;
    int j;
    int k;
    for(i = 0; i < matrix_size_blocks; i++)
    {
        for(j = i; j < matrix_size_blocks; j++)
        {
            R_ii.link_block(Ar.get_block_address(i, j), block_size);
            for(k = j; k < matrix_size_blocks; k++)
            {
                R_ii.VDVt(get_block_address(i, k), get_block_address(j, k), D + block_size * k);
            }
        }
    }
}

double
matrix::residual(matrix &A, matrix &Ar)
{
    int i;
    int j;
    int k;
    for (i = 0; i < matrix_size_blocks; i++)
    {
        for (j = i; j < matrix_size_blocks; j++)
        {
            R_ii.link_block(get_block_address(i, j), block_size);
            R_ii.clear_block();
            for(k = 0; k < matrix_size_blocks; k++)
            {
                if (i <= k)
                {
                    if (k <= j)
                        R_ii.plusAB(A.get_block_address(i, k), Ar.get_block_address(k, j));
                    else
                        R_ii.plusABt(A.get_block_address(i, k), Ar.get_block_address(k, j));
                }
                else if (k <= j)
                    R_ii.plusAtB(A.get_block_address(i, k), Ar.get_block_address(k, j));
                else
                    R_ii.plusAtBt(A.get_block_address(i, k), Ar.get_block_address(k, j));
            }
            if(j == i)
               R_ii.minusE();
        }
    }
    return get_norm();
}




















