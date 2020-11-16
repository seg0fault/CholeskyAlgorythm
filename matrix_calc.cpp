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
                    r0 += fabs(*ptr);
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
                    r0 += fabs(*it);
                block_vector[i] += r0;
            }
        }
        for(i = 0; i < block_size; i++, it += block_size - r_block_size) //обрабатываем маленький блок
        {
            r0 = 0;
            for(j = 0; j < r_block_size; j++, it++)
                r0 += fabs(*it);
            block_vector[i] += r0;
            if(block_vector[i] > max)
                max = block_vector[i];
        }
    }
    epsilon = max * PRECISION;
    return max;
}

int matrix::get_cholesky(matrix& R, double *D, matrix& V)
{
    int i;
    int j;
    int k;
    double *R_ki;
    R_ii.allocate_block();
    for (i = 0; i < matrix_size_blocks; i++)
    {
        if(i != matrix_size_blocks - 1)
        {
            R_ii.get_block(get_block_address(i,i), block_size, block_size);
            R.R_ii.link_block(R.get_block_address(i,i), block_size, block_size);
        }
        else
        {
            R_ii.get_block(get_block_address(i,i), block_size, r_block_size);
            R.R_ii.link_block(R.get_block_address(i,i), block_size, r_block_size);
        }
        for (k = 0; k < i; k++)
        {
            R_ki = R.get_block_address(k,i);
            R_ii.A_RtDR(R_ki, R_ki, D + block_size * k);
        }
        if(R_ii.cholesky_block(R.R_ii, D + i * block_size, epsilon))
        {
            printf("ERROR! BLOCK (%d, %d) invertible\n", i, i);
            return -1;
        }
        R.R_ii.reverse_block_r(R_ii);
        R_ii.put_block(V.get_block_address(i,i));
        for(j = i + 1; j < matrix_size_blocks; j++)
        {
            if(j != matrix_size_blocks - 1)
            {
                R.R_ij.link_block(R.get_block_address(i, j), block_size, block_size);
                R.R_ij.get_block(get_block_address(i, j), block_size, block_size);
            }
            else
            {
                R.R_ij.link_block(R.get_block_address(i, j), block_size, r_block_size);
                R.R_ij.get_block(get_block_address(i, j), block_size, r_block_size);
            }
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
    R_ii.allocate_block();
    int size_i = block_size;
    int size_j;
    for (i = 0; i < matrix_size_blocks; i++)
    {
        size_j = block_size;
        if (i == matrix_size_blocks - 1)
            size_i = r_block_size;
        for (j = i + 1; j < matrix_size_blocks; j++)
        {
            if (j == matrix_size_blocks - 1)
                size_j = r_block_size;
            R_ii.clear_block();
            R_ii.size = size_i;
            R_ii.r_size = size_j;
            for (k = i; k < j; k++)
            {
                R_ii.minusVR(V.get_block_address(i,k), get_block_address(k, j));
            }
            R_ij.link_block(V.get_block_address(i, j), size_i, size_j);
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
    int size_i = block_size;
    int size_j;
    for(i = 0; i < matrix_size_blocks; i++)
    {
        size_j = block_size;
        if (i == matrix_size_blocks - 1)
            size_i = r_block_size;
        for(j = i; j < matrix_size_blocks; j++)
        {
            if (j == matrix_size_blocks - 1)
                size_j = r_block_size;
            R_ii.link_block(Ar.get_block_address(i, j), size_i, size_j);
            for(k = j; k < matrix_size_blocks; k++)
            {
                if(k != matrix_size_blocks - 1)
                    R_ii.VDVt(get_block_address(i, k), get_block_address(j, k), D + block_size * k, block_size);
                else
                    R_ii.VDVt(get_block_address(i, k), get_block_address(j, k), D + block_size * k, r_block_size);            }
        }
    }
}

double
matrix::residual(matrix &A, matrix &Ar)
{
    int i;
    int j;
    int k;
    int size_i = block_size;
    int size_j;
    for (i = 0; i < matrix_size_blocks; i++)
    {
        size_j = block_size;
        if (i == matrix_size_blocks - 1)
            size_i = r_block_size;
        for (j = i; j < matrix_size_blocks; j++)
        {
            if (j == matrix_size_blocks - 1)
                size_j = r_block_size;
            R_ii.link_block(get_block_address(i, j), size_i, size_j);
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




















