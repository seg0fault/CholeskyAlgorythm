#include "block.h"


int block::cholesky_block(block& dest, double *block_vector, double epsilon)
{
    int i;
    int j;
    int k;
    double r_ii;
    double d_ii;
    double r_ij;
    double *cur_start = start;
    double *cur_dest = dest.start;
    double *dest_k;
    for(i = 0; i < r_size; i++)
    {
        r_ii = cur_start[i];
        for(k = 0, dest_k = dest.start; k < i; k++, dest_k += size_it)
            r_ii -= dest_k[i] * dest_k[i] * block_vector[k];
        if(r_ii < 0)
            d_ii = block_vector[i] = -1;
        else
            d_ii = block_vector[i] = 1;
        if((r_ii = cur_dest[i] = sqrt(fabs(r_ii))) < epsilon)
        {
//            printf("%e ", sqrt(fabs(r_ii)));
            return -1;
        }
        for(j = i + 1; j < r_size; j++)
        {
            r_ij = cur_start[j];
//            printf("%lf / %lf\n", r_ij, r_ii);
            for(k = 0, dest_k = dest.start; k < i; k++, dest_k += size_it)
            {
                r_ij -= dest_k[i] * block_vector[k] * dest_k[j];
            }
            cur_dest[j] = r_ij / (r_ii * d_ii);
        }
        cur_start += size_it;
        cur_dest += size_it;
    }
//    dest.print_block();
    return 0;
}

void
block::A_RtDR(double *Rt, double *R, double *D)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < r_size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += Rt[k * size_it + i] * D[k] * R[k * size_it + j];
            }
            start[i * size_it + j] -= sum;
        }
    }
}

void
block::minusVR(double *V, double *R)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < r_size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += V[i * size_it + k] * R[k * size_it + j];
            }
            start[i * size_it + j] -= sum;
        }
    }
}

void
block::multiplyAB(double *A, double *B)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < r_size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += A[i * size_it + k] * B[k * size_it + j];
            }
            start[i * size_it + j] = sum;
        }
    }
}

void block::VDVt(double *Vt, double *V, double *D, int size_small)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < r_size; j++)
        {
            sum = 0;
            for (k = 0; k < size_small; k++)
            {
                sum += Vt[i * size_it + k] * D[k] * V[j * size_it + k];
            }
            start[i * size_it + j] += sum;
        }
    }
}

void
block::plusAB(double *A, double *B)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size_it; i++)
    {
        for (j = 0; j < size_it; j++)
        {
            sum = 0;
            for (k = 0; k < size_it; k++)
            {
                sum += A[i * size_it + k] * B[k * size_it + j];
            }
            start[i * size_it + j] += sum;
        }
    }
}

void
block::plusABt(double *A, double *B)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size_it; i++)
    {
        for (j = 0; j < size_it; j++)
        {
            sum = 0;
            for (k = 0; k < size_it; k++)
            {
                sum += A[i * size_it + k] * B[j * size_it + k];
            }
            start[i * size_it + j] += sum;
        }
    }
}

void
block::plusAtB(double *A, double *B)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size_it; i++)
    {
        for (j = 0; j < size_it; j++)
        {
            sum = 0;
            for (k = 0; k < size_it; k++)
            {
                sum += A[k * size_it + i] * B[k * size_it + j];
            }
            start[i * size_it + j] += sum;
        }
    }
}

void
block::plusAtBt(double *A, double *B)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size_it; i++)
    {
        for (j = 0; j < size_it; j++)
        {
            sum = 0;
            for (k = 0; k < size_it; k++)
            {
                sum += A[k * size_it + i] * B[j * size_it + k];
            }
            start[i * size_it + j] += sum;
        }
    }
}

void
block::minusE()
{
    double *s = start;
    for(int i = 0; i < size; i++, s += size_it + 1)
        *s = *s - 1;
}

void
block::reverse_block_r(block& dest)
{
    int i;
    int j;
    int k;
    double r_ii;
    double sum;
    for (i = r_size - 1; i >= 0; i--)
    {
        r_ii = start[i * size_it + i];
        r_ii = 1 / r_ii;
        dest.start[i * size_it + i] = r_ii;
        for (j = 0; j < i; j++)
            dest.start[i * size_it + j] = 0;
        for (j = i + 1; j < r_size; j++)
        {
            sum = 0;
            for (k = i + 1; k <= j; k++)
            {
                sum += start[i * size_it + k] * dest.start[k * size_it + j];
            }
            dest.start[i * size_it + j] = -sum * r_ii;
        }
    }
}

void
block::DRA(block& R, double *D)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = size - 1; i >= 0; i--)
    {
        for (j = 0; j < r_size; j++)
        {
            sum = 0;
            for (k = 0; k <= i; k++)
            {
                sum += R.start[k * size_it + i] * start[k * size_it + j];
            }
            start[i * size_it + j] = D[i] * sum;
        }
    }
}

/*
void
block::DRA(block& R, double *D)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < r_size; j++)
        {
            sum = 0;
            for (k = i; k <= size; k++)
            {
                sum += R.start[i * size_it + k] * start[k * size_it + j];
            }
            start[i * size_it + j] = D[i] * sum;
        }
    }
}
*/
