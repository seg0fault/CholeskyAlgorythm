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
    for(i = 0; i < size; i++)
    {
        r_ii = cur_start[i];
        for(k = 0, dest_k = dest.start; k < i; k++, dest_k += size)
            r_ii -= dest_k[i] * dest_k[i] * block_vector[k];
        if(r_ii < 0)
            d_ii = block_vector[i] = -1;
        else
            d_ii = block_vector[i] = 1;
        if((r_ii = cur_dest[i] = sqrt(fabs(r_ii))) < epsilon)
        {
            printf("%e ", sqrt(fabs(r_ii)));
            return -1;
        }
        for(j = i + 1; j < size; j++)
        {
            r_ij = cur_start[j];
//            printf("%lf / %lf\n", r_ij, r_ii);
            for(k = 0, dest_k = dest.start; k < i; k++, dest_k += size)
            {
                r_ij -= dest_k[i] * block_vector[k] * dest_k[j];
            }
            cur_dest[j] = r_ij / (r_ii * d_ii);
        }
        cur_start += size;
        cur_dest += size;
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
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += Rt[k * size + i] * D[k] * R[k * size + j];
            }
            start[i * size + j] -= sum;
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
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += V[i * size + k] * R[k * size + j];
            }
            start[i * size + j] -= sum;
        }
    }
}

void block::VDVt(double *Vt, double *V, double *D)
{
    int i;
    int j;
    int k;
    double sum;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += Vt[i * size + k] * D[k] * V[j * size + k];
            }
            start[i * size + j] += sum;
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
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += A[i * size + k] * B[k * size + j];
            }
            start[i * size + j] = sum;
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
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += A[i * size + k] * B[k * size + j];
            }
            start[i * size + j] += sum;
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
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += A[i * size + k] * B[j * size + k];
            }
            start[i * size + j] += sum;
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
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += A[k * size + i] * B[k * size + j];
            }
            start[i * size + j] += sum;
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
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += A[k * size + i] * B[j * size + k];
            }
            start[i * size + j] += sum;
        }
    }
}

void
block::minusE()
{
    double *s = start;
    for(int i = 0; i < size; i++, s += size + 1)
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
    double *dest_curr = dest.start;
    for (i = size - 1; i >= 0; i--)
    {
        r_ii = start[i * size + i];
        r_ii = 1 / r_ii;
        dest.start[i * size + i] = r_ii;
        for (j = 0; j < i; j++)
            dest.start[i * size + j] = 0;
        for (j = i + 1; j < size; j++)
        {
            sum = 0;
            for (k = i + 1; k <= j; k++)
            {
                sum += start[i * size + k] * dest.start[k * size + j];
            }
            dest.start[i * size + j] = -sum * r_ii;
        }
        dest_curr += size;
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
        for (j = 0; j < size; j++)
        {
            sum = 0;
            for (k = 0; k <= i; k++)
            {
                sum += R.start[k * size + i] * start[k * size + j];
            }
            start[i * size + j] = D[i] * sum;
        }
    }
}
