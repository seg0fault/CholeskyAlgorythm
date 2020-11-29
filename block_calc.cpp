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
        if(fabs(r_ii) < epsilon)
            return -1;
        cur_dest[i] = r_ii = sqrt(fabs(r_ii));
        for(j = i + 1; j < r_size; j++)
        {
            r_ij = cur_start[j];
            for(k = 0, dest_k = dest.start; k < i; k++, dest_k += size_it)
            {
                r_ij -= dest_k[i] * block_vector[k] * dest_k[j];
            }
            cur_dest[j] = r_ij / (r_ii * d_ii);
        }
        cur_start += size_it;
        cur_dest += size_it;
    }
    return 0;
}

void
block::A_RtDR(double *Rt, double *R, double *D)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;
  double Dk;
  double *ptrRt;
  double *ptrR;
  double *ptrWrite;
  for (i = 0, ptrWrite = start; i < size; i += 3, ptrWrite += 3 * size_it)
    {
      for (j = 0; j + 2 < r_size; j += 3)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;

          for (k = 0, ptrRt = Rt, ptrR = R; k < size; k ++, ptrRt += size_it, ptrR += size_it)
            {
              Dk = D[k];
              sum0 += ptrRt[i] * Dk * ptrR[j];
              sum1 += ptrRt[i + 1] * Dk * ptrR[j];
              sum2 += ptrRt[i + 2] * Dk * ptrR[j];

              sum3 += ptrRt[i] * Dk * ptrR[j + 1];
              sum4 += ptrRt[i + 1] * Dk * ptrR[j + 1];
              sum5 += ptrRt[i + 2] * Dk * ptrR[j + 1];

              sum6 += ptrRt[i] * Dk * ptrR[j + 2];
              sum7 += ptrRt[i + 1] * Dk * ptrR[j + 2];
              sum8 += ptrRt[i + 2] * Dk * ptrR[j + 2];
            }

          ptrWrite[j] -= sum0;
          ptrWrite[size_it + j] -= sum1;
          ptrWrite[2 * size_it + j] -= sum2;

          ptrWrite[j + 1] -= sum3;
          ptrWrite[size_it + j + 1] -= sum4;
          ptrWrite[2 * size_it + j + 1] -= sum5;

          ptrWrite[j + 2] -= sum6;
          ptrWrite[size_it + j + 2] -= sum7;
          ptrWrite[2 * size_it + j + 2] -= sum8;
        }
      for (; j < r_size; j++)
        {
          sum0 = sum1 = sum2 = 0.;

          for (k = 0, ptrRt = Rt, ptrR = R; k < size; k ++, ptrRt += size_it, ptrR += size_it)
            {
              Dk = D[k];
              sum0 += ptrRt[i] * Dk * ptrR[j];
              sum1 += ptrRt[i + 1] * Dk * ptrR[j];
              sum2 += ptrRt[i + 2] * Dk * ptrR[j];
            }

          ptrWrite[j] -= sum0;
          ptrWrite[size_it + j] -= sum1;
          ptrWrite[2 * size_it + j] -= sum2;
        }
    }
}

void
block::minusVR(double *V, double *R)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;

  double *ptrV, *ptrV1, *ptrV2;
  double *ptrR;
  double *ptrWrite;

  for (i = 0, ptrWrite = start, ptrV = V, ptrV1 = V + size_it, ptrV2 = V + 2 * size_it;
       i < size;
       i += 3, ptrWrite += 3 * size_it, ptrV += 3 * size_it, ptrV1 += 3 * size_it, ptrV2 += 3 * size_it)
    {
      for (j = 0; j + 2 < r_size; j += 3)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;

          for (k = 0, ptrR = R; k < size; k++, ptrR += size_it)
            {
              sum0 += ptrV[k] * ptrR[j];
              sum1 += ptrV1[k] * ptrR[j];
              sum2 += ptrV2[k] * ptrR[j];

              sum3 += ptrV[k] * ptrR[j + 1];
              sum4 += ptrV1[k] * ptrR[j + 1];
              sum5 += ptrV2[k] * ptrR[j + 1];

              sum6 += ptrV[k] * ptrR[j + 2];
              sum7 += ptrV1[k] * ptrR[j + 2];
              sum8 += ptrV2[k] * ptrR[j + 2];
            }

          ptrWrite[j] -= sum0;
          ptrWrite[size_it + j] -= sum1;
          ptrWrite[2 * size_it + j] -= sum2;

          ptrWrite[j + 1] -= sum3;
          ptrWrite[size_it + j + 1] -= sum4;
          ptrWrite[2 * size_it + j + 1] -= sum5;

          ptrWrite[j + 2] -= sum6;
          ptrWrite[size_it + j + 2] -= sum7;
          ptrWrite[2 * size_it + j + 2] -= sum8;
        }

      for (; j < r_size; j++)
        {
          sum0 = sum1 = sum2 = 0.;

          for (k = 0, ptrR = R; k < size; k++, ptrR += size_it)
            {
              sum0 += ptrV[k] * ptrR[j];
              sum1 += ptrV1[k] * ptrR[j];
              sum2 += ptrV2[k] * ptrR[j];
            }

          ptrWrite[j] -= sum0;
          ptrWrite[size_it + j] -= sum1;
          ptrWrite[2 * size_it + j] -= sum2;
        }
    }
}

void
block::multiplyAB(double *A, double *B)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;
  double *ptrA, *ptrA1, *ptrA2;
  double *ptrB;
  double *ptrWrite;
  for (i = 0, ptrWrite = start, ptrA = A, ptrA1 = A + size_it, ptrA2 = A + 2 * size_it;
       i < size;
       i += 3, ptrWrite += 3 * size_it, ptrA += 3 * size_it, ptrA1 += 3 * size_it, ptrA2 += 3 * size_it)
    {
      for (j = 0; j + 2 < r_size; j += 3)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;

          for (k = 0, ptrB = B; k < size; k++, ptrB += size_it)
            {
              sum0 += ptrA[k] * ptrB[j];
              sum1 += ptrA1[k] * ptrB[j];
              sum2 += ptrA2[k] * ptrB[j];

              sum3 += ptrA[k] * ptrB[j + 1];
              sum4 += ptrA1[k] * ptrB[j + 1];
              sum5 += ptrA2[k] * ptrB[j + 1];

              sum6 += ptrA[k] * ptrB[j + 2];
              sum7 += ptrA1[k] * ptrB[j + 2];
              sum8 += ptrA2[k] * ptrB[j + 2];
            }

          ptrWrite[j] = sum0;
          ptrWrite[size_it + j] = sum1;
          ptrWrite[2 * size_it + j] = sum2;

          ptrWrite[j + 1] = sum3;
          ptrWrite[size_it + j + 1] = sum4;
          ptrWrite[2 * size_it + j + 1] = sum5;

          ptrWrite[j + 2] = sum6;
          ptrWrite[size_it + j + 2] = sum7;
          ptrWrite[2 * size_it + j + 2] = sum8;
        }
      for (; j < r_size; j++)
        {
          sum0 = sum1 = sum2 = 0.;

          for (k = 0, ptrB = B; k < size; k++, ptrB += size_it)
            {
              sum0 += ptrA[k] * ptrB[j];
              sum1 += ptrA1[k] * ptrB[j];
              sum2 += ptrA2[k] * ptrB[j];
            }

          ptrWrite[j] = sum0;
          ptrWrite[size_it + j] = sum1;
          ptrWrite[2 * size_it + j] = sum2;
        }
    }
}

void block::VDVt(double *Vt, double *V, double *D, int size_small)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;
  double Dk;
  double *ptrVt, *ptrVt1, *ptrVt2;
  double *ptrV, *ptrV1, *ptrV2;
  double *ptrWrite;
  for (i = 0, ptrWrite = start, ptrVt = Vt, ptrVt1 = Vt + size_it, ptrVt2 = Vt + 2 * size_it;
       i < size;
       i += 3, ptrWrite += 3 * size_it, ptrVt += 3 * size_it, ptrVt1 += 3 * size_it, ptrVt2 += 3 * size_it)
    {
      for (j = 0, ptrV = V, ptrV1 = V + size_it, ptrV2 = V + 2 * size_it;
           j + 2 < r_size;
           j += 3, ptrV += 3 * size_it, ptrV1 += 3 * size_it, ptrV2 += 3 * size_it)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;

          for (k = 0; k < size_small; k++)
            {
              Dk = D[k];
              sum0 += ptrVt[k] * Dk * ptrV[k];
              sum1 += ptrVt1[k] * Dk * ptrV[k];
              sum2 += ptrVt2[k] * Dk * ptrV[k];

              sum3 += ptrVt[k] * Dk * ptrV1[k];
              sum4 += ptrVt1[k] * Dk * ptrV1[k];
              sum5 += ptrVt2[k] * Dk * ptrV1[k];

              sum6 += ptrVt[k] * Dk * ptrV2[k];
              sum7 += ptrVt1[k] * Dk * ptrV2[k];
              sum8 += ptrVt2[k] * Dk * ptrV2[k];
            }

          ptrWrite[j] += sum0;
          ptrWrite[size_it + j] += sum1;
          ptrWrite[2 * size_it + j] += sum2;

          ptrWrite[j + 1] += sum3;
          ptrWrite[size_it + j + 1] += sum4;
          ptrWrite[2 * size_it + j + 1] += sum5;

          ptrWrite[j + 2] += sum6;
          ptrWrite[size_it + j + 2] += sum7;
          ptrWrite[2 * size_it + j + 2] += sum8;
        }

      for (; j < r_size; j++, ptrV += size_it)
        {
          sum0 = sum1 = sum2 = 0.;

          for (k = 0; k < size_small; k++)
            {
              Dk = D[k];
              sum0 += ptrVt[k] * Dk * ptrV[k];
              sum1 += ptrVt1[k] * Dk * ptrV[k];
              sum2 += ptrVt2[k] * Dk * ptrV[k];
            }

          ptrWrite[j] += sum0;
          ptrWrite[size_it + j] += sum1;
          ptrWrite[2 * size_it + j] += sum2;
        }
    }
}

void
block::plusAB(double *A, double *B)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;
  double *ptrA, *ptrA1, *ptrA2;
  double *ptrB;
  double *ptrWrite;
  for (i = 0, ptrWrite = start, ptrA = A, ptrA1 = A + size_it, ptrA2 = A + 2 * size_it;
       i < size_it;
       i += 3, ptrWrite += 3 * size_it, ptrA += 3 * size_it, ptrA1 += 3 * size_it, ptrA2 += 3 * size_it)
    {
      for (j = 0; j < size_it; j += 3)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;

          for (k = 0, ptrB = B; k < size_it; k++, ptrB += size_it)
            {
              sum0 += ptrA[k] * ptrB[j];
              sum1 += ptrA1[k] * ptrB[j];
              sum2 += ptrA2[k] * ptrB[j];

              sum3 += ptrA[k] * ptrB[j + 1];
              sum4 += ptrA1[k] * ptrB[j + 1];
              sum5 += ptrA2[k] * ptrB[j + 1];

              sum6 += ptrA[k] * ptrB[j + 2];
              sum7 += ptrA1[k] * ptrB[j + 2];
              sum8 += ptrA2[k] * ptrB[j + 2];
            }

          ptrWrite[j] += sum0;
          ptrWrite[size_it + j] += sum1;
          ptrWrite[2 * size_it + j] += sum2;

          ptrWrite[j + 1] += sum3;
          ptrWrite[size_it + j + 1] += sum4;
          ptrWrite[2 * size_it + j + 1] += sum5;

          ptrWrite[j + 2] += sum6;
          ptrWrite[size_it + j + 2] += sum7;
          ptrWrite[2 * size_it + j + 2] += sum8;
        }
    }
}

void
block::plusABt(double *A, double *B)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;
  double *ptrA, *ptrA1, *ptrA2;
  double *ptrB, *ptrB1, *ptrB2;
  double *ptrWrite;
  for (i = 0, ptrWrite = start, ptrA = A, ptrA1 = A + size_it, ptrA2 = A + 2 * size_it;
       i < size_it;
       i += 3, ptrWrite += 3 * size_it, ptrA += 3 * size_it, ptrA1 += 3 * size_it, ptrA2 += 3 * size_it)
    {
      for (j = 0, ptrB = B, ptrB1 = B + size_it, ptrB2 = B + 2 * size_it;
           j < size_it;
           j += 3, ptrB += 3 * size_it, ptrB1 += 3 * size_it, ptrB2 += 3 * size_it)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;
          for (k = 0; k < size_it; k++)
            {
              sum0 += ptrA[k] * ptrB[k];
              sum1 += ptrA1[k] * ptrB[k];
              sum2 += ptrA2[k] * ptrB[k];

              sum3 += ptrA[k] * ptrB1[k];
              sum4 += ptrA1[k] * ptrB1[k];
              sum5 += ptrA2[k] * ptrB1[k];

              sum6 += ptrA[k] * ptrB2[k];
              sum7 += ptrA1[k] * ptrB2[k];
              sum8 += ptrA2[k] * ptrB2[k];
            }

          ptrWrite[j] += sum0;
          ptrWrite[size_it + j] += sum1;
          ptrWrite[2 * size_it + j] += sum2;

          ptrWrite[j + 1] += sum3;
          ptrWrite[size_it + j + 1] += sum4;
          ptrWrite[2 * size_it + j + 1] += sum5;

          ptrWrite[j + 2] += sum6;
          ptrWrite[size_it + j + 2] += sum7;
          ptrWrite[2 * size_it + j + 2] += sum8;
        }
    }
}

void
block::plusAtB(double *A, double *B)
{
  int i;
  int j;
  int k;
  double sum0, sum3, sum6;
  double sum1, sum4, sum7;
  double sum2, sum5, sum8;
  double *ptrA;
  double *ptrB;
  double *ptrWrite;
  for (i = 0, ptrWrite = start; i < size_it; i += 3, ptrWrite += 3 * size_it)
    {
      for (j = 0; j < size_it; j += 3)
        {
          sum0 = sum3 = sum6 = 0.;
          sum1 = sum4 = sum7 = 0.;
          sum2 = sum5 = sum8 = 0.;

          for (k = 0, ptrA = A, ptrB = B; k < size_it; k++, ptrA += size_it, ptrB += size_it)
            {
              sum0 += ptrA[i] * ptrB[j];
              sum1 += ptrA[i + 1] * ptrB[j];
              sum2 += ptrA[i + 2] * ptrB[j];

              sum3 += ptrA[i] * ptrB[j + 1];
              sum4 += ptrA[i + 1] * ptrB[j + 1];
              sum5 += ptrA[i + 2] * ptrB[j + 1];

              sum6 += ptrA[i] * ptrB[j + 2];
              sum7 += ptrA[i + 1] * ptrB[j + 2];
              sum8 += ptrA[i + 2] * ptrB[j + 2];
            }

          ptrWrite[j] += sum0;
          ptrWrite[size_it + j] += sum1;
          ptrWrite[2 * size_it + j] += sum2;

          ptrWrite[j + 1] += sum3;
          ptrWrite[size_it + j + 1] += sum4;
          ptrWrite[2 * size_it + j + 1] += sum5;

          ptrWrite[j + 2] += sum6;
          ptrWrite[size_it + j + 2] += sum7;
          ptrWrite[2 * size_it + j + 2] += sum8;
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
    double sum1;
    double sum2;
    for (i = 0; i < size_it; i++)
    {
        for (j = 0; j < size_it; j++)
        {
            sum = 0;
            sum1 = 0;
            sum2 = 0;
            for (k = 0; k < size_it; k += 3)
            {
                sum += A[k * size_it + i] * B[j * size_it + k];
                sum1 += A[(k + 1) * size_it + i] * B[j * size_it + k + 1];
                sum2 += A[(k + 2) * size_it + i] * B[j * size_it + k + 2];
            }
            start[i * size_it + j] += (sum + sum1 + sum2);
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
