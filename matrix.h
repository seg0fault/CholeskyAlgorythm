#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "block.h"
//void *memset(void *s, int c, size_t n);


//Функции автозаполнения
double function_1(int n, int i, int j);
double function_2(int n, int i, int j);
double function_3(int n, int i, int j);
double function_4(int n, int i, int j);

class matrix
{
private:
    double *data = nullptr;
    double *block_vector = nullptr;
    double epsilon = 0;
    block R_ii;
    block R_ij;
    int matrix_size = 0;
    int matrix_size_blocks = 0;
    int block_size = 0;
    int r_block_size = 0;
    int blocks_count = 0;
    int elements_count = 0;
public:
    matrix() = default;
    ~matrix() {matrix_free();}
    int matrix_allocate(int n, int m);
    void matrix_free(){if(data){free(data); free(block_vector);}  data = nullptr; block_vector = nullptr;}

    int init_matrix_file(FILE *fp);
    void init_matrix_formula(double (*f)(int, int, int));

    void print_matrix(int max_print);
    void print_upper(int max_print);
    void print_debug();
    double *get_block_address(int block_i, int block_j);
    double *get_element_address(int i, int j);
    double get_element(int i, int j);
    void put_element(double element, int i, int j);

    double get_norm();

    void get_block_vector(double *src);
    void put_block_vector(double *dest);

    int get_cholesky(matrix& R, double *D, matrix& V); //будем сразу дозаписывать обращенные блоки в V
    void get_reversed_r(matrix& V);
    void get_reversed_a(matrix& Ar, double *D);
    double residual(matrix &A, matrix &Ar);
};


#endif
