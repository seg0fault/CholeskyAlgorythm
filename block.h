#ifndef BLOCK_H
#define BLOCK_H

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

class block
{
private:
    double *start = nullptr; //может указывать на место в matrix, а может иметь собственное выделение
    int size = 0;
    int is_independant = 0;
public:
    block() = default;
    void free_block();
    ~block(){free_block();}
    int allocate_block(int size_);
    void link_block(double *start_, int size_);
    void unlink_block();
    void get_block(double *src);
    double *get_start(){return start;}
    void put_block(double *dest);
    void print_block(); //Дебаг функция
    void clear_block();
    void reverse_block_r(block& dest); //Только верхнетреугольные блоки требуют обращения
    int cholesky_block(block& dest, double *block_vector, double epsilon); //Только верхнетреугольные блоки могут быть разложены
    void A_RtDR(double *Rt, double *R, double *D);
    void VDVt(double *Vt, double *V, double *D);
    void DRA(block& R, double *D);
    void minusVR(double *V, double *R);
    void multiplyAB(double *A, double *B);
    void minusE();
    void plusAB(double *A, double *B);
    void plusAtB(double *A, double *B);
    void plusABt(double *A, double *B);
    void plusAtBt(double *A, double *B);
};


#endif
