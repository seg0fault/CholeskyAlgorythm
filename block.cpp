#include "block.h"

int
block::allocate_block(int size_)
{
    if(start)
    {
        printf("Trying to reallocate a block\n");
        return -1;
    }
    start = (double*) calloc ( size_ * size_, sizeof (double));
    if(!start)
        return -2;
    size = size_;
    is_independant = 1;
    return 0;
}

void
block::link_block(double *start_, int size_)
{
    if(is_independant)
    {
        printf("Trying to link allocated block!\n");
        return;
    }
    is_independant = 0;
    start = start_;
    size = size_;
}

void
block::unlink_block()
{
    start = nullptr;
    size = 0;
}

void
block::free_block()
{
    if(start && is_independant)
    {
        free(start);
        start = nullptr;
        is_independant = 0;
    }
    /*
    else
    {
        printf("Trying to bad free\n");
    }
    */
}

void
block::get_block(double *src)
{
    int i;
    if(!start)
        return;
    for (i = 0; i < size * size; i++)
    {
        start[i] = src[i];
    }
}

void
block::put_block(double *dest)
{
    int i;
    if(!start)
        return;
    for (i = 0; i < size * size; i++)
    {
        dest[i] = start[i];
    }
}

void
block::print_block()
{
    int i;
    int j;
    double *ptr = start;
    if(!start)
    {
        printf("Block not allocated!\n");
        return;
    }
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
//            printf("%10.3e ", ptr[j]);
            printf("%lf ", ptr[j]);
        }
        ptr += size;
        printf("\n");
    }
    printf("\n");
}

void
block::clear_block()
{
    memset(start, 0, size * size * sizeof (double));
}


















