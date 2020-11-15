#include "block.h"

int
block::allocate_block()
{
    if(start)
    {
        printf("Trying to reallocate a block\n");
        return -1;
    }
    start = (double*) calloc ( size_it * size_it, sizeof (double));
    if(!start)
        return -2;
    is_independant = 1;
    return 0;
}

void
block::link_block(double *start_, int size_, int r_size_)
{
    if(is_independant)
    {
        printf("Trying to link allocated block!\n");
        return;
    }
    is_independant = 0;
    start = start_;
    size = size_;
    r_size = r_size_;
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
block::get_block(double *src, int size_, int r_size_)
{
    int i;
    int j;
    int index = 0;
    if(!start)
        return;
    r_size = r_size_;
    size = size_;
    for (i = 0; i < size; i++)
    {
        index = i * size_it;
        for (j = 0; j < r_size; j++, index++)
        {
            start[index] = src[index];
        }
    }
}

void
block::put_block(double *dest)
{
    int i;
    int j;
    int index;
    if(!start)
        return;
    for (i = 0; i < size; i++)
    {
        index = i * size_it;
        for (j = 0; j < r_size; j++, index++)
        {
            dest[index] = start[index];
        }
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
        for(j = 0; j < r_size; j++)
        {
//            printf("%10.3e ", ptr[j]);
            printf("%lf ", ptr[j]);
        }
        ptr += size_it;
        printf("\n");
    }
    printf("\n");
}

void
block::clear_block()
{
    memset(start, 0, size_it * size_it * sizeof (double));
}


















