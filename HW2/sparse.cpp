#include "mmreader.hpp"
#include <time.h>
#include <iostream>
#include <sys/time.h>
#include <unistd.h>

bool
SCsrMatrixfromFile(struct sparse_mtx *A, const char* filePath)
{
    // Check that the file format is matrix market; the only format we can read right now
    // This is not a complete solution, and fails for directories with file names etc...
    // TODO: Should we use boost filesystem?
    std::string strPath( filePath );
    if( strPath.find_last_of( '.' ) != std::string::npos )
    {
        std::string ext = strPath.substr( strPath.find_last_of( '.' ) + 1 );
        if( ext != "mtx" )
        {
            std::cout << "Reading file name error" << std::endl;
            return false;
        }
    }
    else
        return false;

    // Read data from a file on disk into buffers
    // Data is read natively as COO format with the reader
    MatrixMarketReader mm_reader;
    if( mm_reader.MMReadFormat(filePath) )
        return false;

    // JPA: Shouldn't that just be an assertion check? It seems to me that
    // the user have to call clsparseHeaderfromFile before calling this function,
    // otherwise the whole pCsrMatrix will be broken;
    A->nrow = mm_reader.GetNumRows( );
    A->ncol = mm_reader.GetNumCols( );
    A->nnze = mm_reader.GetNumNonZeroes( );

    A->row = (int32_t *)malloc((A->nrow + 1) * sizeof(int32_t));
    A->val = (float *)malloc(A->nnze * sizeof(float));
    A->col = (int32_t *)malloc(A->nnze * sizeof(int32_t));

    if(A->row == NULL || A->col == NULL || A->val == NULL)
    {
        if(A->row == NULL)
            free((void *)A->row);
        if(A->col == NULL)
            free((void *)A->col);
        if(A->val == NULL)
            free((void *)A->val);
        return false;
    }

    //  The following section of code converts the sparse format from COO to CSR
    Coordinate* coords = mm_reader.GetUnsymCoordinates( );

    std::sort( coords, coords + A->nnze, CoordinateCompare );

    int32_t current_row = 1;

    A->row[ 0 ] = 0;

    for (int32_t i = 0; i < A->nnze; i++)
    {
        A->col[ i ] = coords[ i ].y;
        A->val[ i ] = coords[ i ].val;

        while( coords[ i ].x >= current_row )
            A->row[ current_row++ ] = i;
    }

    A->row[ current_row ] = A->nnze;

    while( current_row <= A->nrow )
        A->row[ current_row++ ] = A->nnze;

    return true;
}

void multiply_single(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
    C->nrow = A->nrow;
    C->ncol = B->ncol;
    C->val = (float *)malloc(C->nrow * C->ncol * sizeof(float));

    if(C->val == NULL)
        return;
    // TODO: Implement matrix multiplication with single thread. C=A*B
    for(int32_t i = 0; i < A->nrow; i++)
    {
        int32_t A_col_start = A->row[i];
        int32_t A_col_stop = A->row[i + 1];
        
        for(int32_t j = A_col_start; j < A_col_stop; j++)
        {
            int32_t B_row = A->col[j];

            for(int32_t k = 0; k < B->ncol; k++)
                C->val[i * C->ncol + k] += A->val[j] * B->val[B_row * B->ncol + k];
        }
    }
}
void multiply_pthread(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
    // TODO: Implement matrix multiplication with pthread. C=A*B
}

void multiply_openmp(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
    // TODO: Implement matrix multiplication with openmp. C=A*B
}

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char **argv)
{
    struct sparse_mtx A;
    if(!SCsrMatrixfromFile(&A, argv[1]))
    {
        std::cout << "read failed." << std::endl;
        return 0;
    }
    std::cout << "Matrix: " << argv[1] << std::endl;

    struct dense_mtx B;
    B.nrow = A.ncol;
    B.ncol = atoi(argv[2]);
    if(B.ncol < 0)
    {
        free(A.row);
        free(A.col);
        free(A.val);
        std::cerr << "Invalid argument for the number of columns of B." << std::endl;
    }
    B.val = (float *)malloc(sizeof(float) * B.nrow * B.ncol);

    srand((unsigned int)time(NULL));
    for(int i = 0; i < B.nrow; i++)
    {
        for(int j = 0; j < B.ncol; j++)
        {
            B.val[B.ncol * i + j] = ((float)rand()/(float)(RAND_MAX)) * ((rand() % 2) ? 1.0f : -1.0f);
        }
    }

    struct dense_mtx C1, C2;
    C1.val = NULL;
    C2.val = NULL;

    std::cout << "Single Thread Computation Start" << std::endl;
    uint64_t start = GetTimeStamp();
    multiply_single(&A, &B, &C1);
    uint64_t end = GetTimeStamp();
    std::cout << "Single Thread Computation End: " << end - start  << " us." << std::endl;
    /*
    std::cout << "Pthread Computation Start" << std::endl;
    start = GetTimeStamp();
    multiply_pthread(&A, &B, &C2);
    end = GetTimeStamp();
    std::cout << "Pthread Computation End: " << end - start << " us." << std::endl << std::endl;
    */
    std::cout << "OpenMP Computation Start" << std::endl;
    start = GetTimeStamp();
    multiply_openmp(&A, &B, &C2);
    end = GetTimeStamp();
    std::cout << "OpenMP Computation End: " << end - start << " us." << std::endl << std::endl;

    // TODO: Testing Code by comparing C1 and C2

    free(A.row);
    free(A.col);
    free(A.val);
    free(B.val);
    if(C1.val != NULL)
        free(C1.val);
    if(C2.val != NULL)
        free(C2.val);
    
    return 0;
}
