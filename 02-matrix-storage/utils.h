#pragma once 
#include <stdio.h>
#include <stdlib.h>

void print_sparse(double* A,int n, int m){
    for(int i = 0 ; i<n ; i++)
    {
        for(int j = 0 ; j<m ; j++){
            printf("%lf\t",A[i * m + j]) ;
        }
        printf("\n") ;
    }
}


void print_CSR(double* values, int* row_ptr, int* col_indices, int nnz, int n, int m){

    printf("Values =[") ;
    for (int i = 0; i < nnz; i++)
    {
        printf("%lf ",values[i]) ;
    }
    printf("]\n") ;

    printf("col_indices = [") ;
    for (int i = 0; i < nnz; i++)
    {
        printf("%d ",col_indices[i]) ;
    }
    printf("]\n") ;


    printf("row_ptr = [") ;
    for (int i = 0; i < n + 1; i++)
    {
        printf("%d ",row_ptr[i]) ;
    }
    printf("]\n") ;
}


// takes into input a flatten matrix 
void generate_sparse_matrix(double* A, int n, int m){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            double r = rand() % 100;
            int proba = rand() % 10 ;
            A[i * m + j] = (proba < 7) ? 0.0 : r ;
        }
        
    }
    
}

