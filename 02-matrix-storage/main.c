#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>


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

// convert flatten sparse matrix into a CSR format
// Compressed Sparse Row
// takes into input 
void CSR(double* A, int n, int m, double* value, int* row_ptr, int* col_ind){

    int size = 0;
    int nnz = 0;
    int col_ind_size = 0;
    int row_ptr_size = n + 1; 


    col_ind = (int*)malloc(col_ind_size * sizeof(int));
    value = (double*)malloc(size * sizeof(double));
    row_ptr = (int*)malloc(row_ptr_size * sizeof(int));
    row_ptr[0] = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (A[i * m + j] != 0) {
                size++;
                value = (double*)realloc(value, size * sizeof(double));
                value[size - 1] = A[i * m + j];

                col_ind_size++;
                col_ind = (int*)realloc(col_ind, col_ind_size * sizeof(int));
                col_ind[col_ind_size - 1] = j;

                nnz++;
            }
        }

        row_ptr[i + 1] = nnz;
}


    printf("Value tab = \n") ;
    for (int i = 0; i < nnz; i++)
    {
        printf("%lf\t ",value[i]) ;
    }
    printf("\n") ;

    printf("col_ind_size tab = \n") ;
    for (int i = 0; i < nnz; i++)
    {
        printf("%d\t ",col_ind[i]) ;
    }

    printf("\n") ;
    printf("row_ptr tab = \n") ;
    for (int i = 0; i < nnz; i++)
    {
        printf("%d\t ",row_ptr[i]) ;
    }
    
    
}

void print_sparse(double* A,int n, int m){
    for(int i = 0 ; i<n ; i++)
    {
        for(int j = 0 ; j<m ; j++){
            printf("%f\t",A[i * m + j]) ;
        }
        printf("\n") ;
    }
}



int main(int argc, char ** argv){
    srand(time(NULL)) ;


    // Generate a n * m sparse matrix
    int n = 5 ;
    int m = 4 ;
    double* A = (double*)malloc(n * m * sizeof(double));
    generate_sparse_matrix(A,n,m) ;
    print_sparse(A,n,m) ;

    double B[] = {1, 2, 0, 3, 0, 4, 5, 0, 7};


    int *values, *row_pointers, *col_indices ;

    CSR(A,n,m,values,row_pointers,col_indices);






    

    return 0 ;
}