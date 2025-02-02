#include <time.h>
#include <stdbool.h>
#include "utils.h"
#include <math.h>

// convert flatten sparse matrix into a CSR format
void flat2CSR(double* A,int row, int col, double** values, int** col_indices, int** row_ptr, int* nnz){

    (*nnz) = 0 ;
    for (int i = 0; i < row; i++){
        for(int j = 0 ; j < col ; j++){
            if(A[i * col + j] != 0){
                (*nnz)++ ;
            }
        }
    }

    (*values) = (double*)malloc((*nnz) * sizeof(double));
    (*col_indices) = (int*)malloc((*nnz) * sizeof(int)) ;
    (*row_ptr) = (int*)malloc( ((*nnz) + 1) * sizeof(int)) ;

    (*row_ptr)[0] = 0 ;
    int index = 0 ;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            double current = A[i * col + j] ;
            if(current != 0){
                (*values)[index] = current ;
                (*col_indices)[index] = j  ;
                index++ ;
            }
        }
        (*row_ptr)[i+1] = index ;
    }
    
  

}

void solve_jacobi_CSR(double* values, int* col_indices, int* row_ptr, int nnz, int size, double* RHS, double* X){
    printf("\e[31mJacobi\e[0m\n") ;

    double tmp[size] ;
    int n_iter = 30 ;
    int k = 0 ;
    double epsilon = 0.000001 ;

    while(k < n_iter){
    printf("\e[32mIteration %d\e[0m\n",k) ;
    for (int i = 0; i < size ; i++)
    {
        int start = row_ptr[i] ;
        int end  = row_ptr[i+1] ;
        double sum = 0.0 ;
        double aii = 0.0 ;
        for (int j = start; j < end; j++)
        {
            if(i == col_indices[j]){
                aii = values[j] ;
                continue ;
            }
            sum += values[j] * X[col_indices[j]] ;
            
        }
        tmp[i] = ( RHS[i] - sum ) / aii ;
        printf("X[%i] = %lf\n",i,tmp[i]) ;

    }

    double max = 0.0 ;
    for (int i = 0; i < size; i++)
    {
        double diff = fabs(tmp[i] - X[i]) ;
        if(diff > max)
            max = diff ;
        
    X[i] = tmp[i] ;
    }
    if(max < epsilon)
        break ;
    k++;
    }

}

void solve_gauss_seidel_CSR(double* values, int* col_indices, int* row_ptr, int nnz, int size, double* RHS, double* X){

    printf("\e[31mGauss-Seidel\e[0m\n") ;
    double tmp[size] ;
    int n_iter = 30 ;
    int k = 0 ;
    double epsilon = 0.000001 ;

    while(k < n_iter){
    printf("\e[32mIteration %d\e[0m\n",k) ;

    for (int i = 0; i < size; i++) {
            tmp[i] = X[i];
        }
        
    for (int i = 0; i < size ; i++)
    {
        int start = row_ptr[i] ;
        int end  = row_ptr[i+1] ;
        double sum = 0.0 ;
        double aii = 0.0 ;
        for (int j = start; j < end; j++)
        {
            if(i == col_indices[j]){
                aii = values[j] ;
                continue ;
            }
            sum += values[j] * X[col_indices[j]] ;
            
        }
        X[i] = ( RHS[i] - sum ) / aii ;
        printf("X[%i] = %lf\n",i,X[i]) ;

    }

    double max = 0.0 ;
    for (int i = 0; i < size; i++)
    {
        double diff = fabs(X[i] - tmp[i])  ;
        if(diff > max)
            max = diff ;
    }

    if(max < epsilon)
        break ;
    
    k++;
    
    }
}





int main(int argc, char ** argv){

    srand(time(NULL)) ;
    // Generate a n * m sparse matrix
    int n = 5 ;
    int m = 5 ;
    double* A = (double*)malloc(n * m * sizeof(double));
    generate_sparse_matrix(A,n,m) ;
    //print_sparse(A,n,m) ;

    double D[9] = {1.0, 0.0, 2.0, 
                   0.0, -1.0, 4.0, 
                   3.0, 0.0, 2.0} ;

    double* values = NULL ;
    int* row_pointers = NULL ;
    int* col_indices = NULL ;
    int nnz ;

    flat2CSR(A,n,m,&values,&col_indices,&row_pointers,&nnz) ;

    double B[9] = {5, -2, 3,
                  -3,  9, 1,
                   2, -1, -7} ;

    double RHS[3] = {-1, 2, 3} ;
    double X[3]   = {0.0, 0.0, 0.0} ; 
    int size = 3 ;

    print_sparse(D,size,size) ;

    values = NULL ;
    row_pointers = NULL ;
    col_indices = NULL ;
    nnz = 0 ;

    flat2CSR(B,size,size,&values,&col_indices,&row_pointers,&nnz) ;
    print_CSR(values,row_pointers,col_indices,nnz,size,size) ;
    // Makes sure diagonal elements arent 0 
    solve_jacobi_CSR(values,col_indices,row_pointers,nnz,size,RHS,X) ;



    values = NULL ;
    row_pointers = NULL ;
    col_indices = NULL ;
    nnz = 0 ;
    double B2[9] = {5, -2, 3,
                  -3,  9, 1,
                   2, -1, -7} ;

    double RHS2[3] = {-1, 2, 3} ;
    double X2[3]   = {0.0, 0.0, 0.0} ; 
    flat2CSR(B2,size,size,&values,&col_indices,&row_pointers,&nnz) ;
    solve_gauss_seidel_CSR(values,col_indices,row_pointers,nnz,size,RHS2,X2) ;






    

    return 0 ;
}