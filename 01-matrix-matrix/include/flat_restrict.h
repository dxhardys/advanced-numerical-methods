#pragma once
#include <stdio.h>
// Flatten Matrix Matrix multiply

void fill_matrix_flatr(int n, int m, double *restrict A)
{
    for(int i=0; i<n ;i++)
    {
        for (int j = 0; j < m; j++)
        {
            A[i * m + j] = 1.00 ;
        }
        
    }
}

void print_matrix_flatr(int n, int m, double *restrict A){
    for(int i =0 ; i <n ; i++)
    {
        for(int j=0; j<m;j++)
        {
            printf("%f \t",A[i * m + j]);
        }
    printf("\n");
    }
}

void flatr_dgemm_ijk(int n, int m, int p, double* A, double* B, double* C)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            double sum = 0.0 ;
            for (int k = 0; k < p; k++)
            {
                sum+= A[i * m + k] * B[k * p + j] ;
            }  
            C[i * p + j] = sum ;
        }   
    }
}

// loop interchange + invariant extraction
void flatr_dgemm_ikj(int n, int m, int p, double *restrict A, double *restrict B, double *restrict C)
{
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < p; k++)
        {
            const double _A_ = A[i * m + k];
            for (int j = 0; j < m; j++)
            {
                C[i * p + j] += _A_  * B[k * p + j] ;
            }   
        }   
    }
}

void flatr_dgemm_jki(int n, int m, int p, double *restrict A, double *restrict B, double *restrict C)
{
    for  (int j = 0; j < m; j++)
    {
        for (int k = 0; k < p; k++)
        {
            double const _B_ = B[k * p + j] ;
            for (int i = 0; i < n; i++)
            {
                C[i * p + j] += A[i * m + k] * _B_;
            }   
        }   
    }
}

void flatr_dgemm_jik(int n, int m, int p, double *restrict A, double *restrict B, double *restrict C)
{
    for  (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0 ;
            for (int k = 0; k < p; k++)
            {
                sum += A[i * m + k] * B[k * p + j];
            }  
            C[i * p + j] = sum ; 
        }   
    }
}

void flatr_dgemm_kij(int n, int m, int p, double *restrict A, double *restrict B, double *restrict C)
{
    for (int k = 0; k < p; k++)
    {
        for (int i = 0; i < n; i++)
        {
            double const _A_ = A[i * m + k]  ;
            for  (int j = 0; j < m; j++)
            {
                C[i * p + j] += _A_ * B[k * p + j];
            }  
        }   
    }
}

void flatr_dgemm_kji(int n, int m, int p, double *restrict A, double *restrict B, double *restrict C)
{
    for (int k = 0; k < p; k++)
    {
        for  (int j = 0; j < m; j++)
        {
            double const _B_ = B[k * p + j] ;
            for (int i = 0; i < n; i++)
            {
                C[i * p + j] += A[i * m + k]  * _B_;
            }  
        }   
    }
}