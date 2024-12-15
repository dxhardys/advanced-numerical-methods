# pragma once
#include <stdio.h>
#include <stdlib.h>

void dgemm_ijk(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for(int i =0 ; i <n ; i++)
    {
        for(int j=0; j<p ; j++)
        {
            double sum = 0.0 ;
            for(int k =0; k<m;k++)
            {
                sum += A[i][k] * B[k][j] ;
            }
            C[i][j] = sum ;
        }
    }          
}
void dgemm_ikj(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for(int i =0 ; i <n ; i++)
    {
        for(int k =0; k<m;k++)
        {
            double const _A_ = A[i][k] ;
            for(int j=0; j<p ; j++)
            {
                C[i][j] += _A_ * B[k][j] ;
            }
        }
    }          
}

void dgemm_jik(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for(int j=0; j<p ; j++)
    {
        for(int i=0; i<n ; i++)
        {
            double sum = 0.0 ;
            for(int k =0; k<m;k++)
            {
                sum += A[i][k] * B[k][j] ;
            }
            C[i][j] = sum ;
        }
    }       
}

void dgemm_jki(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for(int j=0; j<p ; j++)
    {
        for(int k =0; k<m; k++)
        {
            double const _B_ = B[k][j] ;
            for(int i =0; i<n; i++)
            {
                C[i][j] += A[i][k] * _B_ ;
            }
        }
    }       
}


void dgemm_kij(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for(int k =0; k<m; k++)
    {
        for(int i =0; i<n; i++)
        {
            double const _A_ = A[i][k] ;
            for(int j =0; j<p; j++)
            {
                C[i][j] +=  _A_ * B[k][j] ;
            }
        }
    }       
}

void dgemm_kji(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for(int k =0; k<m; k++)
    {
        for(int j =0; j<p; j++)
        {
            double const _B_ = B[k][j] ;
            for(int i =0; i<n; i++)
            {
                C[i][j] += A[i][k] * _B_ ;
            }
        }
    }       
}




void print_matrix(int n, int m, double A[n][m]){
    for(int i =0 ; i <n ; i++)
    {
        for(int j=0; j<m;j++)
        {
            printf("%f \t",A[i][j]);
        }
    printf("\n");
    }
}


void fill_matrix(int n, int m, double A[n][m])
{
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < m; k++)
        {
            A[i][k] = 1.0;
        }
    }
}


