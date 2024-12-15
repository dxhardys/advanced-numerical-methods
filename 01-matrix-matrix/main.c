
#include <time.h>
#include <stdlib.h>
#include "include/matrix.h"
#include "include/matrix_flat.h"
#include "include/flat_restrict.h"    


void fill_random(double* A, int n){

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n ; j++)
        {
            A[i * n + j] = rand() ;
        }   
    }
}


int main(int argc, char *argv[])
{
    int n = 100;

    // Parse dimensions of matrices
    if (argc > 1)
    {
        n = atoi(argv[1]);
    }

    int m = n;
    int p = n;

    if (argc > 2)
    {
        m = atoi(argv[2]);
    }

    if (argc > 3)
    {
        p = atoi(argv[3]);
    }

    printf("Working with matrices of size %d x %d and %d x %d\n", n, m, m, p);

    // Allocate memory for matrices
    // Note how we allocate memory for a 2D array in a single block
    // This is done to ensure that the memory is contiguous
    // Can you comment on how data are stored in memory for this 2D array?
    double **A = (double **)malloc(n * m * sizeof(double));
    double **B = (double **)malloc(m * p * sizeof(double));
    double **C = (double **)malloc(n * p * sizeof(double));

    // Initialize matrices
    fill_matrix(n, m, A);
    fill_matrix(m, p, B);
    fill_matrix(n, p, C);


    // Call dgemm_ijk
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    dgemm_ijk(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken dgemm ijk: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    dgemm_ikj(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken dgemm ikj: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    dgemm_jik(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken dgemm jik: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    dgemm_jki(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken dgemm jki: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    dgemm_kij(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken dgemm kij: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    dgemm_kji(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken dgemm kji: %e s\n", time_taken);

    // using flat version 
    double* a = (double*)malloc(n * m * sizeof(double));
    double* b = (double*)malloc(m * p * sizeof(double));
    double* c = (double*)malloc(n * p * sizeof(double));

    fill_matrix_flat(n,m,a);
    fill_matrix_flat(m,p,b);

    


    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flat_dgemm_ijk(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat dgemm kji: %e s\n", time_taken);


    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flat_dgemm_ikj(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat dgemm ikj: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flat_dgemm_jki(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat dgemm jki: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flat_dgemm_jik(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat dgemm jik: %e s\n", time_taken);


    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flat_dgemm_kij(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat dgemm kij: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flat_dgemm_kji(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat dgemm kji: %e s\n", time_taken);

    // Flatten restricted dgemm
    size_t size_aa = (sizeof(double) * n * m) ;
    double* restrict aa = (double*)aligned_alloc(64,size_aa) ;

    size_t size_bb = (sizeof(double) * m * p) ;
    double* restrict bb = (double*)aligned_alloc(64,size_bb) ;

    size_t size_cc = (sizeof(double) * n * p) ;
    double* restrict cc = (double*)aligned_alloc(64,size_cc) ;

    fill_matrix_flatr(n,m,aa);
    fill_matrix_flatr(m,p,bb);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flatr_dgemm_ijk(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat restrict dgemm ijk: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flatr_dgemm_ikj(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat restrict dgemm ikj: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flatr_dgemm_jki(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat restrict dgemm jki: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flatr_dgemm_jik(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat restrict dgemm jik: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flatr_dgemm_kij(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat restrict dgemm kij: %e s\n", time_taken);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    flatr_dgemm_kji(n,m,p,a,b,c);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken flat restrict dgemm kji: %e s\n", time_taken);


    // Exercice 3 
    /*We want to benchmark the function matrix_matrix_multiplication.
    For the case of square matrices m = n = p, generate two random
    matrices A and B and compute their product. Use the time command to
    measure the time it takes to compute the product.*/

    FILE* file = fopen("../results.txt","w");
    if(!file){
        fprintf(stderr,"File not found, make sure it exists");
        return EXIT_FAILURE ;
    } 

    srand(time(NULL)) ;
    int N = 1000 ;
    for (int i = 2 ; i < N; i++)
    {
        printf("Iteration %d x %d : \n",i,i) ;
        double* a_rand = (double*)malloc(i * i * sizeof(double));
        double* b_rand = (double*)malloc(i * i * sizeof(double));
        double* c_rand = (double*)malloc(i * i * sizeof(double));

        fill_random(a_rand,i);
        fill_random(b_rand,i);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        flat_dgemm_ijk( i, i, i,  a_rand, b_rand, c_rand);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
        fprintf(file,"%d %lf \n",i,time_taken);

        free(a_rand);
        free(b_rand);
        free(c_rand);

    }


    
    fclose(file) ;
    // Free memory
    free(A);
    free(B);
    free(C);
    free(a);
    free(b);
    free(c);
    free(aa);
    free(bb);
    free(cc);

    return EXIT_SUCCESS;
}