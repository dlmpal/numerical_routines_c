//
// Created by User on 21/8/2020.
//

#ifndef LINALG_C_ALG_H
#define LINALG_C_ALG_H

#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//============================================================//

void print_matrix_complex(int rows, int cols, double complex matrix[][cols]) ;

void print_vector_complex(int n, double complex vector[n]) ;

double complex determinant_complex(int rows, int cols, double complex matrix[rows][cols]);

void transpose_complex(int rows, int cols, double complex matrix[][cols], double complex transpose_matrix[][rows]) ;

void matrix_mult_complex(int rows1, int cols1, int cols2, double complex matrix1[][cols1], double complex matrix2[][cols2],double complex mult_matrix[rows1][cols2]) ;

void matrix_add_complex(int rows, int cols, double complex matrix1[][cols], double complex matrix2[][cols],double complex sum_matrix[][cols]) ;

void matrix_times_vector_complex(int rows, int cols, double complex matrix[][cols], double complex vector[cols],double complex product[rows]) ;

void vector_times_matrix_complex(int rows, int cols, double complex matrix[rows][cols], double complex vector[rows],double complex product[cols]) ;

void create_identity_matrix_complex(int n, double complex matrix[][n]) ;

void copy_matrix_complex(int rows, int cols, double complex original[][cols], double complex copy[][cols]) ;

void copy_vector_complex(int n, double complex original[n], double complex copy[n]);

void raise_matrix_complex(int n, int power, double complex matrix[][n], double complex raised_matrix[][n]) ;

void matrix_times_scalar_complex(int rows, int cols, double complex matrix[][cols], double complex scalar);

void vector_times_scalar_complex(int dim, double complex vector[dim], double complex scalar) ;

void vector_addition_complex(int dim, double complex vector1[dim], double complex vector2[dim]) ;

void vector_zeroes_complex(int n, double complex vector[n]) ;

void solve_system_complex(int n, double complex A_matrix[][n], double complex B_vector[n], double complex AM[][n],complex double BM[n]);

double complex matrix_element_sum_complex(long rows, long cols, double complex matrix[][cols]) ;

double complex vector_element_sum_complex(long n, double complex vector[n]) ;

void create_dft_matrix(int dim, double complex dft_matrix[][dim]) ;

void create_idft_matrix(int dim, double complex dft_matrix[][dim]) ;

void dft_slow(int dim, double complex f_vector[dim], double complex f_hat_vector[dim]) ;

void idft_slow(int dim, double complex f_hat_vector[dim], double complex f_vector[dim]);

void fft(int dim, double complex f_vector[dim], double complex f_hat_vector[dim]);

#endif //LINALG_C_ALG_H
