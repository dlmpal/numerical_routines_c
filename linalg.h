//
// Created by User on 18/8/2020.
//

#ifndef LINALG_LINALG_H
#define LINALG_LINALG_H
#include <math.h>
#include <stdio.h>

double determinant(int rows, int cols, double matrix[rows][cols]);

void transpose(int rows , int cols , double matrix[][cols], double transpose_matrix[][rows]);

void print_matrix(int rows , int cols , double matrix[][cols]);

void print_vector(int n , double vector[n]);

void matrix_add(int rows , int cols , double matrix1[][cols],double matrix2[][cols],double sum_matrix[][cols]);

void matrix_mult(int rows1 , int cols1 , int cols2 , double matrix1[][cols1] , double matrix2[][cols2],double mult_matrix[rows1][cols2]);

void matrix_times_vector(int rows , int cols , double  matrix[][cols] , double vector[cols] , double product[rows]);

void vector_times_matrix(int rows , int cols , double matrix[rows][cols] , double vector[rows],double product[cols] );

void create_identity_matrix(int n , double matrix[][n]);

void copy_matrix(int rows , int cols , double original[][cols] , double copy[][cols]);

void copy_vector(int n , double original[n] , double copy[n]);

void solve_system(int n , double A_matrix[][n], double B_vector[n] , double AM[][n] , double BM[n]);

void raise_matrix(int n, int power , double matrix[][n],double raised_matrix[][n]);

void invert_matrix(int n , double matrix[][n] ,double inverted_matrix[][n] ,double identity[][n]);

void matrix_times_scalar(int rows , int cols , double matrix[][cols],double scalar);

void vector_times_scalar(int dim , double vector[dim] , double scalar);

void vector_addition(int dim , double vector1[dim],double vector2[dim]);

void vector_zeroes(int n , double vector[n]);

void nabs(double n );

void create_vandermonde_matrix(int rows , int cols , double val_vector[rows] , double v_matrix[][cols]);

double matrix_element_sum(long rows , long cols , double matrix[][cols]  );

double vector_element_sum(long n ,  double vector[n]);

void split_matrix_to_four(int dim , double matrix[][dim] , double a[][dim/2] , double b[][dim/2] , double c[][dim/2],double d[][dim/2]);

void strassen_mult(int dim , double matrix1[][dim] , double matrix2[][dim] , double result[][dim]);



#endif //LINALG_LINALG_H
