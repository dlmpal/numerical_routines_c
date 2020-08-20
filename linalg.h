//
// Created by User on 18/8/2020.
//

#ifndef LINALG_LINALG_H
#define LINALG_LINALG_H
#include <math.h>
#include <stdio.h>
#include <stdio.h>
double determinant(int rows, int cols, double matrix[rows][cols]) {
    double det = 0;

    if (rows != cols || rows == 1 || rows == 0) { // logic check - determinant is defined only for square matrices
        printf("Matrix must be square");
        return 0;
    }
    if (rows == 2) { //base case - ends recursion
        double result = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
        return result;
    } else {
        register int i ;
        register int j ;
        register int k ;
        register int row_count;
        register int col_count;
        for ( i = 0; i < cols; i++) {
             double minor[rows - 1][cols - 1];
            row_count = 0;
            for ( j = 1; j < rows; j++) {
                 col_count = 0;
                for ( k = 0; k < rows; k++) {
                    if (k == i) { continue; }
                    minor[row_count][col_count] = matrix[j][k];
                    col_count++;
                }
                row_count++;
            }

            //printf("i is %d and det is %lf\n",i,det); iteration check
            det += pow(-1, i) * matrix[0][i] * determinant(rows - 1, cols - 1, minor);
        }
        return det;
    }
}

void transpose(int rows , int cols , double matrix[][cols], double transpose_matrix[][rows]){
    register int i ;
    register int j;
    for( i = 0 ; i < rows ; i++){
        for( j = 0 ; j < cols ; j++){
            transpose_matrix[j][i] = matrix[i][j];
        }
    }
}

void print_matrix(int rows , int cols , double matrix[][cols]){
    printf("\n");
    register int i ;
    register int j;
    for( i = 0 ; i < rows ; ++ i ){
        printf("\n");
        for( j = 0 ; j < cols ; ++ j){
            printf("%lf ",matrix[i][j]);
        }
    }
}

void print_vector(int n , double vector[n]){
    printf("\n");
    for(int i = 0 ; i < n ; i++){
        printf("%lf ",vector[i]);
    }
}
void matrix_add(int rows , int cols , double matrix1[][cols],double matrix2[][cols],double sum_matrix[][cols]){
    register int i ;
    register int j ;
    for( i = 0 ; i < rows ; i++){
        for( j = 0 ; j < rows ; j++){
            sum_matrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
}

void matrix_mult(int rows1 , int cols1 , int cols2 , double matrix1[][cols1] , double matrix2[][cols2],double mult_matrix[rows1][cols2]){
    if(rows1<=0 || cols2 <= 0 || cols1 <= 0){
        printf("Wrong format");
    }
    else{
        register int i ;
        register int j;
        register int k;
        double local_sum = 0 ;
        for( i = 0 ; i < rows1 ; i ++){
            for( j = 0 ; j < cols2 ; j++) {
                for ( k = 0; k < cols1; k++) {
                    local_sum += matrix1[i][k] * matrix2[k][j];
                }
                mult_matrix[i][j] = local_sum;
                local_sum = 0 ;
            }
        }
    }
}

void create_identity_matrix(int n , double matrix[][n]){
    register int i ;
    register int j ;
    for( i = 0; i < n ; i++){
        for( j = 0 ; j < n ; j++){
            matrix[i][j] = (i==j) ? 1 : 0 ;
        }
    }
}

void copy_matrix(int rows , int cols , double original[][cols] , double copy[][cols]){
    register int i ; register int j;
    for( i = 0 ; i < rows ; i++){
        for( j = 0 ; j < cols ; j++){
            copy[i][j] = original[i][j];
        }
    }
}

void copy_vector(int n , double original[n] , double copy[n]){
    register int i;
    for( i = 0 ; i < n ; i++){
        copy[i] = original[i];
    }
}

void solve_system(int n , double A_matrix[][n], double B_vector[n] , double AM[][n] , double BM[n]){
    copy_matrix(n,n,A_matrix,AM);
    copy_vector(n,B_vector,BM);
    double scaler ;
    register int fd;
    register int i ;
    register int k ;
    for( fd = 0 ; fd < n ; fd++) {
        scaler = 1 / (AM[fd][fd]);
        for (int j = 0; j < n; j++) {
            AM[fd][j] *= scaler;
        }
        BM[fd] *= scaler;
        double current_scaler;
        for ( i = 0; i < fd; i++) {
            current_scaler = AM[i][fd];
            for ( k = 0; k < n; k++) {
                AM[i][k] = AM[i][k] - current_scaler * AM[fd][k];
            }
            BM[i] = BM[i] - current_scaler * BM[fd];
        }
        for ( i = fd + 1; i < n; i++) {
            current_scaler = AM[i][fd];
            for ( k = 0; k < n; k++) {
                AM[i][k] -= current_scaler * AM[fd][k];
            }
            BM[i] -= current_scaler * BM[fd];
        }
    } // Now the matrix AM is morphed into a nxn identity matrix
    // while BM contains a vector form solution of the system
}

void raise_matrix(int n, int power , double matrix[][n],double raised_matrix[][n]){
    if(power<0){
        printf("Power must be positive or zero ");
    }
    if(power==0){
        create_identity_matrix(n,raised_matrix);
    }
    else{
        create_identity_matrix(n,raised_matrix);
        register int i ;
        for( i = 0; i < power  ; i++){
            double current[n][n];
            copy_matrix(n,n,raised_matrix,current);
            matrix_mult(n,n,n,current,matrix,raised_matrix);
        }
    }
}

void invert_matrix(int n , double matrix[][n] ,double inverted_matrix[][n] ,double identity[][n]){
    if(determinant(n,n,matrix) == 0.0 ){
        printf("Matrix is singular-cannot be inverted");
    }
    //double inverted_matrix[n][n];
    copy_matrix(n,n,matrix,inverted_matrix);
    create_identity_matrix(n,identity);
    print_matrix(n,n,inverted_matrix);
    register int j ;
    register int i ;
    register int k ;
    for(int fd = 0 ; fd < n ; fd++){
        double diag_scaler = 1.0 / matrix[fd][fd];
        for( j = 0 ; j < n  ; j++){
            inverted_matrix[fd][j] *= diag_scaler;
            identity[fd][j] *= diag_scaler;
        }
        double current_scaler;

        for( i = 0 ; i < fd ; i++){
            current_scaler = inverted_matrix[i][fd];
            for( k = 0 ; k < n ; k++){
                inverted_matrix[i][k] -= current_scaler*inverted_matrix[fd][k];
                identity[i][k] -= current_scaler*identity[fd][k];
            }
        }
        for( i = fd+1 ; i < n ; i++){
            current_scaler = inverted_matrix[i][fd];
            for( k = 0; k < n ; k++){
                inverted_matrix[i][k] -= current_scaler*inverted_matrix[fd][k];
                identity[i][k] -= current_scaler*identity[fd][k];
            }
        }

    }

}

void matrix_times_scalar(int rows , int cols , double matrix[][cols],double scalar){
    register int i ;
    register int j;
    for( i = 0 ; i < rows ; ++i){
        for( j = 0 ; j < cols ; j++){
            matrix[i][j] *= scalar;
        }
    }
}
void vector_times_scalar(int dim , double vector[dim] , double scalar){
    register int i ;
    for( i = 0 ; i < dim ; i++){
        vector[i] *= scalar;
    }
}
void vector_addition(int dim , double vector1[dim],double vector2[dim]){
    register int i ;
    for( i = 0 ; i < dim ; i++){
        vector1[i] += vector2[i];
    }
}
void vector_zeroes(int n , double vector[n]){
    for(int i = 0 ; i < n ; i++){
        vector[i] = 0 ;
    }
}
void nabs(double n ){
    n = (n>=0) ? n : -n;
}
long get_num_of_instances(double t_zero , double t_final , double time_step){
    long num_of_instances = (long) (t_final-t_zero)/time_step;
    return num_of_instances;
}

#endif //LINALG_LINALG_H
