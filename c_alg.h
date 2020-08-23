//
// Created by User on 21/8/2020.
//

#ifndef LINALG_C_ALG_H
#define LINALG_C_ALG_H

#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef struct {
    double real, img;
} cnum;

void add_cnum(cnum c1, cnum c2, cnum c3) {
    c3.real = c1.real + c2.real;
    c3.img = c1.img + c2.img;
}

void mult_cnum(cnum c1, cnum c2, cnum c3) {
    //mult_vector = vector([c1.real * c2.real - c1.im * c2.im, c1.im * c2.real + c1.real * c2.im])
    c3.real = c1.real * c2.real - c1.img * c2.img;
    c3.img = c1.img * c2.real + c1.real * c2.img;
}

void c_exp(double theta, cnum c1) {
    //theta in radians
    c1.real = cos(theta);
    c1.img = sin(theta);
}

void cnum_polar(cnum c1, double mag, double theta) {
    mag = sqrt(pow(c1.real, 2) + pow(c1.img, 2));
    theta = atan(c1.img / c1.real);
}

void cnum_cart(cnum c1, double mag, double theta) {
    c1.real = mag * cos(theta);
    c1.img = mag * sin(theta);
}

void div_cnum(cnum c1, cnum c2, cnum c3) {
    double theta1, theta2, mag1, mag2;
    cnum_polar(c1, mag1, theta1);
    cnum_polar(c2, mag2, theta2);
    cnum_cart(c3, mag2 / mag1, theta2 - theta1);
}
//============================================================//

void print_matrix_complex(int rows, int cols, double complex matrix[][cols]) {
    printf("\n");
    register int i;
    register int j;
    for (i = 0; i < rows; ++i) {
        printf("\n");
        for (j = 0; j < cols; ++j) {
            printf(" %lf + i%lf ", creal(matrix[i][j]), cimag(matrix[i][j]));
        }
    }
}

void print_vector_complex(int n, double complex vector[n]) {
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf(" %lf +i%lf ", creal(vector[i]), cimag(vector[i]));
    }
}

double complex determinant_complex(int rows, int cols, double complex matrix[rows][cols]) {
    double complex det = 0 + 0 * I;
    if (rows != cols || rows == 1 || rows == 0) { // logic check - determinant is defined only for square matrices
        printf("Matrix must be square");
        return 0;
    }
    if (rows == 2) { //base case - ends recursion
        double complex result = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
        return result;
    } else {
        register int i;
        register int j;
        register int k;
        register int row_count;
        register int col_count;
        for (i = 0; i < cols; i++) {
            double complex minor[rows - 1][cols - 1];
            row_count = 0;
            for (j = 1; j < rows; j++) {
                col_count = 0;
                for (k = 0; k < rows; k++) {
                    if (k == i) { continue; }
                    minor[row_count][col_count] = matrix[j][k];
                    col_count++;
                }
                row_count++;
            }

            //printf("i is %d and det is %lf\n",i,det); iteration check
            det += pow(-1, i) * matrix[0][i] * determinant_complex(rows - 1, cols - 1, minor);
        }
        return det;
    }
}

void transpose_complex(int rows, int cols, double complex matrix[][cols], double complex transpose_matrix[][rows]) {
    register int i;
    register int j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            transpose_matrix[j][i] = matrix[i][j];
        }
    }
}


void matrix_mult_complex(int rows1, int cols1, int cols2, double complex matrix1[][cols1], double complex matrix2[][cols2],
                    double complex mult_matrix[rows1][cols2]) {
    if (rows1 <= 0 || cols2 <= 0 || cols1 <= 0) {
        printf("Wrong format");
        return;
    } else {
        register int i;
        register int j;
        register int k;
        for (i = 0; i < rows1; i++) {
            for (j = 0; j < cols2; j++) {
                mult_matrix[i][j] = 0. + 0.*I;
                for (k = 0; k < cols1; k++) {
                    mult_matrix[i][j] += matrix1[i][k] * matrix2[k][j];
                }

            }
        }
    }
}


void matrix_mult_cnum(int rows1, int cols1, int cols2, cnum **matrix1, cnum **matrix2,
                      cnum **mult_matrix) {
    if (rows1 <= 0 || cols2 <= 0 || cols1 <= 0) {
        printf("Wrong format");
        return;
    } else {
        register int i;
        register int j;
        register int k;

        for (i = 0; i < rows1; i++) {
            for (j = 0; j < cols2; j++) {
                for (k = 0; k < cols1; k++) {
                    mult_cnum(matrix1[i][k],matrix2[k][j],mult_matrix[i][j]);
                }
            }
        }
    }
}

void matrix_add_complex(int rows, int cols, double complex matrix1[][cols], double complex matrix2[][cols],
                        double complex sum_matrix[][cols]) {
    register int i;
    register int j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < rows; j++) {
            sum_matrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
}

void matrix_times_vector_complex(int rows, int cols, double complex matrix[][cols], double complex vector[cols],
                                 double complex product[rows]) {
    register int i;
    register int j;
    for (i = 0; i < rows; i++) {
        product[i] = 0;
        for (j = 0; j < cols; j++) {
            product[i] += matrix[i][j] * vector[j];
        }
    }
}

void vector_times_matrix_complex(int rows, int cols, double complex matrix[rows][cols], double complex vector[rows],
                                 double complex product[cols]) {
    register int i;
    register int j;
    register double complex current_sum;
    for (j = 0; j < cols; j++) {
        current_sum = 0. + 0. * I;
        for (i = 0; i < rows; i++) {
            current_sum += vector[i] * matrix[i][j];
        }
        product[j] = current_sum;
    }
}

void create_identity_matrix_complex(int n, double complex matrix[][n]) {
    register int i;
    register int j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matrix[i][j] = (i == j) ? 1 + 0 * I : 0 + 0 * I;
        }
    }
}

void copy_matrix_complex(int rows, int cols, double complex original[][cols], double complex copy[][cols]) {
    register int i;
    register int j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            copy[i][j] = original[i][j];
        }
    }
}

void copy_vector_complex(int n, double complex original[n], double complex copy[n]) {
    register int i;
    for (i = 0; i < n; i++) {
        copy[i] = original[i];
    }
}

void raise_matrix_complex(int n, int power, double complex matrix[][n], double complex raised_matrix[][n]) {
    if (power < 0) {
        printf("Power must be positive or zero ");
    }
    if (power == 0) {
        create_identity_matrix_complex(n, raised_matrix);
    } else {
        create_identity_matrix_complex(n, raised_matrix);
        register int i;
        for (i = 0; i < power; i++) {
            double complex current[n][n];
            copy_matrix_complex(n, n, raised_matrix, current);
            matrix_mult_complex(n, n, n, current, matrix, raised_matrix);
        }
    }
}

void matrix_times_scalar_complex(int rows, int cols, double complex matrix[][cols], double complex scalar) {
    register int i;
    register int j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; j++) {
            matrix[i][j] *= scalar;
        }
    }
}

void vector_times_scalar_complex(int dim, double complex vector[dim], double complex scalar) {
    register int i;
    for (i = 0; i < dim; i++) {
        vector[i] *= scalar;
    }
}

void vector_addition_complex(int dim, double complex vector1[dim], double complex vector2[dim]) {
    register int i;
    for (i = 0; i < dim; i++) {
        vector1[i] += vector2[i];
    }
}

void vector_zeroes_complex(int n, double complex vector[n]) {
    for (int i = 0; i < n; i++) {
        vector[i] = 0 + 0 * I;
    }
}

void solve_system_complex(int n, double complex A_matrix[][n], double complex B_vector[n], double complex AM[][n],
                          complex double BM[n]) {
    copy_matrix_complex(n, n, A_matrix, AM);
    copy_vector_complex(n, B_vector, BM);
    double complex scaler;
    register int fd;
    register int i;
    register int k;
    register int j;
    register double complex current_scaler;
    for (fd = 0; fd < n; fd++) {
        scaler = 1 / (AM[fd][fd]);
        for (j = 0; j < n; j++) {
            AM[fd][j] *= scaler;
        }
        BM[fd] *= scaler;
        for (i = 0; i < fd; i++) {
            current_scaler = AM[i][fd];
            for (k = 0; k < n; k++) {
                AM[i][k] = AM[i][k] - current_scaler * AM[fd][k];
            }
            BM[i] = BM[i] - current_scaler * BM[fd];
        }
        for (i = fd + 1; i < n; i++) {
            current_scaler = AM[i][fd];
            for (k = 0; k < n; k++) {
                AM[i][k] -= current_scaler * AM[fd][k];
            }
            BM[i] -= current_scaler * BM[fd];
        }
    } // Now the matrix AM is morphed into a nxn identity matrix
    // while BM contains a vector form solution of the system
}

double complex matrix_element_sum_complex(long rows, long cols, double complex matrix[][cols]) {
    register long i;
    register long j;
    register double complex sum = 0 + 0 * I;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            sum += matrix[i][j];
        }
    }
    return sum;
}

double complex vector_element_sum_complex(long n, double complex vector[n]) {
    register long i;
    register double sum = 0 + 0 * I;
    for (i = 0; i < n; i++) {
        sum += vector[i];
    }
    return sum;
}

void create_dft_matrix(int dim, double complex dft_matrix[][dim]) {
    register int i;
    register int j;
    register double complex current_freq = cexp(I * -2 * M_PI / dim);
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            dft_matrix[i][j] = cpow(current_freq, j * i + 0. * I);
        }
    }
}

void create_idft_matrix(int dim, double complex dft_matrix[][dim]) {
    register int i;
    register int j;
    register double complex current_freq = cexp(I * 2 * M_PI / dim);
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            dft_matrix[i][j] = cpow(current_freq, j * i + 0. * I);
        }
    }
}

void dft_slow(int dim, double complex f_vector[dim], double complex f_hat_vector[dim]) {
    double complex (*dft_matrix)[dim];
    dft_matrix = malloc(sizeof(double complex[dim][dim]));
    create_dft_matrix(dim, dft_matrix);
    matrix_times_vector_complex(dim, dim, dft_matrix, f_vector, f_hat_vector);
}

void idft_slow(int dim, double complex f_hat_vector[dim], double complex f_vector[dim]) {
    double complex (*idft_matrix)[dim];
    idft_matrix = malloc(sizeof(double complex[dim][dim]));
    create_idft_matrix(dim, idft_matrix);
    matrix_times_scalar_complex(dim, dim, idft_matrix, 1. / dim);
    matrix_times_vector_complex(dim, dim, idft_matrix, f_hat_vector, f_vector);
}

double complex dft_matrix[2][2] = {{1., 1.},
                                   {1., -1.}};
void fft(int dim, double complex f_vector[dim], double complex f_hat_vector[dim]) {
    if (dim <= 2) {
        matrix_times_vector_complex(2, 2, dft_matrix, f_vector, f_hat_vector);
        //dft_slow(dim,f_vector,f_hat_vector);
        return;
    }
    register long i;
    double complex *f_odd;
    double complex *f_even;
    double complex *f_odd_hat;
    double complex *f_even_hat;
    f_odd = malloc(sizeof(double complex )*dim/2);
    f_even = malloc(sizeof(double complex )* dim/2);
    f_even_hat = malloc(sizeof(double complex)*dim/2);
    f_odd_hat = malloc(sizeof(double complex)*dim/2);
    register long count_odd;
    register long count_even;
    count_odd = 0;
    count_even = 0;
    for (i = 0; i < dim; i++) {
        if (i % 2) {
            f_odd[count_odd] = f_vector[i];
            count_odd++;
        } else {
            f_even[count_even] = f_vector[i];
            count_even++;
        }
    }
    fft(dim / 2, f_even, f_even_hat);
    fft(dim / 2, f_odd, f_odd_hat);
    long count;
    count = 0;
    for (i = 0; i < dim; i++) {
        if (count == dim / 2) {
            count = 0;
        }
        f_hat_vector[i] = f_even_hat[count] + f_odd_hat[count] * cexp(-2. * I * acos(-1) * i / dim);
        count++;
    }
}

#endif //LINALG_C_ALG_H
