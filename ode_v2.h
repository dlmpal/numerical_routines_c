//
// Created by User on 19/8/2020.
//

#ifndef LINALG_ODE_V2_H
#define LINALG_ODE_V2_H
#include <stdlib.h>
#include "linalg.h"

double *rk78_single_step(int dim , double yk[dim] , double tk , void (*f)(double , double *  , double *),double dt) {
    register double *current_state;
    register double *current_state2;
    double *final_state;
    register double* f1;
    register double* f2;
    register double* f3;
    register double* f4;
    register double* f5;
    register double* f6;
    register int i ;
    current_state = (double*) malloc(dim*sizeof(double));
    current_state2 = (double*) malloc(dim* sizeof(double));
    final_state = (double*) malloc(dim* sizeof(double));
    f1 = (double*) malloc(dim* sizeof(double));
    f2 = (double*) malloc(dim* sizeof(double));
    f3 = (double*) malloc(dim* sizeof(double));
    f4 = (double*) malloc(dim* sizeof(double));
    f5 = (double*) malloc(dim* sizeof(double));
    f6 = (double*) malloc(dim* sizeof(double));
    //
    f(tk,yk,f1);
    vector_times_scalar(dim,f1,dt);
    //
    copy_vector(dim,f1,current_state);
    vector_times_scalar(dim,current_state,.5);
    vector_addition(dim,current_state,yk);
    f(tk + .5*dt,current_state,f2);
    vector_times_scalar(dim,f2,dt);
    //
    copy_vector(dim,f2,current_state);
    vector_addition(dim,current_state,f1);
    vector_times_scalar(dim,current_state,1./4);
    vector_addition(dim,current_state,yk);
    f(tk+.5*dt,current_state,f3);
    vector_times_scalar(dim,f3,dt);
    //
    copy_vector(dim,f3,current_state);
    vector_times_scalar(dim,current_state,2.);
    copy_vector(dim,f2,current_state2);
    vector_times_scalar(dim,current_state2,-1.);
    vector_addition(dim,current_state,current_state2);
    vector_addition(dim,current_state,yk);
    f(tk+dt,current_state,f4);
    vector_times_scalar(dim,f4,dt);
    //
    copy_vector(dim,f4,current_state);
    vector_times_scalar(dim,current_state,1./27);
    copy_vector(dim,f2,current_state2);
    vector_times_scalar(dim,current_state2,10./27);
    vector_addition(dim,current_state,current_state2);
    copy_vector(dim,f1,current_state2);
    vector_times_scalar(dim,current_state2,7./27);
    vector_addition(dim,current_state,current_state2);
    vector_addition(dim,current_state,yk);
    f(tk + (2./3)*dt , current_state,f5);
    vector_times_scalar(dim,f5,dt);
    //
    copy_vector(dim,f5,current_state);
    vector_times_scalar(dim,current_state,-378./625.);
    copy_vector(dim,f4,current_state2);
    vector_times_scalar(dim,current_state2,54./625);
    vector_addition(dim,current_state,current_state2);
    copy_vector(dim,f3,current_state2);
    vector_times_scalar(dim,current_state2,546./625);
    vector_addition(dim,current_state,current_state2);
    copy_vector(dim,f2,current_state2);
    vector_times_scalar(dim,current_state2,-.2);
    vector_addition(dim,current_state,current_state2);
    copy_vector(dim,f1,current_state2);
    vector_times_scalar(dim,current_state2,28./625);
    vector_addition(dim,current_state,current_state2);
    vector_addition(dim,current_state,yk);
    f(tk+.2*dt,current_state,f6);
    vector_times_scalar(dim,f6,dt);
    //
    for(i = 0 ; i < dim ; i++){
        final_state[i] = yk[i] + 1./24 * f1[i] + 5./48 * f4[i] + 27./56 * f5[i] + 125./336 * f6[i]  ;
    }
    free(current_state);
    free(current_state2);
    free(f1);
    free(f2);
    free(f3);
    free(f4);
    free(f5);
    free(f6);
    return final_state;
}

void rk78_solve(int dim, void (*f)(double, double *,double *), double initial_cond[dim], double time_zero, double time_final,
               double time_step, double state_matrix[][dim]) {
    double *yk;
    yk = (double *) malloc(dim * sizeof(double));
    yk = initial_cond;
    double tk = time_zero;
    register long count = 0;
    register int i;
    while (tk <= time_final) {
       for(i = 0 ; i < dim ; i++){
           state_matrix[count][i] = yk[i];
       }
        yk = rk78_single_step(dim, yk, tk, f, time_step);
        tk += time_step;
        count++;
    }
    free(yk);
}
double *rk45_single_step(int dim , double yk[dim] , double tk , void (*f)(double , double *  , double *),double dt) {
    register double *current_state;
    current_state = (double *) malloc(dim * sizeof(double));
    register double *current_state_copy;
    current_state_copy = (double *) malloc(dim * sizeof(double));
    double *final_state;
    final_state = (double *) malloc(dim * sizeof(double));
    copy_vector(dim, yk, final_state);
    f(tk, yk, current_state); // f1
    copy_vector(dim, current_state, current_state_copy);
    vector_times_scalar(dim, current_state_copy, dt / 6);
    vector_addition(dim, final_state, current_state_copy);
    vector_times_scalar(dim, current_state, dt / 2);
    vector_addition(dim, current_state, yk);
    copy_vector(dim, current_state, current_state_copy);
    f(tk + dt / 2, current_state_copy, current_state); // f2
    copy_vector(dim, current_state, current_state_copy);
    vector_times_scalar(dim, current_state_copy, dt / 3);
    vector_addition(dim, final_state, current_state_copy);
    vector_times_scalar(dim, current_state, dt / 2);
    vector_addition(dim, current_state, yk);
    copy_vector(dim, current_state, current_state_copy);
    f(tk + dt / 2, current_state_copy, current_state); //f3
    copy_vector(dim, current_state, current_state_copy);
    vector_times_scalar(dim, current_state_copy, dt / 3);
    vector_addition(dim, final_state, current_state_copy);
    vector_times_scalar(dim, current_state, dt / 2);
    vector_addition(dim, current_state, yk);
    copy_vector(dim, current_state, current_state_copy); //f4
    f(tk + dt, current_state_copy, current_state);
    vector_times_scalar(dim, current_state, dt / 3);
    vector_addition(dim, final_state, current_state);
    free(current_state);
    free(current_state_copy);
    return final_state;
}

void rk45_solve(int dim, void (*f)(double, double *, double *), double *initial_cond, double time_zero, double time_final,
                double time_step, double state_matrix[][dim]) {
        double *yk;
        yk = (double *) malloc(dim * sizeof(double));
        yk = initial_cond;
        double tk = time_zero;
        register long count = 0;
        register int i ;
        while (tk <= time_final) {
            for(i = 0 ; i < dim ; i++){
                state_matrix[count][i] = yk[i];
            }
            yk = rk45_single_step(dim, yk, tk, f, time_step);
            tk += time_step;
            count++;
        }
        free(yk);
    }
#endif //LINALG_ODE_V2_H
