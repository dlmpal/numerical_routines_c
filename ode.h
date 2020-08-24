//
// Created by User on 22/8/2020.
//

#ifndef LINALG_ODE_H
#define LINALG_ODE_H

//TIME STEP : IF DESIRED TIME STEP IS LESS THAN OR EQUAL TO 0.001 THE STATE MATRIX SHOULD
// DECLARED INSIDE OF MAIN AND NOT MALLOCED
// IF LARGER THEN SAID MATRIX MUST BE DECLARED OUTSIDE OF MAIN
// AND MALLOCED INSIDE
// RECOMMENDED THAT TIME STEP IS A MULTIPLE OF 10

void copy_vector_ode(int dim, double original[dim],double copy[dim]);

void ode_45_single_step(int dim  , double *yk, double tk , void (*f)(double , double * , double *) , double dt );

void ode_45_solve(int dim , void  (*f)(double , double * , double *)  , double *initial_condition, double time_zero , double time_final , double time_step , double state_matrix[][dim]);

void ode_78_single_step(int dim  , double *yk, double tk , void (*f)(double , double * , double *) , double dt );

void ode_78_solve(int dim , void  (*f)(double , double * , double *)  , double *initial_condition, double time_zero , double time_final , double time_step , double state_matrix[][dim]);

#endif //LINALG_ODE_H
