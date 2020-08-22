//
// Created by User on 22/8/2020.
//

#ifndef LINALG_ODE_H
#define LINALG_ODE_H

void copy_vector_ode(int dim, double original[dim],double copy[dim]){
    for(int i = 0; i < dim ;i++){
        copy[i] = original[i];
    }
}

void ode_45_single_step(int dim  , double *yk, double tk , void (*f)(double , double * , double *) , double dt ) {
    double f1[dim];
    double f2[dim];
    double f3[dim];
    double f4[dim];
    double current[dim];
    register int i ;
    f(tk,yk,f1);
    for( i = 0 ; i < dim ; i++){
        current[i] = yk[i] + dt/2 * f1[i];
    }
    f(tk+dt/2,current,f2);
    for( i = 0 ; i < dim ; i++){
        current[i] = yk[i] + dt/2 * f2[i];
    }
    f(tk+dt/2,current,f3);
    for( i = 0 ; i < dim ; i ++){
        current[i] = yk[i] + dt/2 * f3[i];
    }
    f(tk+dt , current , f4);
    for( i = 0 ; i < dim ; i++){
        yk[i] += dt/6 * (f1[i] + 2*(f2[i]+f3[i]+f4[i]));
    }

}


void ode_45_solve(int dim , void  (*f)(double , double * , double *)  , double *initial_condition, double time_zero , double time_final , double time_step , double **state_matrix){
    register int i ;
    register long count;
    count = 0;
    double yk[dim];
    double tk = time_zero;
    copy_vector_ode(dim,initial_condition,yk);
    while(tk <= time_final){
        for(i = 0 ; i < dim ; i++){
            state_matrix[count][i] = yk[i];
        }
        ode_45_single_step(dim, yk, tk, f, time_step);
        count++;
        tk += time_step;
    }
}

void ode_78_single_step(int dim  , double *yk, double tk , void (*f)(double , double * , double *) , double dt ) {
    double f1[dim];
    double f2[dim];
    double f3[dim];
    double f4[dim];
    double f5[dim];
    double f6[dim];
    double current[dim];
    register int i ;
    f(tk,yk,f1);
    for(i = 0 ; i < dim ; i++){
        current[i] = dt/2 * f1[i] + yk[i];
    }
    f(tk+dt/2,current,f2);
    for(i = 0 ; i < dim ; i ++){
        current[i] = yk[i] + dt/4 * f1[i] + dt/4 * f2[i];
    }
    f(tk+dt/2,current,f3);
    for(i = 0 ; i < dim ; i++){
        current[i]  = yk[i] - dt*f2[i] + 2.*dt * f3[i];
    }
    f(tk+dt,current,f4);
    for(i = 0 ; i < dim ; i++){
        current[i] = yk[i] + dt * 7./27 *f1[i] + dt*10./27 * f2[i] + dt/27 * f4[i];
    }
    f(tk+dt*2./3,current,f5);
    for(i = 0 ; i < dim ; i++){
        current[i] = yk[i] + dt*28./625*f1[i] - dt/5 * f2[i] + dt*546./625 * f3[i] +dt*54./625*f4[i] -dt*378./625*f5[i];
    }
    f(tk + dt/5,current,f6);
    for(i = 0 ; i < dim ; i++){
        yk[i] += dt/24 * f1[i] + dt*5/48 * f4[i] + dt*27/56 * f5[i] + dt*125/336*f6[i];
    }
}
void ode_78_solve(int dim , void  (*f)(double , double * , double *)  , double *initial_condition, double time_zero , double time_final , double time_step , double state_matrix[][dim]){
    register int i ;
    register long count;
    count = 0;
    double yk[dim];
    double tk = time_zero;
    copy_vector_ode(dim,initial_condition,yk);
    while(tk <= time_final){
        for(i = 0 ; i < dim ; i++){
            state_matrix[count][i] = yk[i];
        }
        ode_78_single_step(dim, yk, tk, f, time_step);
        count++;
        tk += time_step;
    }


}
#endif //LINALG_ODE_H
