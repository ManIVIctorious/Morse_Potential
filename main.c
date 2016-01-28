
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>

int main(int argc, char **argv){

    // Morse Potential: V(x) = D*(1-exp(-alpha*(x-x0)))^2
    // dynamic declaration
    double D = 426.369999999928;    // kJ/mol
    double alpha = 2.511;           // 1/angstrom
    double mu = 1;
    double dx = 0.025;
    double xmin = -2.5;
    double xmax =  2.5;

    // static declaration
    double EvalTerm1, EvalTerm2;
    double E0, E1, E2, E3, E4, E5;
    double N0, N1, N2, N3, N4, N5, lambda;
    double z, x;
    double L0, L1, L2, L3, L4, L5;

    // constants
    double pi       = 3.14159265358979323846;
    double planck   = 6.626070040E-34;          // Js
    double hbar     = 1.05457180013E-34;        // Js           planck/(2*pi)
    double avogadro = 6.02214085774E23;         // 1/mol

    // begin calculation of Eigenvalues E
    // eval = planck*alpha/(2*pi)*sqrt(2*D/mu)       * (n + 0.5)            term1*(n+0.5)
    //      - planck*planck*alpha*alpha/(8*pi*pi*mu) * (n*(n + 1) + 0.25)   term2*(n*n+n+0.25)
    EvalTerm1 =           1.0E10/(2.0*pi)*avogadro * planck*alpha*sqrt(2.0*D/mu);   // kJ/mol
    EvalTerm2 = 1.0E20/(8*pi*pi)*avogadro*avogadro * planck*planck*alpha*alpha/mu;  // kJ/mol

    // get eigenvalues
    E0 = EvalTerm1*(0.5)-EvalTerm2*( 0.25);
    E1 = EvalTerm1*(1.5)-EvalTerm2*( 2.25);
    E2 = EvalTerm1*(2.5)-EvalTerm2*( 6.25);
    E3 = EvalTerm1*(3.5)-EvalTerm2*(12.25);
    E4 = EvalTerm1*(4.5)-EvalTerm2*(20.25);
    E5 = EvalTerm1*(5.5)-EvalTerm2*(30.25);

    // begin calculation of eigenvectors Psi:
    // evec = N(n) * exp(-0.5*z) * z^(lambda-n-0.5) * L(n,x,lambda)
    //
    // dimensionless terms
    lambda = 1/(1.0E10*avogadro) * sqrt(2*mu*D)/(alpha*hbar);         //printf("lambda = %lf\n", lambda);
    N0 = sqrt((  2.0*lambda -    1.0)/gsl_sf_gamma(2.0*lambda));      //printf("N0 = %lf\n", N0);
    N1 = sqrt((  2.0*lambda -    3.0)/gsl_sf_gamma(2.0*lambda - 1));  //printf("N1 = %lf\n", N1);
    N2 = sqrt((  4.0*lambda -   10.0)/gsl_sf_gamma(2.0*lambda - 2));  //printf("N2 = %lf\n", N2);
    N3 = sqrt(( 12.0*lambda -   42.0)/gsl_sf_gamma(2.0*lambda - 3));  //printf("N3 = %lf\n", N3);
    N4 = sqrt(( 48.0*lambda -  216.0)/gsl_sf_gamma(2.0*lambda - 4));  //printf("N4 = %lf\n", N4);
    N5 = sqrt((240.0*lambda - 1320.0)/gsl_sf_gamma(2.0*lambda - 5));  //printf("N5 = %lf\n", N5);

    double evecterm2; 
    double evecterm30, evecterm31, evecterm32, evecterm33, evecterm34, evecterm35;
    for(x = xmin; x <= xmax; x += dx){
        
        z = 2*lambda*exp(-alpha*x); //no dimension
        
        evecterm2 = exp(-0.5*z);
    
        evecterm30 = pow(z, lambda-0.5);
        evecterm31 = pow(z, lambda-1.5);
        evecterm32 = pow(z, lambda-2.5);
        evecterm33 = pow(z, lambda-3.5);
        evecterm34 = pow(z, lambda-4.5);
        evecterm35 = pow(z, lambda-5.5);

        L0 = 1;
        L1 = 2*lambda - 2 - z;
        //Ln=( (2*(n-1) + 1 + 2*lambda-2*n-1 - z)*L(n-1) - (n - 1 + 2*lambda-2*n-1)*L(n-2) )/n;
        L2 = ( (2*(2-1) + 1 + 2*lambda-2*2-1 - z)*L1 - (2 - 1 + 2*lambda-2*2-1)*L0 )/2;
        L3 = ( (2*(3-1) + 1 + 2*lambda-2*3-1 - z)*L2 - (3 - 1 + 2*lambda-2*3-1)*L1 )/3;
        L4 = ( (2*(4-1) + 1 + 2*lambda-2*4-1 - z)*L3 - (4 - 1 + 2*lambda-2*4-1)*L2 )/4;
        L5 = ( (2*(5-1) + 1 + 2*lambda-2*5-1 - z)*L4 - (5 - 1 + 2*lambda-2*5-1)*L3 )/5;

        printf("%lf  ", x);
        printf("%lf  ", D*(1.0 - 2.0*exp(-alpha*x) + exp(-2.0*alpha*x)));
        printf("%lf  ", E0 + N0 * evecterm2 * evecterm30 * L0);
        printf("%lf  ", E1 + N1 * evecterm2 * evecterm31 * L1);
        printf("%lf  ", E2 + N2 * evecterm2 * evecterm32 * L2);
        printf("%lf  ", E3 + N3 * evecterm2 * evecterm33 * L3);
        printf("%lf  ", E4 + N4 * evecterm2 * evecterm34 * L4);
        printf("%lf  ", E5 + N5 * evecterm2 * evecterm35 * L5);
        printf("\n");                                   
    }

    return 0;
}
