
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>

int factorial(int n);

int main(int argc, char **argv){

    // Morse Potential: V(x) = D*(1-exp(-alpha*(x-x0)))^2
    // dynamic variables
    double D = 426.369999999928;    // kJ/mol
    double alpha = 2.511;           // 1/angstrom
    double mu = 1;                  // g/mol
    double dx = 0.025;              // angstrom
    double xmin = -2.5;             // angstrom
    double xmax =  2.5;             // angstrom
    char *outputfile = "/dev/stdout";
    int numberofeigenstates=6;

    // static declaration
    int n;
    double EvalTerm1, EvalTerm2;
    double E[numberofeigenstates];
    double N[numberofeigenstates];
    double powz[numberofeigenstates];
    double L[numberofeigenstates];
    double x, z, expz, lambda;
    FILE *fd = fopen(outputfile, "w");

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
    for(n=0; n<numberofeigenstates; ++n){
        E[n] = EvalTerm1*((double)n + 0.5) - EvalTerm2*((double)n*(double)n + (double)n + 0.25);
    }

    // begin calculation of eigenvectors Psi:
    // evec = N(n) * exp(-0.5*z) * z^(lambda-n-0.5) * L(n,x,lambda)
    lambda = 1/(1.0E10*avogadro) * sqrt(2*mu*D)/(alpha*hbar);        // dimensionless 

    for(n=0; n<numberofeigenstates; ++n){
        N[n] = sqrt(((double)factorial(n)*(2.0*lambda-2.0*(double)n-1.0))/gsl_sf_gamma(2.0*lambda - (double)n));
    }

    fprintf(fd, "#EvalTerm1 = %lf kJ/mol\n", EvalTerm1);
    fprintf(fd, "#EvalTerm2 = %lf kJ/mol\n", EvalTerm2);
    fprintf(fd, "#lambda = %lf\n", lambda);
    for(n=0; n<numberofeigenstates; ++n){
        fprintf(fd, "#N%d = % 26.20lf\t E%d = % 26.20lf\n", n, N[n], n, E[n]);
    }

    // Calculate E+Psi and output data
    //fprintf(fd, "# x\tV(x)\tE0+Psi0\tE1+Psi1\t...En+Psin\n");
    fprintf(fd, "#  x\t\t");
    fprintf(fd, " V(x)\t\t\t");
    for(n=0; n<numberofeigenstates; ++n){
        fprintf(fd, " E%d+Psi%d\t", n,n);
    }
    fprintf(fd, "\n");
    for(x = xmin; x <= xmax; x += dx){
        
        z = 2*lambda*exp(-alpha*x); //no dimension
        
        expz = exp(-0.5*z);
    
        // z^(lambda-n-0.5)
        for(n=0; n<numberofeigenstates; ++n){
            powz[n] = pow(z, lambda - (double)n - 0.5);
        }

        // calculation of Laguerre polynomials through recursion formula: 
        // Ln(z,a)=( (2*(n-1) + 1 + a - z)*L(n-1) - (n - 1 + a)*L(n-2) )/n;
        // L[0] = 1; L[1] = 1+a-z;
        // in the given case:   a = 2*lambda-2*n-1;
        L[0] = 1;
        L[1] = 1 + 2*lambda-2*1-1 - z;
        for(n=2; n<numberofeigenstates; ++n){
            L[n] = ( (2*((double)n-1) + 1 + 2*lambda-2*(double)n-1 - z)*L[n-1] - ((double)n - 1 + 2*lambda-2*(double)n-1)*L[n-2])/(double)n;
        }

        fprintf(fd, "% 8.8lf\t", x);
        fprintf(fd, "% 19.8lf\t", D*(1.0 - 2.0*exp(-alpha*x) + exp(-2.0*alpha*x)));
        for(n=0; n<numberofeigenstates; ++n){
            fprintf(fd, "% 15.8lf\t", E[n] + N[n] * expz * powz[n] * L[n]);
        }
        fprintf(fd, "\n");                                   
    }
    fclose(fd);

    return 0;
}

int factorial(int n){
    int fact = 1;

    while(n > 1){
        fact = fact * n;
        n = n - 1;
    }
    return fact;
}
