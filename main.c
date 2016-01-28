
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>

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

    // static declaration
    double EvalTerm1, EvalTerm2;
    double E0, E1, E2, E3, E4, E5;

    double N0, N1, N2, N3, N4, N5;
    double x, z, expz, lambda;
    double powz0, powz1, powz2, powz3, powz4, powz5;
    double L0, L1, L2, L3, L4, L5;
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

    fprintf(fd, "#EvalTerm1 = %lf kJ/mol\n", EvalTerm1);
    fprintf(fd, "#EvalTerm2 = %lf kJ/mol\n", EvalTerm2);
    fprintf(fd, "#lambda = %lf\n", lambda);
    fprintf(fd, "#N0 = % 22.20lf\n", N0);
    fprintf(fd, "#N2 = % 22.20lf\n", N2);
    fprintf(fd, "#N3 = % 22.20lf\n", N3);
    fprintf(fd, "#N4 = % 22.20lf\n", N4);
    fprintf(fd, "#N5 = % 22.20lf\n", N5);

    // Calculate E+Psi and output data
    //fprintf(fd, "# x\tV(x)\tE0+Psi0\tE1+Psi1\tE2+Psi2\tE3+Psi3\tE4+Psi4\tE5+Psi5\n");
    fprintf(fd, "#  x\t\t");
    fprintf(fd, " V(x)\t\t\t");
    fprintf(fd, " E0+Psi0\t");
    fprintf(fd, " E1+Psi1\t");
    fprintf(fd, " E2+Psi2\t");
    fprintf(fd, " E3+Psi3\t");
    fprintf(fd, " E4+Psi4\t");
    fprintf(fd, " E5+Psi5\t");
    fprintf(fd, "\n");
    for(x = xmin; x <= xmax; x += dx){
        
        z = 2*lambda*exp(-alpha*x); //no dimension
        
        expz = exp(-0.5*z);
    
        // z^(lambda-n-0.5)
        powz0 = pow(z, lambda-0.5);
        powz1 = pow(z, lambda-1.5);
        powz2 = pow(z, lambda-2.5);
        powz3 = pow(z, lambda-3.5);
        powz4 = pow(z, lambda-4.5);
        powz5 = pow(z, lambda-5.5);

        // calculation of Laguerre polynomials through recursion formula: 
        // Ln(z,a)=( (2*(n-1) + 1 + a - z)*L(n-1) - (n - 1 + a)*L(n-2) )/n;
        // L0 = 1; L1 = 1+a-z;
        // in the given case:   a = 2*lambda-2*1-1;
        L0 = 1;
        L1 = 1 + 2*lambda-2*1-1 - z;
        L2 = ( (2*(2-1) + 1 + 2*lambda-2*2-1 - z)*L1 - (2 - 1 + 2*lambda-2*2-1)*L0 )/2;
        L3 = ( (2*(3-1) + 1 + 2*lambda-2*3-1 - z)*L2 - (3 - 1 + 2*lambda-2*3-1)*L1 )/3;
        L4 = ( (2*(4-1) + 1 + 2*lambda-2*4-1 - z)*L3 - (4 - 1 + 2*lambda-2*4-1)*L2 )/4;
        L5 = ( (2*(5-1) + 1 + 2*lambda-2*5-1 - z)*L4 - (5 - 1 + 2*lambda-2*5-1)*L3 )/5;

        fprintf(fd, "% 8.8lf\t", x);
        fprintf(fd, "% 19.8lf\t", D*(1.0 - 2.0*exp(-alpha*x) + exp(-2.0*alpha*x)));
        fprintf(fd, "% 15.8lf\t", E0 + N0 * expz * powz0 * L0);
        fprintf(fd, "% 15.8lf\t", E1 + N1 * expz * powz1 * L1);
        fprintf(fd, "% 15.8lf\t", E2 + N2 * expz * powz2 * L2);
        fprintf(fd, "% 15.8lf\t", E3 + N3 * expz * powz3 * L3);
        fprintf(fd, "% 15.8lf\t", E4 + N4 * expz * powz4 * L4);
        fprintf(fd, "% 15.8lf\t", E5 + N5 * expz * powz5 * L5);
        fprintf(fd, "\n");                                   
    }

    return 0;
}
