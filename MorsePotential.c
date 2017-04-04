
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <getopt.h>

int factorial(int n);

int main(int argc, char **argv){

    // Morse Potential: V(x) = D*(1-exp(-alpha*(x-x0)))^2
    // dynamic variables
    double D = 426.369999999928;    // kJ/mol
    double alpha = 2.511;           // 1/angstrom
    double mu = 1.0;                // g/mol
    double dx = 0.025;              // angstrom
    double xmin = -2.5;             // angstrom
    double xmax =  2.5;             // angstrom
    char *outputfile = "/dev/stdout";
    static int kcal_flag = 0; // set flag to kcal/mol
    static int evecs_flag = 0;
    static int excit_flag = 0;
    double conv_factor = 1.0;
    char *dimension = "kJ/mol";
    int numberofeigenstates = 6;

    int c;
    while(1){
        static struct option long_options[] = {
            /* These options set a flag. */
            {"kcal",  no_argument, &kcal_flag , 1},
            {"evecs", no_argument, &evecs_flag, 1},
            {"excit", no_argument, &excit_flag, 1},
            /* These options don’t set a flag.
               We distinguish them by their indices. */
            {"help",                  no_argument,       0, 'h'},
            {"xmin",                  required_argument, 0, 'a'},
            {"xmax",                  required_argument, 0, 'b'},
            {"x-spacing",             required_argument, 0, 'd'},
            {"well-depth",            required_argument, 0, 'D'},
            {"potential-width",       required_argument, 0, 'w'},
            {"reduced-mass",          required_argument, 0, 'm'},
            {"output-file",           required_argument, 0, 'o'},
            {"number-of-eigenstates", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "a:b:d:D:w:m:o:n:h", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch(c){
            case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;

            case 'h':
                printf("%s\t", argv[0]);
                printf("[--help|-h]");
                printf("[--xmin|-a <value>]");
                printf("[--xmax|-b <value>]");
                printf("\n\t\t\t");
                printf("[--x-spacing|-d <value>]");
                printf("[--well-depth|-D <value>]");
                printf("\n\t\t\t");
                printf("[--potential-width|-w <value>]");
                printf("[--output-file|-o]");
                printf("\n\t\t\t");
                printf("[--number-of-eigenstates|-n]");
                printf("[--kcal]");
                printf("\n\t\t\t");
                printf("[--reduced-mass|-m <value>]");
                printf("[--excit]");
                printf("\n\n");

                printf("-h, --help\t\t\tShow this help dialogue\n");
                printf("-a, --xmin\t\t\tSet minimum x-value\n");
                printf("-b, --xmax\t\t\tSet maximum x-value\n");
                printf("-d, --x-spacing\t\t\tSet spacing between x-values\n");
                printf("-D, --well-depth\t\tSet the D factor for well-depth in kJ/mol\n");
                printf("-w, --potential-width\t\tSet the alpha parameter for potential width\n");
                printf("-m, --reduced-mass\t\tSet reduced mass of involved particles in g/mol\n");
                printf("-o, --output-file\t\tName of output file\n");
                printf("-n, --number-of-eigenstates\tNumber of calculated eigenstates\n");
                printf("    --kcal\t\t\tOutput in kcal/mol\n");
                printf("    --evecs\t\t\tOutput eigenvectors directly (don't add evals)\n");
                printf("    --excit\t\t\tOutput excitation state energy in cm^-1\n");
                printf("\n");
                exit (0);

            case 'a':
                //printf ("option -a with value `%s'\n", optarg);
                xmin = atof(optarg);
                break;

            case 'b':
                xmax = atof(optarg);
                break;

            case 'd':
                dx = atof(optarg);
                break;

            case 'D':
                D = atof(optarg);
                break;

            case 'w':
                alpha = atof(optarg);
                break;

            case 'm':
                mu = atof(optarg);
                break;

            case 'o':
                outputfile = optarg;
                break;

            case 'n':
                numberofeigenstates = atoi(optarg);
                break;

            default:
                abort();
        }
    }

    // static declaration
    int m, n;
    double EvalTerm1, EvalTerm2;
    double E[numberofeigenstates];
    double N[numberofeigenstates];
    double powz[numberofeigenstates];
    double L[numberofeigenstates];
    double x, z, expz, lambda;
    FILE *fd = fopen(outputfile, "w");


    // constants
    double planck     = 6.626070040E-34;          // Js
    double avogadro   = 6.022140857E23;           // 1/mol
    double lightspeed = 299792458;                // m/s

    if(kcal_flag == 1){
        conv_factor = 4.184;
        dimension = "kcal/mol";
    }

    // begin calculation of Eigenvalues E
    // eval = planck*alpha/(2*pi)*sqrt(2*D/mu)       * (n + 0.5)            term1*(n+0.5)
    //      - planck*planck*alpha*alpha/(8*pi*pi*mu) * (n*(n + 1) + 0.25)   term2*(n*n+n+0.25)
    EvalTerm1 =           1.0E10/(2.0*M_PI)*avogadro * planck*alpha*sqrt(2.0*D/mu);   // kJ/mol
    EvalTerm2 = 1.0E20/(8.0*M_PI*M_PI)*avogadro*avogadro * planck*planck*alpha*alpha/mu;  // kJ/mol

    // get eigenvalues
    for(n=0; n<numberofeigenstates; ++n){
        E[n] = EvalTerm1*((double)n + 0.5) - EvalTerm2*((double)n*(double)n + (double)n + 0.25);
    }

    // begin calculation of eigenvectors Psi:
    // evec = N(n) * exp(-0.5*z) * z^(lambda-n-0.5) * L(n,x,lambda)
    lambda = 1.0/(1.0E10*avogadro) * sqrt(2.0*mu*D)/(alpha*planck/(2.0*M_PI));      // dimensionless

    for(n=0; n<numberofeigenstates; ++n){
        N[n] = sqrt(((double)factorial(n)*(2.0*lambda-2.0*(double)n-1.0))/gsl_sf_gamma(2.0*lambda - (double)n));
    }

    fprintf(fd, "# Morse Potential\n");
    fprintf(fd, "#  V(x) = D*(1-exp(-alpha*(x-x0)))^2\n");
    fprintf(fd, "# D:            %+18.12e kJ/mol\n", D);
    fprintf(fd, "# alpha:        %+18.12e 1/angstrom\n", alpha);
    fprintf(fd, "# reduced mass: %+18.12e g/mol\n", mu);
    fprintf(fd, "# EvalTerm1:    %+18.12e %s\n", EvalTerm1/conv_factor, dimension);
    fprintf(fd, "# EvalTerm2:    %+18.12e %s\n", EvalTerm2/conv_factor, dimension);
    fprintf(fd, "# lambda:       %+18.12e\n", lambda);
    for(n=0; n<numberofeigenstates; ++n){
        fprintf(fd, "#\tN%02d = % 26.20lf\t E%02d = % 26.20lf %s\n", n, N[n], n, E[n]/conv_factor, dimension);
    }

    // if excitation flag is set print excitation energies in cm^-1
    if(excit_flag == 1){
        fprintf(fd, "#\n");
        fprintf(fd, "# Excitation energies in cm^-1:\n#");
        for(n = 1; n < numberofeigenstates; ++n){
            for(m = 0; m < n; ++m){
                fprintf(fd, "\tEx%d%d = % 18.12lf", m, n, (E[n]-E[m])/avogadro/planck/lightspeed*10);
            }
            fprintf(fd, "\n#");
        }
        fprintf(fd, "\n#");
    }

    // Calculate E+Psi and output data
    fprintf(fd, "#              x                        V             ");
    for(n=0; n<numberofeigenstates; ++n){
        fprintf(fd, "         E%02d + Psi%02d        ", n, n);
    }
    fprintf(fd, "\n");
    for(x = xmin; x <= xmax; x += dx){

        z = 2.0*lambda*exp(-alpha*x); //no dimension

        expz = exp(-0.5*z);

        // z^(lambda-n-0.5)
        for(n=0; n<numberofeigenstates; ++n){
            powz[n] = pow(z, lambda - (double)n - 0.5);
        }

        // calculation of Laguerre polynomials through recursion formula:
        // Ln(z,a)=( (2*(n-1) + 1 + a - z)*L(n-1) - (n - 1 + a)*L(n-2) )/n;
        // L[0] = 1; L[1] = 1+a-z;
        // in the given case:   a = 2*lambda-2*n-1;
        L[0] = 1.0;
        L[1] = 1.0 + 2.0*lambda-2.0*1.0-1.0 - z;
        for(n=2; n<numberofeigenstates; ++n){
            L[n] = ( (2.0*((double)n-1.0) + 1.0 + 2.0*lambda-2.0*(double)n-1.0 - z)*L[n-1] - ((double)n - 1.0 + 2.0*lambda-2.0*(double)n-1.0)*L[n-2])/(double)n;
        }

        fprintf(fd, "% 26.18lf  ", x);
        fprintf(fd, "% 26.18e  ", D*(1.0 - 2.0*exp(-alpha*x) + exp(-2.0*alpha*x))/conv_factor);

        if(evecs_flag == 1){
            for(n=0; n<numberofeigenstates; ++n){
                fprintf(fd, "% 26.18lf  ", N[n] * expz * powz[n] * L[n]);
            }
        }else{
            for(n=0; n<numberofeigenstates; ++n){
                fprintf(fd, "% 26.18lf  ", E[n]/conv_factor + N[n] * expz * powz[n] * L[n]);
            }
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
