/*Nonlinear Least Squares Curve Fitting Program*/
// Reference: http://mads.lanl.gov/presentations/Leif_LM_presentation_m.pdf

/*Marquardt algorithm from P.R. Bevington,"Data Reduction and Error
Analysis for the Physical Sciences," McGraw-Hill, 1969; Adapted by
Wayne Weimer & David Harris. Jackknife error algorithm of M.S. Caceci,
Anal. Chem. 1989, 61, 2324. Translated to ANSI C by Douglas Harris &
Tim Seufert 7-94 */
#include "stdafx.h"
#include <stdio.h>
#include <Stdlib.h>
#include <math.h>
#include <ctype.h>

#define _CRT_SECURE_NO_WARNINGS
#define maxnpts 50 /* Maximum data pairs -increase if desired */

/*Change nterms to the number of parameters to be fit in your equation*/
/***********************************************************/
#define nterms 2
/* Number of parameters to be fit */
/***********************************************************/



/****** Global Variables ******/
int iteration, n, nfree;
bool use_new_params;
int npts;                               // Number of data pairs
double x[maxnpts], y[maxnpts];          // The data set to be fit
double yfit[maxnpts];                   // Calculated values of y from func()
double params[nterms];                  // The parameters at the beginning of a step before perturbation
double b[nterms];
double c[nterms];                       // The parameter used for the current function call
double final_params[nterms];
double alpha[nterms][nterms];
double **H_lm;                          // The Levenberg-Marquardt Hessian term
double **J;                             // Jacobian of residuals
double **JT;                            // Jacobian transpose
double **H;                             // The approximate Hessian, where H = JT * J
double **residual_vector;               // Vector containing the residuals
double **perturbation_vector;           // The amount by which the parameters should be adjusted for the current step
double lambda;                          // Proportion of gradient search(=0.001 at start)
double chisq;                           // Variance of residuals in curve fit
double chisq_ref_val;
const int BUFF_SIZE = 100;
char filename[20];
char answer[BUFF_SIZE];
FILE *fp;


/****** Function Declarations ******/
void readdata();
void unweightedinput();
void computeChisquare();
void computeJacobian();
void matrixInvert(double **matrix);
void curvefit();
void display();
void matrixAllocate(double ***matrix, int size_i, int size_j);
void matrixFree(double ***matrix, int size_i);
void matrixMultiply(double** A, int A_i, int A_j, double** B, int B_i, int B_j, double** C);
void matrixPrint(double **matirx, int i_size, int j_size);
void arrayPrint(double arry[], int size);
double residual(int i);
void updateResiduals();

#if defined _WIN32
errno_t err;
#endif

int main() {

    matrixAllocate(&J, maxnpts, nterms);
    matrixAllocate(&JT, nterms, maxnpts);
    matrixAllocate(&H, nterms, nterms);
    matrixAllocate(&H_lm, nterms, nterms);
    matrixAllocate(&perturbation_vector, nterms, 1);
    matrixAllocate(&residual_vector, maxnpts, 1);

    int i;
    printf("Least Squares Curve Fitting. You must modify the constant\n");
    printf("'nterms' and the fuction 'Func' for new problems.\n");
    readdata();
    printf("\n\nEnter initial guesses for parameters:\n\n");
    printf("\t(Note: Parameters cannot be exactly zero.)\n");
    for (i = 0; i < nterms; i++) {
        while(params[i] == 0.0) {
            printf("Parameter #%d =   ", i + 1);
            fgets(answer, BUFF_SIZE, stdin);
            params[i] = atof(answer);
        }
    }

    printf("\nInitial parameters:\n");
    arrayPrint(params, nterms);
    lambda = 0.001;
    iteration = 0;

    do {
        curvefit();
        iteration++;
        display();
        iteration = 0;
        printf("\n\tAnother iteration (Y/N)? ");
        fgets(answer, BUFF_SIZE, stdin);
    } while (answer[0] != 'N' && answer[0] != 'n');
    return 0;
}

// Displays the data entered
void print_data() {
    int i;
    for (i = 0; i < npts; i++) {
        printf("%d\tx = %- #12.8f\ty = %- #12.8f\t\n", i + 1, x[i], y[i]);
    }
}

                        /*******************************/
double func(double p_x) /* The function you are fitting*/
{                       /*******************************/

    // Use parameters after perturbation
    if (use_new_params) {
        for (int i = 0; i < nterms; i++) {
            c[i] = b[i];
        }
    }
    // Use parameters prior to this step's perturbation
    else {
        for (int i = 0; i < nterms; i++) {
            c[i] = params[i];
        }
    }

    /********************************************/
    /*      Enter the function to be fit:       */
    /********************************************/    
    //double value = c[0] * p_x * p_x + c[1] * p_x + c[2]; /*Ax^2 + Bx + C*/
    double value = c[0] * cos(c[1] * p_x) + c[1] * sin(c[0] * p_x);
    //printf("\nfunc(x) -- A*cos(BX) + B*sin(AX) -- x: %f a: %f  b: %f  =  value: %f\n", p_x, c[0], c[1], value);
    return ( value );
}

void readdata() {
    int n = 0;

    // Prompt for data entry type
    do {
        printf("\nDo you want to enter x,y values or read them from a file?\n");
        printf("\tType E for enter and F for File: ");
        fgets(answer, BUFF_SIZE, stdin);
        answer[0] = toupper(answer[0]);
    } while (answer[0] != 'E' && answer[0] != 'F');

    // Read from file
    if (answer[0] == 'F') {
        do {
            printf("\nPlease enter the name of the data file: ");

#if defined _WIN32
            gets_s(filename, BUFF_SIZE);
            printf("\n");
            err = fopen_s(&fp, filename, "rb");
            if (err != 0) {
                printf("Fatal error: could not open file %s\n", filename);
                exit(1);
            }
#else
            fgets(filename, BUFF_SIZE, stdin);
            printf("\n");
            fp = fopen(filename, "rb");
            if (fp == NULL) {
                printf("Fatal error: could not open file %s\n", filename);
                exit(1);
            }
#endif

            for (n = 0; !feof(fp); n++) {
                fread(&x[n], sizeof( double ), 1, fp);
                fread(&y[n], sizeof( double ), 1, fp);
            }
            fclose(fp);
            npts = n - 1;
            printf("This set contains %d points\n", npts);
            print_data();
            printf("\nIs this data correct (Y/N)?");
            fgets(answer, BUFF_SIZE, stdin);
        } while (answer[0] != 'Y' && answer[0] != 'y');
    }
    // Enter data manually
    else {
        do {
            unweightedinput();
            print_data();
            printf("Is this data correct(Y/N)?");
            fgets(answer, BUFF_SIZE, stdin);
        } while (answer[0] != 'y' && answer[0] != 'Y');
        printf("Enter name of file to save the data in: ");        
#if defined _WIN32
        gets_s(filename, BUFF_SIZE);
        err = fopen_s(&fp, filename, "wb");
        if (err != 0) {
            printf("Fatal error: could not open file %s\n", filename);
            exit(1);
        }
#else
        fgets(filename, BUFF_SIZE, stdin);
        fp = fopen(filename, "wb");
        if (fp == NULL) {
            printf("Fatal error: could not open file %s\n", filename);
            exit(1);
        }
#endif

        for (n = 0; n < npts; n++) {
            fwrite(&x[n], sizeof( double ), 1, fp);
            fwrite(&y[n], sizeof( double ), 1, fp);
        }

        fclose(fp);
        printf("Data saved in file %s\n", filename);
    }
}

/* Enter equal weight data */
void unweightedinput() {
    int i, n;
    printf("List the data in the order: x y, with one set on each\n");
    printf("line and a space (not a comma) between the numbers.\n");
    printf("Type END to end input\n");
    for (n = 0;; n++) {
        fgets(answer, BUFF_SIZE, stdin);
        if (answer[0] == 'E' || answer[0] == 'e') {
            break;
        }
        // Convert first part of string input
        x[n] = atof(answer);
        i = 0;
        while (answer[i] != ' ' && answer[i] != '\0') {
            i++;
        }
        // Convert second half of string input
        y[n] = atof(answer + i);
    }
    npts = n;
}

// Sum of square of differences between measured and calculated y values
void computeChisquare() {
    updateResiduals();
    int i;
    chisq = 0;
    for (i = 0; i < npts; i++){        
        chisq += residual(i) * residual(i);
        //printf("y[i]: %f -- yfit[i]: %f -- chisq: %f\n", y[i], yfit[i], chisq);
    }
    //chisq /= nfree;
    printf("\nchisq: %f\n\n", chisq);
}


void computeJacobian() {

    // Populate Jacobian matrix
    for (int j = 0; j < nterms; j++) {
        double param_temp = params[j];
        double delta = fabs(params[j] / 100000);
        params[j] = param_temp + delta;
        for (int i = 0; i < npts; i++) {
            printf("i: %d m: %d\n", i, j);
            J[i][j] = (func(x[i]) - yfit[i]) / delta;
        }
        params[j] = param_temp;
    }

    printf("\nJacobian matrix:\n");
    matrixPrint(J, npts, nterms);

    // Populate Jacobian transpose matrix
    for(int i = 0; i < nterms; i++){
        for(int j = 0; j < npts; j++){
            JT[i][j] = J[j][i];
        }
    }

    printf("\nJacobian transpose matrix:\n");
    matrixPrint(JT, nterms, npts);

    matrixMultiply(JT, nterms, npts, J, npts, nterms, H);

    printf("\nHessian matrix:\n");
    matrixPrint(H, nterms, nterms);
}

void matrixAllocate(double ***matrix, int size_i, int size_j){
    printf("Allocating memory for [%d] x [%d] matrix\n\n", size_i, size_j);
    *matrix = (double **)malloc(size_i * sizeof(double*));
    for(int i = 0; i < size_i; i++){
        (*matrix)[i] = (double *)malloc(size_j * sizeof( double ));
    }    
}

void matrixFree(double ***matrix, int size_j){
    for (int i = 0; i < size_j; ++i) {
      free(*matrix[i]);
    }
    free(*matrix);
}

void matrixMultiply(double** A, int A_i, int A_j, double** B, int B_i, int B_j, double** C){
    double sum = 0;

    for (int c = 0 ; c < A_i ; c++ ) {
        for (int d = 0 ; d < B_j ; d++ ) {
            for (int k = 0 ; k < B_i ; k++ ) {
                sum = sum + A[c][k] * B[k][d];
            }
            C[c][d] = sum;
            sum = 0;
        }
    }
}

// Inverts the matrix array[][]
// Pivoting reduces rounding error
void matrixInvert(double **matrix) {
    int i, j, k, ik[nterms], jk[nterms];
    double rsave, amax;

    for (k = 0; k < nterms; k++) {

        amax = 0.0;

        for (i = k; i < nterms; i++) {
            for (j = k; j < nterms; j++) {
                if (fabs(amax) <= fabs(matrix[i][j])) {
                    amax = matrix[i][j];
                    ik[k] = i;
                    jk[k] = j;
                }
            }
        }

        i = ik[k];

        if (i > k) {
            for (j = 0; j < nterms; j++) {
                rsave = matrix[k][j];
                matrix[k][j] = matrix[i][j];
                matrix[i][j] = -1 * rsave;
            }
        }

        j = jk[k];

        if (j>k) {
            for (i = 0; i < nterms; i++) {
                rsave = matrix[i][k];
                matrix[i][k] = matrix[i][j];
                matrix[i][j] = -1 * rsave;
            }
        }
        for (i = 0; i < nterms; i++) {
            if (i != k) {
                matrix[i][k] = -1 * (matrix[i][k]) / amax;
            }
        }
        for (i = 0; i < nterms; i++) {
            for (j = 0; j < nterms; j++) {
                if (j != k && i != k) {
                    matrix[i][j] = matrix[i][j] + matrix[i][k] * (matrix[k][j]);
                }
            }
        }
        for (j = 0; j < nterms; j++) {
            if (j != k) {
                matrix[k][j] = matrix[k][j] / amax;
            }
        }
        matrix[k][k] = 1 / amax;
    }
    for (k = nterms - 1; k > -1; k--) {
        j = ik[k];
        if (j > k) {
            for (i = 0; i < nterms; i++) {
                rsave = matrix[i][k];
                matrix[i][k] = -1 * (matrix[i][j]);
                matrix[i][j] = rsave;
            }
        }
        i = jk[k];
        if (i > k) {
            for (j = 0; j < nterms; j++) {
                rsave = matrix[k][j];
                matrix[k][j] = -1 * (matrix[i][j]);
                matrix[i][j] = rsave;
            }
        }
    }
}

// Curve fitting algorithm
void curvefit() {
    nfree = npts - nterms;

    double **gradient_vector;
    matrixAllocate(&gradient_vector, nterms, 1);
    
    // Clear b and gradient vectors
    for (int j = 0; j < nterms; j++) {
        b[j] = gradient_vector[j][0] = 0;
        for (int k = 0; k <= j; k++) {
            alpha[j][k] = 0;
        }
    }

    // For the first fit attempt, make func() use the old parameters
    use_new_params = false;
    
    // Find y values for current parameter values
    computeChisquare();
    chisq_ref_val = chisq;

    // Find the Jacobian
    printf("\nComputing Jacobian matrix\n");
    computeJacobian();

    // Populate the residual_vector
    for(int i = 0; i < npts; i++){
        residual_vector[i][0] = residual(i);
    }

    printf("Residual vector:\n");
    matrixPrint(residual_vector, npts, 1);

    matrixMultiply(JT, nterms, npts, residual_vector, npts, 1, gradient_vector);

    printf("Gradient vector:\n");
    matrixPrint(gradient_vector, nterms, 1);

    int loop_count = 0;
    int loop_max = 20;
    
    // Keep looping until new chisq is less than old chisq
    do {
        loop_count++;

        printf("H:\n");
        matrixPrint(H, nterms, nterms);

        // Copy the Hessian matrix into the LM Hessian, H_lm
        for (int j = 0; j < nterms; j++) {
            for (int i = 0; i < nterms; i++) {
                H_lm[i][j] = H[i][j];
            }
        }

        // Add the lambda*diag(Hessian) terms
        for (int i = 0; i < nterms; i++) {
            H_lm[i][i] += H[i][i] * lambda;
        }

        printf("H_lm:\n");
        matrixPrint(H_lm, nterms, nterms);

        matrixInvert(H_lm);

        printf("Inverse H_lm\n");
        matrixPrint(H_lm, nterms, nterms);

        // Get the difference vector, delta = (H + lambda*Diag(H))^-1) * gradient
        matrixMultiply(H_lm, nterms, nterms, gradient_vector, nterms, 1, perturbation_vector);

        printf("Parameter offset\n");
        matrixPrint(perturbation_vector, nterms, 1);

        // Generate new parameters from perturbation vector
        for (int i = 0; i < nterms; i++) {
            b[i] = params[i] - perturbation_vector[i][0];
        }

        printf("New parameters\n");
        arrayPrint(b, nterms);
        
        // Find the new chi squared value using the new parameters
        use_new_params = true;        
        computeChisquare();
        
        // If this step wasn't successful, try again with a higher lambda value
        if (( chisq_ref_val - chisq ) < 0) {
            lambda *= 10;
            printf("Failed to lower chisq, trying bigger lambda: %f\n", lambda);
        }
    // Loop until the new parameters produce a lower chi squared value than we started with
    } while(chisq > chisq_ref_val && loop_count < loop_max);

    if(loop_count < loop_max){
        printf("Step accepted, saving new params, decreasing lambda value\n");
    }
    else{
        printf("Failed to converge. :-(");
    }
    printf("Old chisq: %f  new chisq: %f", chisq_ref_val, chisq);

    // Update the parameters list with the new parameter set
    for (int j = 0; j < nterms; j++) {
        params[j] = b[j];
    }

    // Step is accepted
    lambda /= 10;

    matrixFree(&gradient_vector, 1);
}

void updateParameters() {

}

void updateResiduals(){
    for (int i = 0; i < npts; i++) {
        yfit[i] = func(x[i]);
    }
}

double residual(int i){
    return y[i] - yfit[i];
}

// Prints result of curve fit
void display() {
    int i;
    printf("\nIteration #%d\n", iteration);
    for (i = 0; i < nterms; i++) {
        printf("Params[%3dl = %-#12.8f\n", i, params[i]);
        final_params[i] = params[i];
    }
    printf("Sum of squares of residuals = %- #12.8f", chisq);
}   

void matrixPrint(double **matrix, int i_size, int j_size) {
    for (int i = 0; i < i_size; i++) {
        for (int j = 0; j < j_size; j++) {
            printf("%f, ", matrix[i][j]);
        }
        printf("\n\n");
    }
}

void arrayPrint(double this_array[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%f, ", this_array[i]);  
    }
    printf("\n\n");
}
