//
//  main.c
//  MC simple option
//
//  Created by Leon Kwok on 24/10/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double sigma, r, E, T, M;
double ds, dt, sd;

int main() {
    clock_t begin, end;
    begin = clock();
    
    double u1, u2;
    double x, y, sum, smin, smax;
    double MC1(double, double);
    double MC2(double, double);
    double kmax;
    double exact_price_function(double);       //exact solution of option price
    double error;
    
    srand(time(NULL));
    
    FILE *op;
    op=fopen("option_price.csv", "w");
    
    r= 0.05;        //interest rate
    sigma= 0.2;     //volatility
    E= 10;          //exercise price
    T= 0.5;         //6 months but T is in the unit of years
    sd=sigma*sqrt(T);
    ds= 0.01;
    smax= 11;
    smin= 9;
    kmax=(smax-smin)/ds;
    dt= 0.0005;
    M=T/dt;
    printf("M=%lf\n", M);
    int N=100000;
    
    //printf("x=%lf\t y=%lf\n", x, y);
    double *S = malloc(N*sizeof(double));
    double *C = malloc(round(kmax)*sizeof(double));
    double *exact_C = malloc(round(kmax)*sizeof(double));
    for (int k=0; k<kmax+1; k++) {
        for (int i=0; i<N; i++) {
            S[i]=(smin)+k*ds;
                for (int j=0; j<M; j++) {
                    u1= (double)rand()/(double)RAND_MAX;
                    u2= (double)rand()/(double)RAND_MAX;
                    x=sqrt(-2*log(u1))*cos(2*M_PI*u2);
                    S[i]=MC1(S[i], x);
                }
        }
        sum=0;
        for (int i=0; i<N; i++) {
            if (S[i]>=E) {
                S[i]=(S[i]-E)*exp(-r*T);
                sum+=S[i];
            }
        }
        C[k]=(sum/N);
        exact_C[k]=exact_price_function(smin+k*ds);
        error=(C[k]-exact_C[k])/exact_C[k]*100;
        //printf("C=%lf\t exact_C=%lf\t error=%lf\n", C[k], exact_C[k], error);
        fprintf(op, "%lf\t %lf\t %lf\n", smin+k*ds, C[k], exact_C[k]);
    }
    
    free(S);
    fclose(op);
    
    end = clock();
    printf ("time used to evaluate program: %f microsecs\n", (double)(end-begin)/(double)CLOCKS_PER_SEC);
    return 0;
}

double MC1(double S0, double dw){
    double y;
    y=S0*exp((r-0.5*sigma*sigma)*dt+sigma*sqrt(dt)*dw);
    return y;
}

double MC2(double S0, double dw){
    double y, ds;
    ds= S0*(r*dt+sigma*sqrt(dt)*dw);
    y=S0+ds;
    return y;
}

double exact_price_function(double S){
    double y, d1, d2;       //integrate from negative infinity to x
    d1= (log(S/E)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));       //log in c is natural log
    d2= (log(S/E)+(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    y= S*0.5*(erf(d1/sqrt(2))+1)-E*exp(-r*T)*0.5*(erf(d2/sqrt(2))+1);       //N() function has 1/sqrt(2pi), so y is different from erf()
    return y;
}
