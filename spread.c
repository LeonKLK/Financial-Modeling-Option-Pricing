//
//  main.c
//  MC spread option
//
//  Created by Leon Kwok on 11/11/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double sigma1, sigma2, r, E, T, M, rho;
double ds, dt, sd;

int main() {
    clock_t begin, end;
    begin = clock();
    
    double u1, u2;
    double x, y, sum, iS1, iS2;
    double MC1(double, double, double);
    double MC2(double, double, double);
    double kmax;
    double kirk_function1(double, double, double, double);       //kirk solution of option price
    double kirk_function2(double, double, double, double);
    double error1, error2;
    
    srand(time(NULL));
    
    FILE *op;
    op=fopen("option_price.csv", "w");
    
    r= 0.05;
    sigma1=0.3;
    sigma2=0.2;
    E= 20;
    rho=0.6;
    T= 1;
    ds= 0.5;
    kmax=20/ds;
    dt= 0.0005;
    M=T/dt;
    printf("M=%lf\n", M);
    int N=100000;
    
    double *S1 = malloc(N*sizeof(double));
    double *S2 = malloc(N*sizeof(double));
    double *C = malloc(round(kmax)*sizeof(double));
    double *exact_C1 = malloc(round(kmax)*sizeof(double));
    double *exact_C2 = malloc(round(kmax)*sizeof(double));
    iS1=100;
    iS2=90;
    for (int k=0; k<kmax+1; k++) {
        for (int i=0; i<N; i++) {
            S1[i]=iS1+k*ds;
            S2[i]=iS2;
                for (int j=0; j<M; j++) {
                    u1= (double)rand()/(double)RAND_MAX;
                    u2= (double)rand()/(double)RAND_MAX;
                    x=sqrt(-2*log(u1))*cos(2*M_PI*u2);
                    y=rho*x+sqrt(1-rho*rho)*sqrt(-2*log(u1))*sin(2*M_PI*u2);
                    S1[i]=MC1(S1[i], x, sigma1);
                    S2[i]=MC1(S2[i], y, sigma2);
                }
            //printf("S1=%lf S2=%lf\n", S1[i], S2[i]);
            //printf("diff=%lf\n", S2[i]-S1[i]);
        }
        sum=0;
        for (int i=0; i<N; i++) {
            if ((S1[i]-S2[i])>=E) {
                sum = sum+(S1[i]-S2[i]-E)*exp(-r*T);
                //printf("S2=%lf S1=%lf diff=%lf\n", S2[i], S1[i], S2[i]-S1[i]);
                //printf("sum=%lf\n", sum);
            }
        }
        //printf("%lf\n", sum);
        C[k]=(sum/N);
        //exact_C1[k]=kirk_function1(iS1, iS2+k*ds, sigma1, sigma2);
        exact_C2[k]=kirk_function2(iS1+k*ds, iS2, sigma1, sigma2);
        //error1=(C[k]-exact_C1[k])/exact_C1[k]*100;
        error2=(C[k]-exact_C2[k])/exact_C2[k]*100;
        printf("C=%lf\t exact_C2=%lf\t error2=%lf\n", C[k], exact_C2[k], error2);
        fprintf(op, "%lf\t %lf\t %lf\t %lf\t %lf\n",iS1, iS2+k*ds, C[k], exact_C2[k], error2);
        //fprintf(op, "%lf\t %lf\t %lf\n", (E-1)+k*ds, C[k], exact_C[k]);
    }
    
    free(S1);
    free(S2);
    free(C);
    free(exact_C1);
    //free(exact_C2);
    fclose(op);
    
    end = clock();
    printf ("time used to evaluate program: %f secs\n", (double)(end-begin)/(double)CLOCKS_PER_SEC);
    return 0;
}

double MC1(double S0, double dw, double sigma){
    double y;
    y=S0*exp((r-0.5*sigma*sigma)*dt+sigma*sqrt(dt)*dw);
    return y;
}

double MC2(double S0, double dw, double sigma){
    double y, ds;
    ds= S0*(r*dt+sigma*sqrt(dt)*dw);
    y=S0+ds;
    return y;
}

double kirk_function1(double S1, double S2, double sigma1, double sigma2){
    double y, d1, d2, sigmak;
    sigmak=sqrt(sigma1*sigma1-2*S2*rho*sigma1*sigma2/(S2+E)+(S2/(S2+E))*(S2/(S2+E))*sigma2*sigma2);
    d1= (log(S1/(S2+E))+0.5*sigmak*sigmak*T)/(sigmak*sqrt(T));
    d2= d1-sigmak*sqrt(T);
    y= exp(-r*T)*(S1*0.5*(erf(d1/sqrt(2))+1)-(S2+E)*0.5*(erf(d2/sqrt(2))+1));
    //N() function has 1/sqrt(2pi), so y is different from erf()
    //more specific, 0.5*(erf(d1/sqrt(2))+1)=N(d1)
    return y;
}

double kirk_function2(double S1, double S2, double sigma1, double sigma2){
    double y, d1, d2, sigmak;
    sigmak=sqrt(sigma1*sigma1-2*rho*sigma1*sigma2*S2/(S2+E*exp(-r*T))+(S2/(S2+E*exp(-r*T)))*(S2/(S2+E*exp(-r*T)))*sigma2*sigma2);
    d1= (log(S1/(S2+E*exp(-r*T)))+0.5*sigmak*sigmak*T)/(sigmak*sqrt(T));
    d2= d1-sigmak*sqrt(T);
    y= S1*0.5*(erf(d1/sqrt(2))+1)-(S2+E*exp(-r*T))*0.5*(erf(d2/sqrt(2))+1);
    //N() function has 1/sqrt(2pi), so y is different from erf()
    //more specific, 0.5*(erf(d1/sqrt(2))+1)=N(d1)
    return y;
}

