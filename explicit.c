//
//  main.c
//  explicit finite difference method
//
//  Created by Leon Kwok on 23/7/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double ini_con(double, double);     //I.C.
double Nmin(double);    //B.C.
double Nmax(double);    //B.C.
double OP(double, double, double);      //calculate option price from u
double SP(double, double);              //stock price
double exact_price_function(double);       //exact solution of option price
double exact_delta_function(double);       //analytic delta
double exact_gamma_function(double);
double exact_theta_function(double, double);

double sigma, r, E, T, M;
int N;
double a, b, k;     //a,b are for calculating real option price, and k is for N_max B.C.
double alpha;       //different globle variable, for iteraction calculation

int main(){
    double x, dx, dt, sd;
    double xmin, xmax;
    clock_t begin, end;
    
    begin = clock();
    FILE *option_price;
    FILE *exact_option_price;
    FILE *payoff;
    FILE *delta;
    FILE *exact_delta;
    FILE *gamma;
    FILE *exact_gamma;
    //FILE *theta;
    //FILE *exact_theta;
    
    option_price=fopen("C_price.csv", "w");
    exact_option_price=fopen("exact_sol.csv", "w");
    payoff=fopen("payoff.csv", "w");
    delta=fopen("delta.csv", "w");
    exact_delta=fopen("exact_delta.csv", "w");
    gamma=fopen("gamma.csv","w");
    exact_gamma=fopen("exact_gamma.csv", "w");
    //theta=fopen("theta.csv", "w");
    //exact_theta=fopen("exact_theta.csv", "w");
    
    r= 0.05;        //interest rate
    sigma= 0.2;     //volatility
    E= 10;          //exercise price
    T= 0.5;         //6 months but T is in the unit of years
    dx= 0.01;
    dt= 0.000051;
    alpha= dt/(dx*dx);
    sd= sigma*sqrt(T);
    x=-3;
    xmin=-5*sd;
    xmax=5*sd;
    k=r/(0.5*sigma*sigma);
    a= -0.5*(k-1);
    b= -0.25*(k+1)*(k+1);
    
    printf("Alpha=%lf\n", alpha);
    
    N=fabs(x)/dx*2;
    M= 0.5*sigma*sigma*T/dt;
    printf("N=%d\n ", N);
    printf("M=%lf\n ", M);
    double old_u[N+1], new_u[N+1];
    double inter_u[N+1];
    double price_c[N+1], SP_x[N+1], exact_c[N+1], error_c[N+1], payoff_x[N+1];
    double delta_c[N+1], exact_delta_c[N+1], error_delta[N+1];
    double gamma_c[N+1], exact_gamma_c[N+1], error_gamma[N+1];
    double theta_c[N+1], exact_theta_c[N+1], error_theta[N+1];
    
    for (int i=0; i<=N; ++i) {      //putting initial condition
        old_u[i]= ini_con(x, x);
        //printf("x=%lf %lf\n", x, old_u[i]);
        x=x+dx;
    }
    
    for (int j=0; j<(M-1); ++j) {
        for (int i=1; i<=N-1 ; ++i) {
            new_u[i]= alpha*old_u[i-1]+(1-2*alpha)*old_u[i]+alpha*old_u[i+1];
        }
        new_u[0]= Nmin(-3);
        new_u[N]= Nmax(3);
        for (int i=0; i<=N; ++i) {
            inter_u[i]=old_u[i];       //storing the value of u at t=T-dt
            old_u[i]=new_u[i];
        }
    }
    
    x=-3;
    for (int i=0; i<=N ; ++i) {
        SP_x[i] =SP(x, E);
        price_c[i]=OP(x, M*dt, old_u[i]);   //calculate option price c from u
        exact_c[i]=exact_price_function(SP_x[i]);
        error_c[i]=(price_c[i]-exact_c[i])/exact_c[i]*100;
        inter_u[i]=OP(x, (M-1)*dt, inter_u[i]);     //calculate option price c from u at t=T-dt
        if (x>=0) {
            payoff_x[i]=SP_x[i]-E;
        } else {
            payoff_x[i]=0;
        }
        if (SP_x[i]>E*exp(xmin) && SP_x[i]<E*exp(xmax)) {
        //printf("x=%lf S=%lf C=%lf exact_c=%lf error=%lf\n", x,SP_x[i], price_c[i], exact_c[i], error_c[i]);
        }
        x= x+dx;
    }
    //whole calculation is done here.
    
    x=-3;
    for (int i=0; i<=N; ++i) {
        if (SP_x[i]>E*exp(xmin) && SP_x[i]<E*exp(xmax)) {
            fprintf(option_price,"%lf\t %lf\t %lf\n", SP_x[i], price_c[i], error_c[i]);   //We extract the option price within 3 sd of x
            fprintf(exact_option_price,"%lf\t %lf\n", SP_x[i], exact_c[i]);
            fprintf(payoff,"%lf\t %lf\n", SP_x[i], payoff_x[i]);
            //printf("x=%lf S=%lf C=%lf exact_c=%lf error=%lf\n", x,SP_x[i], price_c[i], exact_c[i], error_c[i]);
            exact_delta_c[i]= exact_delta_function(SP_x[i]);        //we only calculate delta within 3 sd of x
            delta_c[i]=(price_c[i+1]-price_c[i-1])/(SP_x[i+1]-SP_x[i-1]);     //central difference with second order accuracy
            error_delta[i]=-(exact_delta_c[i]-delta_c[i])/exact_delta_c[i]*100;
            //printf("S=%lf delta=%lf exact_delta=%lf error=%lf\n ",SP_x[i], delta_c[i], exact_delta_c[i], error_delta[i]);
            fprintf(delta,"%lf\t %lf\t %lf\n", SP_x[i], delta_c[i], error_delta[i]);
            fprintf(exact_delta, "%lf\t %lf\n", SP_x[i], exact_delta_c[i]);
        }
    x= x+dx;
    }
    
    x=-3;
    for (int i=0; i<=N; ++i) {
        if (SP_x[i]>E*exp(xmin) && SP_x[i]<E*exp(xmax)) {
            //gamma_c[i]=(price_c[i+1]-2*price_c[i]+price_c[i-1])/pow(SP_x[i+1]-SP_x[i], 2);      //order 2
            gamma_c[i]=(delta_c[i+1]-delta_c[i-1])/(SP_x[i+1]-SP_x[i-1]);
            //gamma_c[i]=(-price_c[i+2]+16*price_c[i+1]-30*price_c[i]+16*price_c[i-1]-price_c[i-2])/(12*pow(SP_x[i+1]-SP_x[i], 2));       //order 4
            exact_gamma_c[i]=exact_gamma_function(SP_x[i]);
            error_gamma[i]=(gamma_c[i]-exact_gamma_c[i])/exact_gamma_c[i]*100;
            //printf("C=%lf S=%lf gamma=%lf exact_gamma=%lf error=%lf\n ",price_c[i], SP_x[i], gamma_c[i], exact_gamma_c[i], error_gamma[i]);
            fprintf(gamma,"%lf\t %lf\t %lf\n", SP_x[i], gamma_c[i], error_gamma[i]);
            fprintf(exact_gamma,"%lf\t %lf\n", SP_x[i], exact_gamma_c[i]);
        }
    x= x+dx;
    }
    
    /*x=-3;
    for (int i=0; i<=N; ++i) {
        if (SP_x[i]>E*exp(xmin) && SP_x[i]<E*exp(xmax)) {
            exact_theta_c[i]=exact_theta_function(SP_x[i], 0);
            theta_c[i]=-(price_c[i]-inter_u[i])/(dt/(0.5*sigma*sigma));
            //error_theta[i]=fabs((exact_theta_c[i]-theta_c[i])/exact_theta_c[i])*100;
            fprintf(exact_theta, "%lf\t %lf\n", SP_x[i], exact_theta_c[i]);
            fprintf(theta, "%lf\t %lf\n", SP_x[i], theta_c[i]);
            //printf("S=%lf c=%lf u_p=%lf theta=%lf exact_theta=%lf error=%lf\n ",SP_x[i],price_c[i], inter_u[i], theta_c[i], exact_theta_c[i], error_theta[i]);
        }
    x= x+dx;
    }*/
    
    fclose(option_price);
    fclose(exact_option_price);
    fclose(payoff);
    fclose(delta);
    fclose(exact_delta);
    fclose(gamma);
    fclose(exact_gamma);
    //fclose(theta);
    //fclose(exact_theta);
    end = clock();
    
    printf ("time used to evaluate program: %f secs\n", (double)(end- begin)/(double)CLOCKS_PER_SEC);
    return 0;
}

double ini_con(double x_1, double x_2)
{
    double y;
    if (x_1>=0 && x_2>=0) {
        y=exp(0.5*(k+1)*x_1)-exp(0.5*(k-1)*(x_2));
    } else {
        y=0;
    }
    return y;
}

double Nmin(double x_min)
{
    return 0;
}

double Nmax(double x)
{
    double y;
    y= exp(0.5*(k+1)*x);
    return y;
}

double OP(double x, double tau, double u){
    double y;
    y= E*exp(a*x+b*tau)*u;
    return y;
}

double SP(double x, double E){
    double y;
    y= E*exp(x);
    return y;
}

double exact_price_function(double S){
    double y, d1, d2;       //integrate from negative infinity to x
    d1= (log(S/E)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));       //log in c is natural log
    d2= (log(S/E)+(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    y= S*0.5*(erf(d1/sqrt(2))+1)-E*exp(-r*T)*0.5*(erf(d2/sqrt(2))+1);       //N() function has 1/sqrt(2pi), so y is different from erf()
    return y;
}

double exact_delta_function(double S){
    double y, d1;
    d1= (log(S/E)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    y= 0.5*(erf(d1/sqrt(2))+1);
    return y;
}

double exact_gamma_function(double S){
    double y, d1;
    d1= (log(S/E)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    y= exp(-d1*d1/2)/(sqrt(2*M_PI*T)*S*sigma);
    return y;
}

double exact_theta_function(double S, double t){
    double y, d1, d2;
    d1= (log(S/E)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));       //log in c is natural log
    d2= (log(S/E)+(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    y=-S*sigma*exp(-d1*d1/2)/(sqrt(8*M_PI*(T-t)))-r*E*exp(-r*(T-t))*0.5*(erf(d2/sqrt(2))+1);
    return y;
}
