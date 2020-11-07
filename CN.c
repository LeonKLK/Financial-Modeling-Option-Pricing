//
//  main.c
//  CN_method
//
//  Created by Leon Kwok on 2/10/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double ini_con(double, double);
double Nmin(double);
double Nmax(double);
double OP(double, double, double);
double SP(double, double);
double exact(double);       //function declaration
double exact_delta_function(double);            //analytic delta
double exact_gamma_function(double);
double exact_theta_function(double, double);

double a, b, k;
double sigma, r, E, T, M, alpha;
double x, dx, dt;
double xmin, xmax, sd;

int main() {
    
    clock_t begin, end;
    begin = clock();
    FILE *option_price;
    FILE *exact_option_price;
    FILE *payoff;
    FILE *delta;
    FILE *exact_delta;
    FILE *gamma;
    FILE *exact_gamma;
    
    option_price=fopen("C_price.csv", "w");
    exact_option_price=fopen("exact_sol.csv", "w");
    payoff=fopen("payoff.csv", "w");
    delta=fopen("delta.csv", "w");
    exact_delta=fopen("exact_delta.csv", "w");
    gamma=fopen("gamma.csv","w");
    exact_gamma=fopen("exact_gamma.csv", "w");
    
    int N;
    
    r= 0.05;
    sigma= 0.2;
    E= 10;
    T= 0.5; //6 months but T is in the unit of years
    
    sd= sigma*sqrt(T);
    xmin=-10*sd;
    xmax=10*sd;
    dx= 0.01;
    dt= 0.00005;
    
    alpha= dt/(dx*dx);
    N= (xmax-xmin)/dx;
    M= 0.5*sigma*sigma*T/dt;
    
    double b_m[N+1], y[N+1], q[N+1], u[N+1], price_c[N+1], SP_x[N+1], exact_c[N+1], error_c[N+1];;
    double payoff_x[N+1];
    double delta_c[N+1], exact_delta_c[N+1], error_delta[N+1];
    double gamma_c[N+1], exact_gamma_c[N+1], error_gamma[N+1];
    
    k=r/(0.5*sigma*sigma);
    a= -0.5*(k-1);
    b= -0.25*(k+1)*(k+1);
    
    x=xmin;
    for (int i=0; i<=N; ++i) {      //putting initial condition
        u[i]= ini_con(x, x);
        x=x+dx;
    }
    
    //setting z[i]
    for (int i=1; i<=N-1 ; ++i) {
        b_m[i]= (1-alpha)*u[i]+alpha/2*(u[i-1]+u[i+1]);
    }
    b_m[1]=b_m[1]+alpha/2*Nmin(xmin);
    b_m[N-1]=b_m[N-1]+alpha/2*Nmax(xmax);
    
    //setting y[i]
    y[1]= 1+2*alpha;
    for (int i=2; i<=N-1; ++i) {
        y[i]=1+alpha-(alpha*alpha/4)/y[i-1];
    }
    
    //setting q[i]
    q[1]=b_m[1];
    for (int i=2; i<=N-1; ++i) {
        q[i]=b_m[i]+alpha/2*q[i-1]/y[i-1];
    }
    
    //setting u[i]
    u[N-1]= q[N-1]/y[N-1];
    for (int i=N-2 ; i>=1; --i) {
        u[i]= (q[i]+alpha/2*u[i+1])/y[i];
    }
    u[0]= Nmin(xmin);
    u[N]= Nmax(xmax);
    //initialization is done.
    
    for (int j=0; j<M-1 ; ++j) {
        x= xmin;
        for (int i=1; i<=N-1; ++i) {
            x=x+dx;
            b_m[i]= (1-alpha)*u[i]+alpha/2*(u[i-1]+u[i+1]);
            //printf("x=%lf %lf\n", x, b_m[i]);
        }
        b_m[1]=b_m[1]+alpha/2*Nmin(xmin);
        b_m[N-1]=b_m[N-1]+alpha/2*Nmax(xmax);
        //b[i] setting done
        
        //setting q[i]
        q[1]=b_m[1];
        for (int i=2; i<=N-1; ++i) {
            q[i]=b_m[i]+alpha/2*q[i-1]/y[i-1];
        }
        
        u[N-1]= q[N-1]/y[N-1];
        for (int i=N-2 ; i>=1; --i) {
            u[i]= (q[i]+alpha/2*u[i+1])/y[i];
        }
    }
    
    x=xmin;
    for (int i=0; i<=N ; ++i) {     //calculate option price c from u
        SP_x[i] =SP(x, E);
        price_c[i]=OP(x, M*dt, u[i]);
        exact_c[i]=exact(SP_x[i]);
        error_c[i]=fabs(price_c[i]-exact_c[i])/exact_c[i]*100;
        if (x>=0) {
            payoff_x[i]=SP_x[i]-E;
        } else {
            payoff_x[i]=0;
        }
        printf("x=%lf S=%lf C=%lf exact_c=%lf error=%lf\n", x, SP_x[i], price_c[i], exact_c[i], error_c[i]);
        x= x+dx;
    }
    
    //The following code of "fprint" corresponds for printing data file.
    x=xmin;
    for (int i=0; i<=N; ++i) {
        fprintf(option_price,"%lf\t %lf\t %lf\n", SP_x[i], price_c[i], error_c[i]);   //We extract the option price within 3 sd of x
        fprintf(exact_option_price,"%lf\t %lf\n", SP_x[i], exact_c[i]);
        //printf("x=%lf S=%lf C=%lf exact_c=%lf error=%lf\n", x,SP_x[i], price_c[i], exact_c[i], error_c[i]);
        exact_delta_c[i]= exact_delta_function(SP_x[i]);        //we only calculate delta within 3 sd of x
        delta_c[i]=(price_c[i+1]-price_c[i-1])/(SP_x[i+1]-SP_x[i-1]);     //central difference with second order accuracy
        error_delta[i]=-(exact_delta_c[i]-delta_c[i])/exact_delta_c[i]*100;
        //printf("S=%lf delta=%lf exact_delta=%lf error=%lf\n ",SP_x[i], delta_c[i], exact_delta_c[i], error_delta[i]);
        fprintf(delta,"%lf\t %lf\t %lf\n", SP_x[i], delta_c[i], error_delta[i]);
        fprintf(exact_delta, "%lf\t %lf\n", SP_x[i], exact_delta_c[i]);
    x= x+dx;
    }
    
    x=xmin;
    for (int i=0; i<=N; ++i) {
        gamma_c[i]=(delta_c[i+1]-delta_c[i-1])/(SP_x[i+1]-SP_x[i-1]);
        exact_gamma_c[i]=exact_gamma_function(SP_x[i]);
        error_gamma[i]=(gamma_c[i]-exact_gamma_c[i])/exact_gamma_c[i]*100;
        //printf("C=%lf S=%lf gamma=%lf exact_gamma=%lf error=%lf\n ",price_c[i], SP_x[i], gamma_c[i], exact_gamma_c[i], error_gamma[i]);
        fprintf(gamma,"%lf\t %lf\t %lf\n", SP_x[i], gamma_c[i], error_gamma[i]);
        fprintf(exact_gamma,"%lf\t %lf\n", SP_x[i], exact_gamma_c[i]);
    x= x+dx;
    }
    
    fclose(option_price);
    fclose(exact_option_price);
    fclose(payoff);
    fclose(delta);
    fclose(exact_delta);
    fclose(gamma);
    fclose(exact_gamma);
    
    end = clock();
    printf ("time used to evaluate program: %f microsecs\n", (double)(end- begin)/(double)CLOCKS_PER_SEC);
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

double exact(double S){
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
