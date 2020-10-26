//
// Created by leanne on 10/25/20.
//

#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cassert>
#include <numeric>
using namespace std;


// Returns factorial of n
static int fact(int n)
{
    int result = 1;

    for (int i=1; i<=n; ++i) {
      result *= i;
    }

    return result;
}


// Function definition
static int nCr(int n, int r)
{
    return fact(n) / (fact(r) * fact(n - r));
}

static double my_fct(double t, double a, double k){
    return (t>= a ? pow(t-a, k) : 0);
}


template <typename A>
static double ar(int r, double t, const A &a)
{
    return a(r, t);
}

double volterra::rect(double t, double low, double up){
    return ((t>=low)  && (t< up)? 1 : 0);
}


//double volterra::impul(double t){
//    return (t==0 ? 1 : 0);
//}
double
volterra::beta(
        int n,
        int r,
        std::function<double(int r, double t)> a,
        double t
)
{
    double value = 0;
    int l;

    if (n=0, r=0){
        return 1;
    }
    else if (n=0, r>=1){
        return 0;
    }
    else if (n=1, r>=0){
        return ar(r,t,a);
    }
    else {
        for (l = 0; l <= n; l++)
        {
            value += beta(n,r,a,t) + ar(l,t,a)*beta(n-l,r,a,t);
        }

        // lacking an argument, also need to incorporate kernel
        return value;
    }
}


void volterra::norms(double* x, int n, double& norm1, double& norminf)
{
    norm1 = fabs(x[0]);
    norminf = norm1;
    for(int i=1; i<n; i++)
    {
        x[i] = fabs(x[i]);
        norm1 += x[i];
        if(x[i]> norminf)
            norminf = x[i];
    }
}


double volterra::findError(double h1,double h2)
{
    return fabs(h1-h2);
}


double volterra::series(int n, double t){
// calculating the value of (n-1)!

    int n1fact = fact(n-1);
    double sum = 0.0;

// loop to display the series
    for ( int i=0; i<n+1; i++) {

        int mPow = pow(-1,i);
        int comb = nCr(n, i);
        double myfct = my_fct(t,i,n-1);
        sum += mPow * comb * myfct/ n1fact;
    }
    return sum;
}
