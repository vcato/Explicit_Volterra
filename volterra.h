//
// Created by leanne on 10/25/20.
//

#ifndef EXPLICIT_VOLTERRA_VOLTERRA_H
#define EXPLICIT_VOLTERRA_VOLTERRA_H

#include <string>
#include <functional>


namespace volterra {
    double impul(double t);
    double rect(double t, double low, double up);
    double series(int n, double t);

    double findError(double h1,double h2);

    double
    beta(
            int n,
            int r,
            std::function<double(int r, double t)> a,
            double t
    );

    void norms(double* x, int n, double& norm1, double& norminf);

}

#endif //EXPLICIT_VOLTERRA_VOLTERRA_H
