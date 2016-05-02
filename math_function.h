#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include<Eigen/Eigen>

using namespace Eigen;

class math_function
{
public:
    math_function();
    static double lagrange(double fx, double* form, int size);          //Lagrange Interpolation Polynomial
    static void kalman_filter(MatrixXd &Xt, const MatrixXd &F,
                              MatrixXd &Pt, const MatrixXd &Q,
                        const MatrixXd &H,  const MatrixXd &R,
                        const MatrixXd &Z,  const MatrixXd &I);
};

#endif // MATH_FUNCTION_H
