#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H


class math_function
{
public:
    math_function();
    static double lagrange(double fx, double* form, int size);          //Lagrange Interpolation Polynomial

};

#endif // MATH_FUNCTION_H
