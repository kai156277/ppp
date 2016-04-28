#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H


class math_function
{
public:
    math_function();
    static double lagrange(double fx, double* form, int size);          //Lagrange Interpolation Polynomial
    static void sunPostion(double juliday, double *posCTS);
    static void moonPostion(double juliday,double *posCTS);

private:
    static void getSunPosCIS(double juliday, double *posCIS);
    static void getMoonPosCIS(double juliday);
    static void CIS2CTS(const spaceRectangular &posCIS,double juliday);

};

#endif // MATH_FUNCTION_H
