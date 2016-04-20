#include "math_function.h"

math_function::math_function()
{

}

double math_function::lagrange(double fx, double *form, int size)
{
    size = size / 2;
    double y = 0;
    for(int i = 0; i<size; i++)
    {
        double x = 1;
        for(int j = 0; j<size; j++)
        {
            if(i == j)
            {
                continue;
            }
            x *= (fx - form[j]) / (form[i] - form[j]);
        }
        x = x * form[size+i];
        y += x;
    }
    return y;
}

