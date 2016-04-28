#include <iostream>
#include "math_function.h"

math_function::math_function()
{

}

/*
 * Lagrange插值函数
 * PARAMETER   I/O   TYPE       DESCRIPTION
 *    fx       IN    double     插值求解点
 *    *form    IN    *double    查找表（前一半为提供的x，后一半为y）
 *    size     IN    int        查找表的长度
 *    return   OUT   double     插值结果
*/
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
