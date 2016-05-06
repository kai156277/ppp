#include <iostream>
#include "math_function.h"


math_function::math_function()
{

}

/* -----------------------------------------------------------------------------
 * Lagrange插值函数
 * PARAMETER   I/O   TYPE       DESCRIPTION
 *    fx       IN    double     插值求解点
 *    *form    IN    *double    查找表（前一半为提供的x，后一半为y）
 *    size     IN    int        查找表的长度
 *    return   OUT   double     插值结果
 * -----------------------------------------------------------------------------*/
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

/* -----------------------------------------------------------------------------
 * Xt：t历元下的状态向量
 * F ：状态转移矩阵
 * Pt：状态噪声的协方差矩阵
 * Q ：动态噪声的协方差矩阵
 * H ：系数阵
 * R ：观测噪声的协方差阵
 * Z ：实际观测值
 * I ：单位阵
 * -----------------------------------------------------------------------------*/

void math_function::kalman_filter(      MatrixXd &Xt, const MatrixXd &F,
                                        MatrixXd &Pt, const MatrixXd &Q,
                                  const MatrixXd &H,  const MatrixXd &R,
                                  const MatrixXd &Z,  const MatrixXd &I)
{
    MatrixXd HT = H.transpose();
    MatrixXd FT = F.transpose();

    /*更新时间：预测--------------------------------------------------------------*/
    MatrixXd X = F * Xt;
    MatrixXd P = F * Pt* FT + Q;

    /*计算增益矩阵K--------------------------------------------------------------*/
    MatrixXd Kt_1= (H * P * HT + R).inverse();
    MatrixXd Kt  = P * HT * Kt_1;

    /*更新观测值：改正------------------------------------------------------------*/
    Xt = X + Kt * (Z - H * X);
    Pt = (I - Kt * H) * P;

}
