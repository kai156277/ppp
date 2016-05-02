#ifndef PPP_MODEL_H
#define PPP_MODEL_H

#include "ppp_data.h"
#include <QVector>
#include <Eigen/Eigen>

using namespace Eigen;

class ppp_model
{
public:
    ppp_model();
    void generic_model(const ppp_file &file, result_file &result);
private:
    void kalman_prepare();
    void building_kalman(MatrixXd &Xt, MatrixXd &F, MatrixXd &Pt, MatrixXd &Q,
                        MatrixXd &H,  MatrixXd &R, MatrixXd &Z,
                        const Vector4d &X1 , const Matrix4d &Q1,
                        const MatrixXd &B1 , const ppp_epoch &epoch);
    QVector<bool> cycle_slip(const ppp_epoch &epoch1,const ppp_epoch &epoch2);
    void point_positioning(const ppp_epoch &epoch, Vector4d &X, Matrix4d &Q, MatrixXd &B);

    static double receiver_x;
    static double receiver_y;
    static double receiver_z;

    const static double Pi;

    const static double c;
    /*GPS载波频率(Hz)------------------------------------------------------------*/
    const static double f1;
    const static double f2;
    const static double f5;
    /*GPS载波的波长--------------------------------------------------------------*/
    const static double lambda1;
    const static double lambda2;
};

#endif // PPP_MODEL_H
