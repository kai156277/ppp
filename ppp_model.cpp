#include "ppp_model.h"
#include "math_function.h"
#include <Eigen/Eigen>
#include <cmath>

using namespace Eigen;

const double ppp_model::Pi = 3.141592653;
const double ppp_model::c = 299792458;
/*GPS载波频率--------------------------------------------------------------------*/
const double ppp_model::f1 = 154 * 10.23 * 1e6;
const double ppp_model::f2 = 120 * 10.23 * 1e6;
const double ppp_model::f5 = 115 * 10.23 * 1e6;
/*GPS载波波长--------------------------------------------------------------------*/
const double ppp_model::lambda1 = ppp_model::c / ppp_model::f1;
const double ppp_model::lambda2 = ppp_model::c / ppp_model::f2;
double ppp_model::receiver_x = 0;
double ppp_model::receiver_y = 0;
double ppp_model::receiver_z = 0;

ppp_model::ppp_model()
{

}

void ppp_model::generic_model(const ppp_file &file, result_file &result)
{
    receiver_x = file.receiver_x;
    receiver_y = file.receiver_y;
    receiver_z = file.receiver_z;
    /*-------------------------------计算初始值----------------------------------*/
    Vector4d X1;
    Matrix4d Q1;
    MatrixXd B1,Xt, F, Pt, Q, H,  R, Z;
    point_positioning(file.file[0],X1,Q1,B1);
    Matrix4d I;
    I.setIdentity(4,4);

    building_kalman(Xt,F,Pt,Q,H,R,Z,X1,Q1,B1,file.file[0]);

    for(int i = 0; i<400; i++)
    {
        /*1.周跳探测-------------------------------------------------------------*/
        QVector<bool> cycle = cycle_slip(file.file[i],file.file[i+1]);

        /*2.kalman filter*/
        math_function::kalman_filter(Xt,F,Pt,Q,H,R,Z,I);
        ppp_X X;
        X.dx = Xt(0);   X.dy = Xt(1);   X.dz = Xt(2);   X.dt = Xt(3);
        result.file.push_back(X);
    }
}

void ppp_model::building_kalman(MatrixXd &Xt, MatrixXd &F, MatrixXd &Pt, MatrixXd &Q,
                               MatrixXd &H, MatrixXd &R, MatrixXd &Z,
                               const Vector4d &X1, const Matrix4d &Q1,const MatrixXd &B1,
                               const ppp_epoch & epoch)
{
    int sate_num = epoch.sate_info.size();
    int num = 5 + sate_num;
    /*1.状态向量-----------------------------------------------------------------*/
    VectorXd X(num);
    for(int i = 0; i<4; i++)
    {
        X << X1(0);
    }
    X << 0;    //天顶对流层延迟
    for(int i = 0; i<sate_num; i++)
    {
        X << 0;
    }
    Xt = X;


    /*2.状态转移矩阵-------------------------------------------------------------*/
    F.setIdentity(num,num);

    /*3.状态噪声的协方差矩阵------------------------------------------------------*/
    Pt.setZero(num,num);
    double sigma2 = 900.0 / Q1(0,0) + Q1(1,1) + Q1(2,2);
    for(int i = 0; i<4; i++)
    {
        Pt(i,i) = sigma2 * Q(i,i);
    }

    /*4. 系统状态噪声的协方差阵---------------------------------------------------*/
    Q.setZero(num,num);
    Q(3,3) = 900;
    Q(4,4) = 1e-8;

    /*5.观测矩阵----------------------------------------------------------------*/
    H = B1;

    /*6.观测噪声的协方差矩阵*/
    R.setZero(sate_num,sate_num);
    double r_sigma2 = 4e-6;
    for(int i = 0; i<sate_num; i++)
    {
        double radian = epoch.sate_info[i].elevation / 180 *Pi;
        if(epoch.sate_info[i].elevation > 30)
        {
            R(i,i) = r_sigma2 / sin(radian);
        }
        else
        {
            R(i,i) = r_sigma2 / (sin(radian)*sin(radian));
        }
    }

    /*7.实际观测值*/
    VectorXd v(sate_num);
    for(int i = 0; i<sate_num; i++)
    {
        v(i) = epoch.sate_info[i].windup + epoch.sate_info[i].offsetL1 + epoch.sate_info[i].offsetL2 +
               epoch.sate_info[i].trop_delay + epoch.sate_info[i].tide_effect +
               epoch.sate_info[i].stat_ant_height + epoch.sate_info[i].sate_antenna +
               epoch.sate_info[i].sagnac + epoch.sate_info[i].relativity ;
    }
    Z = H * Xt + v;

}

QVector<bool> ppp_model::cycle_slip(const ppp_epoch &epoch1, const ppp_epoch &epoch2)
{
    QVector <bool> cycle(epoch2.sate_info.size(),true);
    bool list = epoch1.PRN == epoch2.PRN;
    if(list == true)
    {
        for ( int i =  0; i<epoch1.sate_info.size(); i++)
        {
            double L11 = epoch1.sate_info[i].L1;        //历元一卫星的L1
            double L12 = epoch2.sate_info[i].L1;        //历元二卫星的L1;
            double L21 = epoch1.sate_info[i].L2;        //历元一卫星的L2
            double L22 = epoch2.sate_info[i].L2;        //历元二卫星的L2;
            double P11 = epoch1.sate_info[i].P1;        //历元一卫星的P1
            double P12 = epoch2.sate_info[i].P1;        //历元二卫星的P1;
            double P21 = epoch1.sate_info[i].P2;        //历元一卫星的P2
            double P22 = epoch2.sate_info[i].P2;        //历元二卫星的P2;
            double slip1 = (L12 - L11) - ( f1/f2 )*(L22 - L21);       //电离层残差法
            double L1 = (lambda1*L11) - (lambda2*L21) + (P11-P21);        //双频P码
            double L2 = (lambda1*L12) - (lambda2*L22) + (P12-P22);
            double slip2 = L2 - L1;
            double Nw1 = (L11-L21) - ((f1-f2)/(f1+f2))*((P11/lambda1)+(P21/lambda2));//MW组合法
            double Nw2 = (L12-L22) - ((f1-f2)/(f1+f2))*((P12/lambda1)+(P22/lambda2));
            double slip3 = Nw2 -Nw1;
            if((fabs(slip1)>0.3) || (fabs(slip2)>3.5) || (fabs(slip3)>3.5))
            {
                cycle[i] = false;
            }
        }
        return cycle;
    }
    else
    {
        for( int i = 0; i<epoch2.sate_info.size(); i++)
        {
            QString PRN = epoch2.PRN[i];
            int index = epoch1.PRN.indexOf(PRN);
            double L11 = epoch1.sate_info[index].L1;        //历元一卫星的L1
            double L12 = epoch2.sate_info[i].L1;            //历元二卫星的L1;
            double L21 = epoch1.sate_info[index].L2;        //历元一卫星的L2
            double L22 = epoch2.sate_info[i].L2;            //历元二卫星的L2;
            double P11 = epoch1.sate_info[index].P1;        //历元一卫星的P1
            double P12 = epoch2.sate_info[i].P1;            //历元二卫星的P1;
            double P21 = epoch1.sate_info[index].P2;        //历元一卫星的P2
            double P22 = epoch2.sate_info[i].P2;            //历元二卫星的P2;
            double slip1 = (L12 - L11) - ( f1/f2 )*(L22 - L21);       //电离层残差法
            double L1 = (lambda1*L11) - (lambda2*L21) + (P11-P21);        //双频P码
            double L2 = (lambda1*L12) - (lambda2*L22) + (P12-P22);
            double slip2 = L2 - L1;
            double Nw1 = (L11-L21) - ((f1-f2)/(f1+f2))*((P11/lambda1)+(P21/lambda2));//MW组合法
            double Nw2 = (L12-L22) - ((f1-f2)/(f1+f2))*((P12/lambda1)+(P22/lambda2));
            double slip3 = Nw2 -Nw1;
            if((fabs(slip1)>0.3) || (fabs(slip2)>3.5) || (fabs(slip3)>3.5))
            {
                cycle[i] = false;
            }
        }
    }
}

void ppp_model::point_positioning(const ppp_epoch &epoch, Vector4d &X, Matrix4d &Q, MatrixXd &B)
{
    int num = epoch.sate_info.size();
    VectorXd r(num);          //近似伪距
    VectorXd ionized(num);    //电离层
    VectorXd t(num);          //钟差
    VectorXd w(num);          //其他误差
    VectorXd l(num);
    MatrixXd A(num,4);        //系数阵


    for(int j = 0; j<num ;j++) //r
    {
        /*1.近似伪距-------------------------------------------------------------*/
        r(j) = epoch.sate_info[j].distance;

        /*2.系数阵---------------------------------------------------------------*/
        A(j,0) = ( epoch.sate_info[j].position_x - receiver_x ) / r(j);
        A(j,1) = ( epoch.sate_info[j].position_y - receiver_y ) / r(j);
        A(j,2) = ( epoch.sate_info[j].position_z - receiver_z ) / r(j);
        A(j,3) = -1;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(j) = ( epoch.sate_info[j].P1*f1*f1 - epoch.sate_info[j].P2*f2*f2 ) / ( f1*f1 - f2*f2 );

        /*4.钟差----------------------------------------------------------------*/
        t(j) = epoch.sate_info[j].clock;

        /*5.其他误差-------------------------------------------------------------*/
        w(j) = epoch.sate_info[j].relativity + epoch.sate_info[j].sagnac + epoch.sate_info[j].trop_delay
             + epoch.sate_info[j].tide_effect + epoch.sate_info[j].stat_ant_height + epoch.sate_info[j].sate_antenna
             + epoch.sate_info[j].offsetL1 + epoch.sate_info[j].offsetL2 + epoch.sate_info[j].windup;
    }

    l = r - ionized - t + w;

    /*-----------------------------计算参数的改正数-------------------------------*/
    X = (A.transpose() * A).inverse() * A.transpose() * l;

    /*-------------------------------Q 协方差阵----------------------------------*/
    Q = A.transpose() * A;

    /*--------------------------------B系数阵------------------------------------*/
    B = A;
}

