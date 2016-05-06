#include "ppp_model.h"
#include "math_function.h"
#include <Eigen/Eigen>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

const double ppp_model::Pi = 3.141592653;
const double ppp_model::c = 299792458;
/*GPS载波频率--------------------------------------------------------------------*/
const double ppp_model::f1 = 154 * 10.23 * 1e6;
const double ppp_model::f2 = 120 * 10.23 * 1e6;
const double ppp_model::f5 = 115 * 10.23 * 1e6;
/*GPS载波波长--------------------------------------------------------------------*/
const double ppp_model::lambda = c / (f1 + f2);
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
    fstream f4("C:\\PPP\\readFile\\1\\kalman.txt",ios::app);

    receiver_x = file.receiver_x;
    receiver_y = file.receiver_y;
    receiver_z = file.receiver_z;
    int size = file.file.size() - 1;
    /*-------------------------------计算初始值----------------------------------*/
    VectorXd X1;
    MatrixXd B1,Xt, F, Pt, Q, H,R, Z,I;
    point_positioning(file.file[0],X1,B1,Z);
    building_kalman(Xt,F,Pt,Q,H,R,Z,I,X1,B1,Z,file.file[0]);
    math_function::kalman_filter(Xt,F,Pt,Q,H,R,Z,I);
    bool flag = true;
    for(int i = 0; i<size; i++)
    {
        /*1.周跳探测-------------------------------------------------------------*/
        QVector<bool> cycle = cycle_slip(file.file[i],file.file[i+1]);
        if((cycle.indexOf(true) >= 0)||cycle.size()!=file.file[i].sate_number)
        {
            /*存在周跳-----------------------------------------------------------*/
            building_kalman(cycle,file.file[i+1],file.file[i].PRN,Xt,H,Pt,R,Z,F,Q,I);
        }
        else
        {
            /*不存在周跳---------------------------------------------------------*/
            building_kalman(file.file[i+1],H,Z,R);
            flag = false;
        }

        /*2.kalman filter*/
        math_function::kalman_filter(Xt,F,Pt,Q,H,R,Z,I);
        f4 << file.file[i].year << "-" << file.file[i].month << "-" << file.file[i].day << "-"
           << file.file[i].hour << "-" << file.file[i].minute << "-" <<file.file[i].second
           << setw(15) << Xt.transpose() << endl;


        std::cout << i << "\t";
    }
    f4.close();
}

void ppp_model::building_kalman(MatrixXd &Xt, MatrixXd &F, MatrixXd &Pt, MatrixXd &Q,
                               MatrixXd &H, MatrixXd &R, MatrixXd &Z, MatrixXd &I,
                               const VectorXd &X1, const MatrixXd &B1,
                               const MatrixXd &L1, const ppp_epoch & epoch)
{
    int sate_num = epoch.sate_info.size();
    int num = 5 + sate_num;
    /*1.状态向量-----------------------------------------------------------------*/
    Xt = X1;

    /*2.状态转移矩阵-------------------------------------------------------------*/
    F.setIdentity(num,num);

    /*3.状态向量的协方差矩阵------------------------------------------------------*/
    //Pt = Nbb_ ;

    Pt.setZero(num,num);
    for(int i = 0; i<3; i++)
    {
        Pt(i,i) = 1e5;
    }
    Pt(3,3) = 1e6;
    Pt(4,4) = 900;
    for(int i = 0; i<sate_num; i++)
    {
        Pt(5+i,5+i) = 1e20;
    }


    /*4.动态噪声的协方差阵--------------------------------------------------------*/
    Q.setZero(num,num);
    Q(3,3) = 900;
    Q(4,4) = 1e-8;

    /*5.观测矩阵-----------------------------------------------------------------*/
    H = B1;

    /*6.观测值的协方差矩阵--------------------------------------------------------*/
    R.setZero(sate_num*2,sate_num*2);
    for(int i = 0; i<sate_num; i++)
    {
        double radian = epoch.sate_info[i].elevation / 180 *Pi;
        if(epoch.sate_info[i].elevation > 30)
        {
            R(i,i) = 0.009 / sin(radian);
            R(i+sate_num,i+sate_num) = 4e-6 / sin(radian);
        }
        else
        {
            R(i,i) = 0.009 / (sin(radian)*sin(radian));
            R(i+sate_num,i+sate_num) = 4e-6 / (sin(radian)*sin(radian));
        }
    }

    /*7.实际观测值*/
    Z = L1;

    /*8.单位阵*/
    I.setIdentity(num,num);

}

void ppp_model::building_kalman(const ppp_epoch &epoch, MatrixXd &B, MatrixXd &L, MatrixXd &R)
{
    int num = epoch.sate_info.size();
    int dnum = num * 2;
    VectorXd r(dnum);          //近似伪距
    VectorXd ionized(dnum);    //电离层
    VectorXd w(dnum);          //其他误差
    VectorXd l(dnum);
    B.setZero(dnum,5+num);

    /*1.构造B矩阵----------------------------------------------------------------*/
    /*-----------------------------构造伪距的观测方程-----------------------------*/
    for(int i = 0; i<num ;i++)
    {
        /*1.近似伪距-------------------------------------------------------------*/
        r(i) = epoch.sate_info[i].distance;

        /*2.系数阵---------------------------------------------------------------*/
        B(i,0) = ( epoch.sate_info[i].position_x - receiver_x ) / r(i);
        B(i,1) = ( epoch.sate_info[i].position_y - receiver_y ) / r(i);
        B(i,2) = ( epoch.sate_info[i].position_z - receiver_z ) / r(i);
        B(i,3) = -1;
        B(i,4) = epoch.sate_info[i].trop_map;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(i) = epoch.sate_info[i].ionized_pseudo;


        /*4.其他误差-------------------------------------------------------------*/
        w(i) = epoch.sate_info[i].clock + epoch.sate_info[i].relativity + epoch.sate_info[i].sagnac
             - epoch.sate_info[i].trop_delay + epoch.sate_info[i].tide_effect +
               epoch.sate_info[i].stat_ant_height + epoch.sate_info[i].sate_antenna;
    }

    /*-----------------------------构造载波的观测方程-----------------------------*/
    for(int i = num; i<dnum ;i++)
    {
        int j = i - num;
        /*1.近似伪距-------------------------------------------------------------*/
        r(i) = r(j);

        /*2.系数阵---------------------------------------------------------------*/
        B(i,0) = B(j,0);
        B(i,1) = B(j,1);
        B(i,2) = B(j,2);
        B(i,3) = -1;
        B(i,4) = epoch.sate_info[j].trop_map;
        B(i,5+i-num) = lambda;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(i) = epoch.sate_info[j].ionized_carrier;


        /*4.其他误差-------------------------------------------------------------*/
        w(i) = w(j);
    }
    /*2.构造L矩阵----------------------------------------------------------------*/
    l = r - ionized  - w;
    L = l;

    /*4.观测值的协方差矩阵--------------------------------------------------------*/
    R.setZero(dnum,dnum);
    for(int i = 0; i<num; i++)
    {
        double radian = epoch.sate_info[i].elevation / 180 *Pi;
        if(epoch.sate_info[i].elevation > 30)
        {
            R(i,i) = 0.009 / sin(radian);
            R(i+num,i+num) = 4e-6 / sin(radian);
        }
        else
        {
            R(i,i) = 0.009 / (sin(radian)*sin(radian));
            R(i+num,i+num) = 4e-6 / (sin(radian)*sin(radian));
        }
    }

}

void ppp_model::building_kalman(const QVector<bool> &cycle, const ppp_epoch &epoch,
                                const QStringList &epoch1,
                                MatrixXd &X, MatrixXd &B, MatrixXd &Pt, MatrixXd &R,
                                MatrixXd &L, MatrixXd &F, MatrixXd &Q , MatrixXd &I)
{
    /*--------------------------------------------------------------------------
     * epoch1 当前历元的卫星信息
     * cycle 下一历元的卫星周跳信息
     * epoch 下一历元的卫星信息
     * -------------------------------------------------------------------------*/
    int num = epoch.sate_info.size();
    int dnum = num * 2;
    VectorXd r(dnum);          //近似伪距
    VectorXd ionized(dnum);    //电离层
    VectorXd w(dnum);          //其他误差
    VectorXd l(dnum);
    B.setZero(dnum,5+num);

    /*1.构造B矩阵----------------------------------------------------------------*/
    /*-----------------------------构造伪距的观测方程-----------------------------*/
    for(int i = 0; i<num ;i++)
    {
        /*1.近似伪距-------------------------------------------------------------*/
        r(i) = epoch.sate_info[i].distance;

        /*2.系数阵---------------------------------------------------------------*/
        B(i,0) = ( epoch.sate_info[i].position_x - receiver_x ) / r(i);
        B(i,1) = ( epoch.sate_info[i].position_y - receiver_y ) / r(i);
        B(i,2) = ( epoch.sate_info[i].position_z - receiver_z ) / r(i);
        B(i,3) = -1;
        B(i,4) = epoch.sate_info[i].trop_map;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(i) = epoch.sate_info[i].ionized_pseudo;


        /*4.其他误差-------------------------------------------------------------*/
        w(i) = epoch.sate_info[i].clock + epoch.sate_info[i].relativity + epoch.sate_info[i].sagnac
             - epoch.sate_info[i].trop_delay + epoch.sate_info[i].tide_effect +
               epoch.sate_info[i].stat_ant_height + epoch.sate_info[i].sate_antenna;
    }

    /*-----------------------------构造载波的观测方程-----------------------------*/
    for(int i = num; i<dnum ;i++)
    {
        int j = i - num;
        /*1.近似伪距-------------------------------------------------------------*/
        r(i) = r(j);

        /*2.系数阵---------------------------------------------------------------*/
        B(i,0) = B(j,0);
        B(i,1) = B(j,1);
        B(i,2) = B(j,2);
        B(i,3) = -1;
        B(i,4) = epoch.sate_info[j].trop_map;
        B(i,5+i-num) = lambda;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(i) = epoch.sate_info[j].ionized_carrier;


        /*4.其他误差-------------------------------------------------------------*/
        w(i) = w(j);
    }

    /*2.构造L矩阵----------------------------------------------------------------*/
    l = r - ionized  - w;
    L = l;

    /*3.状态向量的协方差矩阵------------------------------------------------------*/
    QVector<double> Pt0(num,0);
    QVector<double> Pt1(5,0);
    for(int i = 0; i<cycle.size();i++)
    {
        if(cycle[i] == true)
        {
            Pt0[i] = 1e20;
        }
        else
        {
            int index = epoch1.indexOf(epoch.PRN[i]);
            Pt0[i] = Pt(5+index,5+index);

        }
    }
    for(int i = 0; i<5; i++)
    {
        Pt1[i] = Pt(i,i);
    }
    Pt.setZero(5+num,5+num);
    for(int i = 0; i<5; i++)
    {
        Pt(i,i) = Pt1[i];
    }
    for(int i = 0; i<num; i++)
    {
        Pt(5+i,5+i) = Pt0[i];
    }

    /*4.观测值的协方差矩阵--------------------------------------------------------*/
    R.setZero(dnum,dnum);
    for(int i = 0; i<num; i++)
    {
        double radian = epoch.sate_info[i].elevation / 180 *Pi;
        if(epoch.sate_info[i].elevation > 30)
        {
            R(i,i) = 0.009 / sin(radian);
            R(i+num,i+num) = 4e-6 / sin(radian);
        }
        else
        {
            R(i,i) = 0.009 / (sin(radian)*sin(radian));
            R(i+num,i+num) = 4e-6 / (sin(radian)*sin(radian));
        }
    }

    /*5.状态向量-----------------------------------------------------------------*/
    MatrixXd X2;
    X2 = (B.transpose() * B).inverse() * B.transpose() * l;
    QVector<double> X0(num,0);
    QVector<double> X1(5,0);
    for(int i = 0; i<cycle.size();i++)
    {
        if(cycle[i] == true)
        {
            X0[i] = X2(5+i,0);
        }
        else
        {
            int index = epoch1.indexOf(epoch.PRN[i]);
            X0[i] = X(5+index,0);

        }
    }
    for(int i = 0 ; i<5; i++)
    {
        X1[i] = X(i,0);
    }
    X.setZero(5+num,1);
    for(int i = 0; i<5; i++)
    {
        X(i,0) = X1[i];
    }
    for(int i = 0 ; i<cycle.size(); i++)
    {
        X(5+i,0) = X0[i];
    }



    /*6.状态转移矩阵-------------------------------------------------------------*/
    F.setIdentity(5+num,5+num);

    /*7.动态噪声的协方差阵--------------------------------------------------------*/
    Q.setZero(5+num,5+num);
    Q(3,3) = 900;
    Q(4,4) = 1e-8;

    /*8.单位阵------------------------------------------------------------------*/
    I.setIdentity(5+num,5+num);
}

QVector<bool> ppp_model::cycle_slip(const ppp_epoch &epoch1, const ppp_epoch &epoch2)
{
    /*--------------------------------------------------------------------------
     * false,代表没有发生周跳
     * true,代表发生周跳
     * -------------------------------------------------------------------------*/
    QVector <bool> cycle(epoch2.sate_info.size(),false);
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
                cycle[i] = true;
            }
        }
        return cycle;
    }
    else if(list!=true )
    {
        for( int i = 0; i<epoch2.sate_info.size(); i++)
        {
            QString PRN = epoch2.PRN[i];
            int index = epoch1.PRN.indexOf(PRN);
            if(index == -1)
            {
                cycle[i] = true;
            }
            else
            {
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
                    cycle[i] = true;
                }
            }
        }
    }
}

void ppp_model::point_positioning(const ppp_epoch &epoch, VectorXd &X,MatrixXd &B, MatrixXd &L)
{
    /*---------------------------伪距/载波联合单点定位-----------------------------
     * 1.将两个定位方程联合用于一次解出 dx,dy,dz,dt,T和每颗卫星的N
     * 2.误差方程为：V = A * X - L
     * 3.V阵的结构为：
     *      [Vp,Vc]T
     * 4.A阵的结构为：
     *      B   T   0    此行为伪距的系数
     *      B   T   C    此行为载波的系数
     *      B为dx,dy,dz,dt的系数阵
     *      T为天顶对流层的系数阵
     *      C为整周模糊度的系数阵
     * 5.X阵的结构为：
     *      [dx dy dz dt Tztd N1 N2 ...Nn]T
     * -------------------------------------------------------------------------*/
    int num = epoch.sate_info.size();
    int dnum = num * 2;

    VectorXd r(dnum);          //近似伪距
    VectorXd ionized(dnum);    //电离层
    VectorXd w(dnum);          //其他误差
    VectorXd l(dnum);
    MatrixXd Nbb_;
    B.setZero(dnum,5+num);

    /*-----------------------------构造伪距的观测方程-----------------------------*/
    for(int i = 0; i<num ;i++)
    {
        /*1.近似伪距-------------------------------------------------------------*/
        r(i) = epoch.sate_info[i].distance;

        /*2.系数阵---------------------------------------------------------------*/
        B(i,0) = ( epoch.sate_info[i].position_x - receiver_x ) / r(i);
        B(i,1) = ( epoch.sate_info[i].position_y - receiver_y ) / r(i);
        B(i,2) = ( epoch.sate_info[i].position_z - receiver_z ) / r(i);
        B(i,3) = -1;
        B(i,4) = epoch.sate_info[i].trop_map;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(i) = epoch.sate_info[i].ionized_pseudo;


        /*4.其他误差-------------------------------------------------------------*/
        w(i) = epoch.sate_info[i].clock + epoch.sate_info[i].relativity + epoch.sate_info[i].sagnac
             - epoch.sate_info[i].trop_delay + epoch.sate_info[i].tide_effect +
               epoch.sate_info[i].stat_ant_height + epoch.sate_info[i].sate_antenna;
    }

    /*-----------------------------构造载波的观测方程-----------------------------*/
    for(int i = num; i<dnum ;i++)
    {
        int j = i - num;
        /*1.近似伪距-------------------------------------------------------------*/
        r(i) = r(j);

        /*2.系数阵---------------------------------------------------------------*/
        B(i,0) = B(j,0);
        B(i,1) = B(j,1);
        B(i,2) = B(j,2);
        B(i,3) = -1;
        B(i,4) = epoch.sate_info[j].trop_map;
        B(i,5+i-num) = lambda;

        /*3.电离层延迟-----------------------------------------------------------*/
        ionized(i) = epoch.sate_info[j].ionized_carrier;


        /*4.其他误差-------------------------------------------------------------*/
        w(i) = w(j);
    }
    l = r - ionized  - w;

    /*--------------------------计算参数的方差-协方差阵----------------------------*/
    Nbb_ = (B.transpose() * B).inverse();

    /*-----------------------------计算参数的改正数-------------------------------*/
    X = Nbb_ * B.transpose() * l;

    L = l;

}

void ppp_model::row_swap(MatrixXd &a, int a_row, const MatrixXd &b, int b_row,int col)
{
    for(int i  = 0; i<col; i++)
    {
        a(a_row,i) = b(b_row,i);
    }
}

