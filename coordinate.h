#ifndef COORDINATE_H
#define COORDINATE_H


class Coordinate
{
public:
    Coordinate();
    ~Coordinate();
    Coordinate &setENU(double sx,double sy,double sz,double rx,double ry,double rz);
    Coordinate &setxyn(double xf,double yf,double nf,double j);
    Coordinate &setXYZ(double Xf,double Yf,double Zf);
    Coordinate &setBLH(double Bf,double Lf,double Hf);
    Coordinate &setB(double d, double m, double s);
    Coordinate &setL(double d, double m, double s);
    Coordinate &setH(double h);
    Coordinate &ENUparameter();
    Coordinate &XYZtoBLH();
    Coordinate &BLHtoXYZ();
    Coordinate &xyntoBLH();
    Coordinate &BLHtoxyn();
    Coordinate &XYZtoxyn();
    Coordinate &xyntoXYZ();
    double J;               //j高斯投影的分度带
    double x,y,n;           //高斯投影坐标及带号
    double X,Y,Z;           //大地直角坐标系
    double B,L,H;           //大地经纬坐标系
    double Sx, Sy, Sz, Rx, Ry, Rz;
    double azimuth, elevation, distance;
private:
    static double a;        //椭球长半轴
    static double E;        //椭球第一偏心率的平方
    static double PI;       //圆周率
    static double A1;
    static double B1;
    static double C1;
    static double D1;
    static double a1;
    static double b1;
    static double c1;
    static double d1;
    static double b2;
    static double c2;
    static double d2;       //A1-d2子午线弧长计算公式所需常量
    static double K0;
    static double K1;
    static double K2;
    static double K3;       //A1-K3高斯坐标反算所需常量
    Coordinate &RotationalMatrices(double Senu[]);
};


#endif // COORDINATE_H
