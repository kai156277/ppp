#ifndef COORDINATE_SYSTEM_H
#define COORDINATE_SYSTEM_H


class BLH_coordinate
{
public:
    BLH_coordinate();
    ~BLH_coordinate();

    void setBLH_radians(double B, double L, double H);
    void setBLH_angle(double B, double L, double H);

    double get_latitude_radians() const;
    double get_latitude_angle() const;
    double get_longitude_radians() const;
    double get_longitude_angle() const;
    double getHeight() const;
private:
    /* -------------------------------------------------------------------------
     * latitude and longitude is angle,
     * north latitude is positive number,
     * south latitude is negative number,
     * east longitude is 0~180,
     * west longitude is 180~360,
     * height is geodetic height.
     * -------------------------------------------------------------------------*/
    double latitude;
    double longitude;
    double height;

    const static double Pi;
};

class XYZ_coordinate
{
public:
    XYZ_coordinate();
    XYZ_coordinate(double x ,double y ,double z);
    ~XYZ_coordinate();
    void setXYZ(double x, double y, double z);
    double getX() const;
    double getY() const;
    double getZ() const;

private:
    /* -------------------------------------------------------------------------
     * earth-rectangular coordinate,
     * X,Y,Z is meter.
     * -------------------------------------------------------------------------*/
    double X;
    double Y;
    double Z;
};

class ENU_coordiante
{
public:
    ENU_coordiante();
    ~ENU_coordiante();
    void setENU_meter(double e,double n,double u);
    void setENU_angle(double e,double a,double d);
    double getSate_e() const;
    double getSate_n() const;
    double getSate_u() const;
    double getElevation() const;
    double getAzimuth() const;
    double getDistance() const;

private:
    /* -------------------------------------------------------------------------
     * sate_e,sate_n,sate_u is meter,
     * elevation,azimuth is angle,
     * distance is meter.
     * -------------------------------------------------------------------------*/
    double sate_e;
    double sate_n;
    double sate_u;
    double elevation;
    double azimuth;
    double distance;
    static const double Pi;
};

class xyn_coordiante
{
public:
    xyn_coordiante();
    ~xyn_coordiante();
    void set_xyn(double x,double y,double n,double J);
    double getX() const;
    double getY() const;
    double getN() const;
    double getJ() const;

private:
    /* -------------------------------------------------------------------------
     * Gauss projection coordinates,
     * x,y is meter,
     * n is reel number,
     * J is 6 angle or 3 angle.
     * -------------------------------------------------------------------------*/
    double x;
    double y;
    double n;
    double J;
};

class coordinate_transform
{
public:
    coordinate_transform();
    ~coordinate_transform();
    static XYZ_coordinate BLHtoXYZ(const BLH_coordinate &BLH);
    static XYZ_coordinate xyntoXYZ(const xyn_coordiante &xyn);
    static BLH_coordinate XYZtoBLH(const XYZ_coordinate &XYZ);
    static BLH_coordinate xyntoBLH(const xyn_coordiante &xyn);
    static xyn_coordiante BLHtoxyn(const BLH_coordinate &BLH);
    static xyn_coordiante XYZtoxyn(const XYZ_coordinate &XYZ);
    static ENU_coordiante sateENU(const XYZ_coordinate &station_XYZ,const XYZ_coordinate &satellite_XYZ);
private:
    static const double a;        //椭球长半轴
    static const double E;        //椭球第一偏心率的平方
    static const double PI;       //圆周率
    static const double A1;
    static const double B1;
    static const double C1;
    static const double D1;
    static const double a1;
    static const double b1;
    static const double c1;
    static const double d1;
    static const double b2;
    static const double c2;
    static const double d2;       //A1-d2子午线弧长计算公式所需常量
    static const double K0;
    static const double K1;
    static const double K2;
    static const double K3;       //A1-K3高斯坐标反算所需常量
};

#endif // COORDINATE_SYSTEM_H
