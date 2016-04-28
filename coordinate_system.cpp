#include <cmath>

#include "coordinate_system.h"
#include "iostream"

using namespace std;

const double coordinate_transform::a = 6378137;
const double coordinate_transform::E = 0.00669437999014132;
const double coordinate_transform::PI = 3.141592653;
const double coordinate_transform::A1 = 1 + 3.0/4*E + 45.0/64*pow(E,2) + 175.0/256*pow(E,3) + 11025.0/16384*pow(E,4) + 43659.0/65536*pow(E,5);
const double coordinate_transform::B1 =     3.0/4*E + 15.0/16*pow(E,2) + 525.0/512*pow(E,3) +   2205.0/2048*pow(E,4) + 72765.0/65536*pow(E,5);
const double coordinate_transform::C1 =               15.0/64*pow(E,2) + 105.0/256*pow(E,3) +   2205.0/4096*pow(E,4) + 10395.0/16384*pow(E,5);
const double coordinate_transform::D1 =                                   35.0/512*pow(E,3) +    315.0/2048*pow(E,4) + 31185.0/13072*pow(E,5);
const double coordinate_transform::a1 = A1*a*(1-E);
const double coordinate_transform::b1 = -1.0 / 2*B1*a*(1-E);
const double coordinate_transform::c1 =  1.0 / 4*C1*a*(1-E);
const double coordinate_transform::d1 = -1.0 / 6*D1*a*(1-E);
const double coordinate_transform::b2 = -b1 / a1;
const double coordinate_transform::c2 = -c1 / a1;
const double coordinate_transform::d2 = -d1 / a1;
const double coordinate_transform::K0 = a1;
const double coordinate_transform::K1 = 2*b1 + 4*c1 + 6*d1;
const double coordinate_transform::K2 = -8 * c1 - 32 * d1;
const double coordinate_transform::K3 = 32 * d1;
const double BLH_coordinate::Pi = 3.141592653;
const double ENU_coordiante::Pi = 3.141592653;


coordinate_transform::coordinate_transform()
{

}

coordinate_transform::~coordinate_transform()
{

}

XYZ_coordinate coordinate_transform::BLHtoXYZ(const BLH_coordinate &BLH)
{
    double N;
    double B = BLH.get_latitude_radians();
    double L = BLH.get_longitude_radians();
    double H = BLH.getHeight();
    N=a/(sqrt(1-E*sin(B)*sin(B)));
    double X=(N+H)*cos(B)*cos(L);
    double Y=(N+H)*cos(B)*sin(L);
    double Z=(N*(1-E)+H)*sin(B);
    XYZ_coordinate XYZ(X,Y,Z);
    return XYZ;
}

XYZ_coordinate coordinate_transform::xyntoXYZ(const xyn_coordiante &xyn)
{
    BLH_coordinate BLH = coordinate_transform::xyntoBLH(xyn);
    return coordinate_transform::BLHtoXYZ(BLH);
}

BLH_coordinate coordinate_transform::XYZtoBLH(const XYZ_coordinate &XYZ)
{
    double X = XYZ.getX();
    double Y = XYZ.getY();
    double Z = XYZ.getZ();
    double N;
    double i;
    double B=atan(Z/(sqrt(X*X+Y*Y)));
    for(i=0;i<=20;i++)
    {
        N=a/(sqrt(1-E*sin(B)*sin(B)));
        B=atan((Z+N*E*sin(B))/sqrt(X*X+Y*Y));
    }
    double L=atan2(Y,X);
    double H=(sqrt(X*X+Y*Y)/cos(B))-N;
    BLH_coordinate BLH;
    BLH.setBLH_radians(B,L,H);
    return BLH;
}

BLH_coordinate coordinate_transform::xyntoBLH(const xyn_coordiante &xyn)
{
    double x = xyn.getX();
    double y = xyn.getY();
    int i;
    double bi[4];
    double ci[4];
    double di[4];
    double B0,Bf,NNf,Vf,tf,Nf,l,k1,k2,k3;
    double L0;
    if(xyn.getJ()==6)
        L0=6*xyn.getN()-3;
    else
        L0=3*xyn.getN();
    bi[0]=b2;
    ci[0]=c2;
    di[0]=d2;
    for(i=1;i<=3;i++)
    {
        bi[i]=bi[0]+bi[0]*ci[i-1]-1.5*bi[0]*bi[i-1]*bi[i-1]-2*ci[0]*bi[i-1];
        ci[i]=ci[0]+bi[0]*bi[i-1];
        di[i]=di[0]+bi[0]*ci[i-1]+0.5*bi[0]*bi[i-1]*bi[i-1]+2*ci[0]*bi[i-1];
    }
    k1=2*bi[3]+4*ci[3]+6*di[3];
    k2=8*ci[3]+32*di[3];
    k3=32*di[3];
    B0=x/a1;
    Bf=B0+cos(B0)*(k1*sin(B0)-k2*pow(sin(B0),3)+k3*pow(sin(B0),5));
    NNf=(E/(1-E))*pow(cos(Bf),2);
    Vf=sqrt(1+NNf);
    tf=tan(Bf);
    Nf=a/sqrt(1-E*pow(sin(Bf),2));
    double B=Bf-0.5*Vf*Vf*tf*(y/Nf)*(y/Nf)+(1.0/24)*(5+3*tf*tf+NNf-9*NNf*tf*tf)*Vf*Vf*tf*pow((y/Nf),4)-(1.0/720)*(61+90*tf*tf+45*pow(tf,4))*Vf*Vf*tf*pow((y/Nf),6);
    l=(1/cos(Bf))*(y/Nf)-(1.0/6)*(1+2*tf*tf+NNf)*(1/cos(Bf))*pow((y/Nf),3)+(1.0/120)*(5+28*tf*tf+24*pow(tf,4)+6*NNf+8*NNf*tf*tf)*(1/cos(Bf))*pow((y/Nf),5);
    l=l/PI*180;
    double L=(L0+l)/180*PI;
    BLH_coordinate BLH;
    BLH.setBLH_radians(B,L,0);
    return BLH;
}

xyn_coordiante coordinate_transform::BLHtoxyn(const BLH_coordinate &BLH)
{
    xyn_coordiante xyn;
    int Ld;
    double L0;
    double l,t,m,NN,N,x0;
    double B = BLH.get_latitude_radians();
    double L = BLH.get_longitude_radians();
    double x,y,n;
    Ld = L*180/PI/3600;
    if(Ld%6==0)
        n=Ld/6+2;
    else
        n=Ld/6+1;
    n=20;
    L0=(6*n-3)*PI/180;
    l=L-L0;
    t=tan(B);
    m=l*cos(B);
    NN=E/(1-E)*cos(B)*cos(B);
    N=a/sqrt(1-E*sin(B)*sin(B));
    x0=K0*B+cos(B)*(K1*sin(B)+K2*pow(sin(B),3)+K3*pow(sin(B),5));
    x=x0+0.5*N*t*m*m+1.0/24*(5-t*t+9*NN+4*NN*NN)*N*t*pow(m,4)+1.0/720*(61-58*t*t+pow(t,4))*N*t*pow(m,6);
    y=N*m+1.0/6*(1-t*t+NN)*N*pow(m,3)+1.0/120*(5-18*t*t+pow(t,4)+14*NN-58*NN*t*t)*N*pow(m,6);
    xyn.set_xyn(x,y,n,6);
    return xyn;
}

xyn_coordiante coordinate_transform::XYZtoxyn(const XYZ_coordinate &XYZ)
{
    BLH_coordinate BLH = coordinate_transform::XYZtoBLH(XYZ);
    return coordinate_transform::BLHtoxyn(BLH);
}

ENU_coordiante coordinate_transform::sateENU(const XYZ_coordinate &station_XYZ, const XYZ_coordinate &satellite_XYZ)
{
    BLH_coordinate station_BLH = XYZtoBLH(station_XYZ);
    double B = station_BLH.get_latitude_radians();
    double L = station_BLH.get_longitude_radians();

    double M[3][3] = {0};
    double xyz[3] = {0};
    double satellite_xyz[3] = {satellite_XYZ.getX(),
                               satellite_XYZ.getY(),
                               satellite_XYZ.getZ()};
    double station_xyz[3] = {station_XYZ.getX(),
                             station_XYZ.getY(),
                             station_XYZ.getZ()};
    double sate_ead[3] = {0};
    /*XYZ到ENU的旋转矩阵*/
    M[0][0] = -sin( L );
    M[0][1] = cos( L );
    M[0][2] = 0;
    M[1][0] = -sin( B ) * cos( L );
    M[1][1] = -sin( B ) * sin( L );
    M[1][2] = cos( B );
    M[2][0] = cos( B ) * cos( L );
    M[2][1] = cos( B ) * sin( L );
    M[2][2] = sin( B );

    for(int i =0; i<3 ;i++)
    {
        xyz[i] = satellite_xyz[i] - station_xyz[i];
    }
    for(int i = 0; i<3; i++ )
    {
        for(int j = 0; j<3; j++)
        {
            sate_ead[i] += M[i][j] * xyz[j];
        }
    }

    ENU_coordiante sate_ENU;
    sate_ENU.setENU_meter(sate_ead[0],sate_ead[1],sate_ead[2]);
    return sate_ENU;
}



BLH_coordinate::BLH_coordinate()
    :latitude(0),longitude(0),height(0)
{

}

BLH_coordinate::~BLH_coordinate()
{

}

void BLH_coordinate::setBLH_radians(double B, double L, double H)
{
    height = H;
    if(B > Pi/2.0 || B < -Pi/2.0)
    {
        cerr << "You want set a error latitude(radins)!" << endl
             << "That value is " << B << endl;
        exit(1);
    }
    else
    {
        latitude = B * 180 / Pi;
    }

    if(L > Pi || L < 0)
    {
        cerr << "Your want set a error longitude(radians)!" << endl
             << "That value is " << L << endl;
        exit(1);
    }
    else
    {
        longitude = L * 180 / Pi;
    }
}

void BLH_coordinate::setBLH_angle(double B, double L, double H)
{
    height = H;
    if(B < -90 || B > 90)
    {
        cerr << "You want set a error latitude(angle)!" << endl
             << "That value is " << B << endl;
        exit(1);
    }
    else
    {
        latitude = B;
    }

    if(L < 0 || L > 360)
    {
        cerr << "You want set a error longitude(angle)!" << endl
             << "Tant value is " << L << endl;
        exit(1);
    }
    else
    {
        longitude = L;
    }
}

double BLH_coordinate::get_latitude_radians() const
{
    return latitude / 180 * Pi;
}

double BLH_coordinate::get_longitude_radians() const
{
    return longitude / 180 * Pi ;
}

double BLH_coordinate::getHeight() const
{
    return height;
}

double BLH_coordinate::get_latitude_angle() const
{
    return latitude;
}

double BLH_coordinate::get_longitude_angle() const
{
    return longitude;
}

XYZ_coordinate::XYZ_coordinate()
    :X(0),Y(0),Z(0)
{

}

XYZ_coordinate::XYZ_coordinate(double x, double y, double z)
{
    setXYZ(x,y,z);
}

XYZ_coordinate::~XYZ_coordinate()
{

}

void XYZ_coordinate::setXYZ(double x, double y, double z)
{
    X = x;
    Y = y;
    Z = z;
}

double XYZ_coordinate::getX() const
{
    return X;
}

double XYZ_coordinate::getY() const
{
    return Y;
}

double XYZ_coordinate::getZ() const
{
    return Z;
}



ENU_coordiante::ENU_coordiante()
    :sate_e(0),sate_n(0),sate_u(0),
      azimuth(0),elevation(0),distance(0)
{

}

ENU_coordiante::~ENU_coordiante()
{

}

void ENU_coordiante::setENU_meter(double e, double n, double u)
{
    sate_e = e;
    sate_n = n;
    sate_u = u;

    distance = sqrt( pow(sate_e,2) + pow(sate_n,2) + pow(sate_u,2));
    elevation = atan( sate_u / (sqrt( pow(sate_e,2) + pow(sate_n,2) )) ) * 180 / Pi;
    azimuth = atan( sate_e / sate_n ) *180 / Pi;
    if(sate_n < 0)
    {
        azimuth += 180;
    }
    else if ((sate_e < 0) && (sate_n) >= 0)
    {
        azimuth += 360;
    }
}

void ENU_coordiante::setENU_angle(double e, double a, double d)
{
    if(e > 90)
    {
        cerr << "You want set a error elevation!" << endl
             << "That value is " << e << endl;
        exit(1);
    }
    else
    {
        elevation = e;
    }

    if(a < 0)
    {
        cerr << "You want set a error azimuth!" << endl
             << "That value is " << a << endl;
        exit(1);
    }
    else
    {
        azimuth = a;
    }
    if(d < 0)
    {
        cerr << "You want set a error distance!" << endl
             << "That value is " << d << endl;
    }
    else
    {
        distance = d;
    }
    sate_u = distance * sin(elevation * Pi / 180);
    double line = distance * cos(elevation * Pi / 180);
    sate_e = line * sin(azimuth * Pi / 180);
    sate_n = line * cos(azimuth * Pi / 180);
}

double ENU_coordiante::getSate_e() const
{
    return sate_e;
}
double ENU_coordiante::getSate_n() const
{
    return sate_n;
}
double ENU_coordiante::getSate_u() const
{
    return sate_u;
}
double ENU_coordiante::getElevation() const
{
    return elevation;
}
double ENU_coordiante::getAzimuth() const
{
    return azimuth;
}
double ENU_coordiante::getDistance() const
{
    return distance;
}

xyn_coordiante::xyn_coordiante()
{

}

xyn_coordiante::~xyn_coordiante()
{

}

void xyn_coordiante::set_xyn(double x, double y, double n, double J)
{
    this->x = x;
    this->y = y;
    this->n = n;
    if( J == 3 )
    {
        this->J = 3;
    }
    else if( J == 6 )
    {
        this->J = 6;
    }
    else
    {
        cerr << "You want set a cerr J in Gauss projection coordinates!" << endl
             << "That value is " << J << endl;
        exit(1);
    }
}

double xyn_coordiante::getX() const
{
    return x;
}
double xyn_coordiante::getY() const
{
    return y;
}
double xyn_coordiante::getN() const
{
    return n;
}
double xyn_coordiante::getJ() const
{
    return J;
}



