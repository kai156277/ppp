#include "coordinate.h"
#include <math.h>

double Coordinate::a = 6378137;
double Coordinate::E = 0.00669437999014132;
double Coordinate::PI = 3.141592653;
double Coordinate::A1 = 1 + 3.0/4*E + 45.0/64*pow(E,2) + 175.0/256*pow(E,3) + 11025.0/16384*pow(E,4) + 43659.0/65536*pow(E,5);
double Coordinate::B1 =     3.0/4*E + 15.0/16*pow(E,2) + 525.0/512*pow(E,3) +   2205.0/2048*pow(E,4) + 72765.0/65536*pow(E,5);
double Coordinate::C1 =               15.0/64*pow(E,2) + 105.0/256*pow(E,3) +   2205.0/4096*pow(E,4) + 10395.0/16384*pow(E,5);
double Coordinate::D1 =                                   35.0/512*pow(E,3) +    315.0/2048*pow(E,4) + 31185.0/13072*pow(E,5);
double Coordinate::a1 = A1*a*(1-E);
double Coordinate::b1 = -1.0 / 2*B1*a*(1-E);
double Coordinate::c1 =  1.0 / 4*C1*a*(1-E);
double Coordinate::d1 = -1.0 / 6*D1*a*(1-E);
double Coordinate::b2 = -b1 / a1;
double Coordinate::c2 = -c1 / a1;
double Coordinate::d2 = -d1 / a1;
double Coordinate::K0 = a1;
double Coordinate::K1 = 2*b1 + 4*c1 + 6*d1;
double Coordinate::K2 = -8 * c1 - 32 * d1;
double Coordinate::K3 = 32 * d1;
Coordinate::Coordinate()
    :B(0),L(0),H(0),
     X(0),Y(0),Z(0),
     x(0),y(0),n(0),J(6),
     Sx(0),Sy(0),Sz(0),
     Rx(0),Ry(0),Rz(0),
     azimuth(0), elevation(0), distance(0)
{

}

Coordinate::~Coordinate()
{

}

Coordinate &Coordinate::setENU(double sx, double sy, double sz, double rx, double ry, double rz)
{
    Sx = sx;
    Sy = sy;
    Sz = sz;
    Rx = rx;
    Ry = ry;
    Rz = rz;
    return *this;
}


Coordinate &Coordinate::setxyn(double xf, double yf, double nf,double j)
{
    x = xf;
    y = yf;
    n = nf;
    J = j;
    return *this;
}

Coordinate &Coordinate::setXYZ(double Xf, double Yf, double Zf)
{
    X = Xf;
    Y = Yf;
    Z = Zf;
    return *this;
}

Coordinate &Coordinate::setBLH(double Bf, double Lf, double Hf)
{
    B = Bf;
    L = Lf;
    H = Hf;
    return *this;
}

Coordinate &Coordinate::setB(double d, double m, double s)
{
    B=((d*3600+m*60+s)/3600)*PI/180;
    return *this;
}

Coordinate &Coordinate::setL(double d, double m, double s)
{
    L=((d*3600+m*60+s)/3600)*PI/180;
    return *this;
}

Coordinate &Coordinate::setH(double h)
{
    H = h;
    return *this;
}

Coordinate &Coordinate::XYZtoBLH()
{
    double N;
    double i;
    B=atan(Z/(sqrt(X*X+Y*Y)));
    for(i=0;i<=20;i++)
    {
        N=a/(sqrt(1-E*sin(B)*sin(B)));
        B=atan((Z+N*E*sin(B))/sqrt(X*X+Y*Y));
    }
    L=atan2(Y,X);
    H=(sqrt(X*X+Y*Y)/cos(B))-N;
    return *this;

}

Coordinate &Coordinate::XYZtoxyn()
{
    XYZtoBLH();
    BLHtoxyn();
    return *this;
}

Coordinate &Coordinate::BLHtoXYZ()
{
    double N;
    N=a/(sqrt(1-E*sin(B)*sin(B)));
    X=(N+H)*cos(B)*cos(L);
    Y=(N+H)*cos(B)*sin(L);
    Z=(N*(1-E)+H)*sin(B);
    return *this;

}

Coordinate &Coordinate::BLHtoxyn()
{

    int Ld;
    double L0;
    double l,t,m,NN,N,x0;

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
    return *this;

}

Coordinate &Coordinate::xyntoBLH()
{
    int i;
    double bi[4];
    double ci[4];
    double di[4];
    double B0,Bf,NNf,Vf,tf,Nf,l,k1,k2,k3;
    double L0;
    if(J==6)
        L0=6*n-3;
    else
        L0=3*n;
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
    B=Bf-0.5*Vf*Vf*tf*(y/Nf)*(y/Nf)+(1.0/24)*(5+3*tf*tf+NNf-9*NNf*tf*tf)*Vf*Vf*tf*pow((y/Nf),4)-(1.0/720)*(61+90*tf*tf+45*pow(tf,4))*Vf*Vf*tf*pow((y/Nf),6);
    l=(1/cos(Bf))*(y/Nf)-(1.0/6)*(1+2*tf*tf+NNf)*(1/cos(Bf))*pow((y/Nf),3)+(1.0/120)*(5+28*tf*tf+24*pow(tf,4)+6*NNf+8*NNf*tf*tf)*(1/cos(Bf))*pow((y/Nf),5);
    l=l/PI*180;
    L=(L0+l)/180*PI;
    return *this;
}

Coordinate &Coordinate::xyntoXYZ()
{
    xyntoBLH();
    BLHtoXYZ();
    return *this;
}


Coordinate &Coordinate::ENUparameter()
{
    double Senu[3]={0,0,0}; //Senu[0] = x,Senu[1] = y,Senu[2]=z
    setXYZ(Rx,Ry,Rz).XYZtoBLH().RotationalMatrices( Senu );
    distance = sqrt( pow(Senu[0],2) + pow(Senu[1],2) + pow(Senu[2],2));
    elevation = atan( Senu[2] / (sqrt( pow(Senu[0],2) + pow(Senu[1],2) )) ) * 180 / PI;
    azimuth = atan( Senu[0] / Senu[1] ) *180 /PI;
    if(Senu[1] < 0)
    {
        azimuth += 180;
    }
    else if ((Senu[0] < 0) && (Senu[1]) >= 0)
    {
        azimuth += 360;
    }
    return *this;
}


Coordinate &Coordinate::RotationalMatrices(double Senu[])
{
    void Matrices31Multiplication(const double a[][3],const double b[],double c[]);
    double M[3][3] = {0};
    double xyz[3] = {0};
    double Sxyz[3] = {Sx,Sy,Sz};
    double Rxyz[3] = {Rx,Ry,Rz};
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
        xyz[i] = Sxyz[i] - Rxyz[i];
    }
    Matrices31Multiplication(M,xyz,Senu);
    return *this;
}

void Matrices31Multiplication(const double M[][3], const double xyz[], double Senu[])
{
    for(int i = 0; i<3; i++ )
    {
        for(int j = 0; j<3; j++)
        {
            Senu[i] += M[i][j] * xyz[j];
        }
    }
}
