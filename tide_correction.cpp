#include "tide_correction.h"

tide_correction::tide_correction()
{

}

double tide_correction::tideEffect(Vector3d tide, neuParameter neuPos)
{
    double elev  = neuPos.elevation;
    double cosel = cos(elev);
    double sinel = sin(elev);
    double cosaz = cos(neuPos.azimuth);
    double sinaz = sin(neuPos.azimuth);
    double tide1 = tide[0]*cosel*cosaz + tide[1]*cosel*sinaz + tide[2]*sinel;
    return tide1;
}

Vector3d tide_correction::solidTide(XYZ_coordinate receiver, XYZ_coordinate sunPos, XYZ_coordinate moonPos, GPS_Time t)
{
    const double H20 = 0.6078;
    const double L20 = 0.0847;
    const double H3  = 0.292;
    const double L3  = 0.015;
    const double GM  = 3.986004415e14;
    const double GMS = 1.3271250e20;
    const double GMM = 4.9027890e12;
    const double AE  = 6378.137e3;
    //SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR
    Vector3d XSTA,XSUN,XMON;
    XSTA[0] = receiver.X; XSTA[1] = receiver.Y; XSTA[2] = receiver.Z;
    XSUN[0] = sunPos.X;   XSUN[1] = sunPos.Y;   XSUN[2] = sunPos.Z;
    XMON[0] = moonPos.X;  XMON[1] = moonPos.Y;  XMON[2] = moonPos.Z;
    double SCS  = XSTA.dot(XSUN);
    double RSTA = XSTA.norm();
    double RSUN = XSUN.norm();
    double SCM  = XSTA.dot(XMON);
    double RMON = XMON.norm();
    double SCSUN = SCS / RSTA / RSUN;
    double SCMON = SCM / RSTA / RMON;

    //COMPUTATION OF NEW H2 AND L2
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
    double H2     = H20 - 0.0006*(1.0-1.5*COSPHI*COSPHI);
    double L2     = L20 + 0.0002*(1.0-1.5*COSPHI*COSPHI);
    //P2-TERM
    double P2SUN = 3.0*(H2/2.0-L2)*SCSUN*SCSUN - H2/2.0;
    double P2MON = 3.0*(H2/2.0-L2)*SCMON*SCMON - H2/2.0;
    //P3-TERM
    double P3SUN = 2.5*(H3-3.0*L3)*SCSUN*SCSUN*SCSUN + 1.5*(L3-H3)*SCSUN;
    double P3MON = 2.5*(H3-3.0*L3)*SCMON*SCMON*SCMON + 1.5*(L3-H3)*SCMON;
    //TERM IN DIRECTION OF SUN/MOON VECTOR
    double X2SUN = 3.0 * L2 * SCSUN;
    double X2MON = 3.0 * L2 * SCMON;
    double X3SUN = 1.5 * L3 * (5.0*SCSUN*SCSUN-1.0);
    double X3MON = 1.5 * L3 * (5.0*SCMON*SCMON-1.0);
    //FACTORS FOR SUN/MOON
    double temp1 = AE / RSUN;
    double temp2 = AE / RMON;
    double FAC2SUN = GMS / GM*AE*temp1*temp1*temp1;
    double FAC2MON = GMM / GM*AE*temp2*temp2*temp2;
    double FAC3SUN = FAC2SUN * temp1;
    double FAC3MON = FAC2MON * temp2;

    Vector3d DXTIDE;
    //TOTAL DISPLACEMENT
    for (int i=0;i<3;i++)
    {
        DXTIDE[i] = FAC2SUN*( X2SUN*XSUN[i]/RSUN  + P2SUN*XSTA[i]/RSTA )
                  + FAC2MON*( X2MON*XMON[i]/RMON + P2MON*XSTA[i]/RSTA )
                  + FAC3SUN*( X3SUN*XSUN[i]/RSUN  + P3SUN*XSTA[i]/RSTA )
                  + FAC3MON*( X3MON*XMON[i]/RMON + P3MON*XSTA[i]/RSTA );
    }
    //CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I AND L_2^(0)I )
    //FIRST, FOR THE DIURNAL BAND
    Vector3d XCOSTA;
    XCOSTA = ST1DIU(XSTA,XSUN,XMON,FAC2SUN,FAC2MON);
    DXTIDE += XCOSTA;
    //SECOND, FOR THE SEMI-DIURNAL BAND
    XCOSTA = ST1SEM(XSTA,XSUN,XMON,FAC2SUN,FAC2MON);
    DXTIDE += XCOSTA;
    //CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )
    XCOSTA = ST1L1(XSTA,XSUN,XMON,FAC2SUN,FAC2MON);
    DXTIDE += XCOSTA;
    //CONSIDER CORRECTIONS FOR STEP 2
    // CORRECTIONS FOR THE DIURNAL BAND:
    // FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES
    julianTime JJ;
    gregorianTime gc;
    timeTransform::GPSTToJulian(t,JJ);
    timeTransform::GPSTToGregorian(t,gc);
    double T = (JJ.day-2451545.)/36525.0;
    //AND THE HOUR IN THE DAY
    double FHR = gc.hour + gc.minute/60 + (gc.second)/3600;

    XCOSTA = STEP2DIU(XSTA,FHR,T);
    DXTIDE += XCOSTA;
    //CORRECTIONS FOR THE LONG-PERIOD BAND:
    XCOSTA = ST2LON(XSTA,T);
    DXTIDE += XCOSTA;

    XYZ_coordinate rec; rec.X = XSTA[0];   rec.Y = XSTA[1];   rec.Z = XSTA[2];
    neuCoordinate res;
    geodetic recGeodetic;
    coordinateTransform::XYZtoBLH(rec,recGeodetic);

    double B = recGeodetic.B *= (pi/180);
    double L = recGeodetic.L *= (pi/180);

    Matrix3d M;
    M << -sin(B)*cos(L),-sin(B)*sin(L),cos(B),
         -sin(L),cos(L),0,
          cos(B)*cos(L), cos(B)*sin(L),sin(B);
    DXTIDE = M * DXTIDE;

    return DXTIDE;
}

Vector3d tide_correction::oceanTide(const ocean &ocean, GPS_Time t)
{
    double DEG_TO_RAD=0.017453;
    double NUM_COMPONENTS = 3;
    double NUM_HARMONICS = 11;
    MatrixXd arguments = getArg(t);
    Vector3d oLoading;
    for (int i=0;i<NUM_COMPONENTS;i++)
    {
        double temp = 0.0;
        for (int k=0;k<NUM_HARMONICS;k++)
        {
            temp = temp + ocean.data[i][k]*cos(arguments(0,k)-ocean.data[i+3][k]*DEG_TO_RAD);
        }
        oLoading[i] = temp;
    }
    oLoading[1] = -oLoading[1];
    oLoading[2] = -oLoading[2];
    Vector3d otide;
    otide[0] = oLoading[2];
    otide[1] = oLoading[1];
    otide[2] = oLoading[0];
    return otide;
}

Vector3d tide_correction::poleTide(const erpDataRecord &erp, XYZ_coordinate receiver, GPS_Time t)
{
    double xdisp = erp.xp;
    double ydisp = erp.yp;
    geodetic blh;
    coordinateTransform::XYZtoBLH(receiver,blh);
    double B = blh.B / 180.0 * pi;
    double L = blh.L / 180.0 * pi;
    //double H = blh.H;
    gregorianTime obstime(2000,1,1,12,0,0);
    julianTime j2000;
    julianTime t1;
    timeTransform::gregorianToJulian(obstime,j2000);
    timeTransform::GPSTToJulian(t,t1);
    double timedif = (t1.day - j2000.day) / 365.25;
    double xpbar = 0.054 + timedif*0.00083;
    double ypbar = 0.357 + timedif*0.00395;
    double m1 = xdisp-xpbar;
    double m2 = ypbar-ydisp;

    double sin2lat = sin(2.0*B);
    double cos2lat = cos(2.0*B);
    double sinlat  = sin(B);
    double sinlon  = sin(L);
    double coslon  = cos(L);
    Vector3d ptide;
    ptide[2] = -0.033 * sin2lat * (m1*coslon+m2*sinlon);
    ptide[1] = +0.009 * sinlat * (m1*sinlon-m2*coslon);
    ptide[0] = -0.009 * cos2lat * (m1*coslon+m2*sinlon);
    return ptide;
}

Vector3d tide_correction::ST1DIU(const Vector3d &XSTA, const Vector3d &XSUN, const Vector3d &XMON, double F2SUN, double F2MON)
{
    const double DHI = -0.0025;
    const double DLI = -0.0007;
    double RSTA =XSTA.norm();
    double SINPHI = XSTA[2] / RSTA;
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1]) / RSTA;
    double COS2PHI = COSPHI*COSPHI - SINPHI*SINPHI;
    double SINLA = XSTA[1] / COSPHI / RSTA;
    double COSLA = XSTA[0] / COSPHI / RSTA;
    double RMON = XMON.norm();
    double RSUN = XSUN.norm();

    double temp1 = RSUN*RSUN;
    double temp2 = RMON*RMON;
    double DRSUN = -3.0*DHI*SINPHI*COSPHI*F2SUN*XSUN[2]*(XSUN[0]*SINLA-XSUN[1]*COSLA)/temp1;
    double DRMON = -3.0*DHI*SINPHI*COSPHI*F2MON*XMON[2]*(XMON[0]*SINLA-XMON[1]*COSLA)/temp2;
    double DNSUN = -3.0*DLI*COS2PHI*F2SUN*XSUN[2]*(XSUN[0]*SINLA-XSUN[1]*COSLA)/temp1;
    double DNMON = -3.0*DLI*COS2PHI*F2MON*XMON[2]*(XMON[0]*SINLA-XMON[1]*COSLA)/temp2;
    double DESUN = -3.0*DLI*SINPHI*F2SUN*XSUN[2]*(XSUN[0]*COSLA+XSUN[1]*SINLA)/temp1;
    double DEMON = -3.0*DLI*SINPHI*F2MON*XMON[2]*(XMON[0]*COSLA+XMON[1]*SINLA)/temp2;

    double DR = DRSUN + DRMON;
    double DN = DNSUN + DNMON;
    double DE = DESUN + DEMON;

    Vector3d XCOSTA;
    XCOSTA[0] = DR*COSLA*COSPHI - DE*SINLA - DN*SINPHI*COSLA;
    XCOSTA[1] = DR*SINLA*COSPHI + DE*COSLA - DN*SINPHI*SINLA;
    XCOSTA[2] = DR*SINPHI + DN*COSPHI;
    return XCOSTA;
}

Vector3d tide_correction::ST1SEM(const Vector3d &XSTA, const Vector3d &XSUN, const Vector3d &XMON, double F2SUN, double F2MON)
{
    const double DHI = -0.0022;
    const double DLI = -0.0007;
    double RSTA = XSTA.norm();
    double SINPHI = XSTA[2] / RSTA;
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1]) / RSTA;
    double SINLA = XSTA[1] / COSPHI / RSTA;
    double COSLA = XSTA[0] / COSPHI / RSTA;
    double CTWOLA = COSLA*COSLA - SINLA*SINLA;
    double STWOLA = 2.0 * COSLA * SINLA;
    double RMON = XMON.norm();
    double RSUN = XSUN.norm();

    double temp1 = RSUN * RSUN;
    double temp2 = RMON * RMON;
    double DRSUN = -3.0/4.0*DHI*COSPHI*COSPHI*F2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*STWOLA-2.0*XSUN[0]*XSUN[1]*CTWOLA)/temp1;
    double DRMON = -3.0/4.0*DHI*COSPHI*COSPHI*F2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*STWOLA-2.0*XMON[0]*XMON[1]*CTWOLA)/temp2;
    double DNSUN = 3.0/2.0*DLI*SINPHI*COSPHI*F2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*STWOLA-2.0*XSUN[0]*XSUN[1]*CTWOLA)/temp1;
    double DNMON = 3.0/2.0*DLI*SINPHI*COSPHI*F2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*STWOLA-2.0*XMON[0]*XMON[1]*CTWOLA)/temp2;
    double DESUN = -3.0/2.0*DLI*COSPHI*F2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*CTWOLA+2.0*XSUN[0]*XSUN[1]*STWOLA)/temp1;
    double DEMON = -3.0/2.0*DLI*COSPHI*F2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*CTWOLA+2.0*XMON[0]*XMON[1]*STWOLA)/temp2;

    double DR = DRSUN + DRMON;
    double DN = DNSUN + DNMON;
    double DE = DESUN + DEMON;
    Vector3d XCOSTA;
    XCOSTA[0] = DR*COSLA*COSPHI - DE*SINLA - DN*SINPHI*COSLA;
    XCOSTA[1] = DR*SINLA*COSPHI + DE*COSLA - DN*SINPHI*SINLA;
    XCOSTA[2] = DR*SINPHI + DN*COSPHI;
    return XCOSTA;
}

Vector3d tide_correction::ST1L1(const Vector3d &XSTA, const Vector3d &XSUN, const Vector3d &XMON, double F2SUN, double F2MON)
{
    Vector3d XCOSTA;
    const double L1D = 0.00120;
    const double L1SD = 0.00240;
    double RSTA = XSTA.norm();
    double SINPHI = XSTA[2] / RSTA;
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
    double SINLA = XSTA[1]/COSPHI/RSTA;
    double COSLA = XSTA[0]/COSPHI/RSTA;
    double RMON = XMON.norm();
    double RSUN = XSUN.norm();

    double temp1 = RSUN * RSUN;
    double temp2 = RMON * RMON;
    //FOR THE DIURNAL BAND
    double L1 = L1D;
    double DNSUN = -L1*SINPHI*SINPHI*F2SUN*XSUN[2]*(XSUN[0]*COSLA+XSUN[1]*SINLA)/temp1;
    double DNMON = -L1*SINPHI*SINPHI*F2MON*XMON[2]*(XMON[0]*COSLA+XMON[1]*SINLA)/temp2;
    double DESUN = L1*SINPHI*(COSPHI*COSPHI-SINPHI*SINPHI)*F2SUN*XSUN[2]*(XSUN[0]*SINLA-XSUN[1]*COSLA)/temp1;
    double DEMON = L1*SINPHI*(COSPHI*COSPHI-SINPHI*SINPHI)*F2MON*XMON[2]*(XMON[0]*SINLA-XMON[1]*COSLA)/temp2;

    double DE = 3.0*(DESUN+DEMON);
    double DN = 3.0*(DNSUN+DNMON);
    XCOSTA[0] = -DE*SINLA-DN*SINPHI*COSLA;
    XCOSTA[1] = DE*COSLA-DN*SINPHI*SINLA;
    XCOSTA[2] = DN*COSPHI;
    //FOR THE SEMI-DIURNAL BAND
    L1 = L1SD;
    double CTWOLA = COSLA*COSLA - SINLA*SINLA;
    double STWOLA = 2.0*COSLA*SINLA;
    DNSUN = -L1/2.0*SINPHI*COSPHI*F2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*CTWOLA+2.0*XSUN[0]*XSUN[1]*STWOLA)/temp1;
    DNMON = -L1/2.0*SINPHI*COSPHI*F2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*CTWOLA+2.0*XMON[0]*XMON[1]*STWOLA)/temp2;
    DESUN = -L1/2.0*SINPHI*SINPHI*COSPHI*F2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*STWOLA-2.0*XSUN[0]*XSUN[1]*CTWOLA)/temp1;
    DEMON = -L1/2.0*SINPHI*SINPHI*COSPHI*F2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*STWOLA-2.0*XMON[0]*XMON[1]*CTWOLA)/temp2;
    DE = 3.0 * (DESUN+DEMON);
    DN = 3.0 * (DNSUN+DNMON);

    XCOSTA[0] = XCOSTA[0] - DE*SINLA - DN*SINPHI*COSLA;
    XCOSTA[1] = XCOSTA[1] + DE*COSLA - DN*SINPHI*SINLA;
    XCOSTA[2] = XCOSTA[2] + DN*COSPHI;
    return XCOSTA;
}

Vector3d tide_correction::STEP2DIU(const Vector3d &XSTA, double FHR, double T)
{
    Vector3d XCORSTA;
    double DATDI[31][9]={
    {-3.0, 0.0, 2.0, 0.0, 0.0,-0.010,-0.010, 0.00 , 0.00 },
    {-3.0, 2.0, 0.0, 0.0, 0.0,-0.010,-0.010, 0.00 , 0.00 },
    {-2.0, 0.0, 1.0,-1.0, 0.0,-0.020,-0.010, 0.00 , 0.00 },
    {-2.0, 0.0, 1.0, 0.0, 0.0,-0.080,-0.050,-0.010,-0.020},
    {-2.0, 2.0,-1.0, 0.0, 0.0,-0.020,-0.010, 0.00 , 0.00 },
    {-1.0, 0.0, 0.0,-1.0, 0.0,-0.100,-0.050, 0.00 ,-0.020},
    {-1.0, 0.0, 0.0, 0.0, 0.0,-0.510,-0.260,-0.020,-0.120},
    {-1.0, 2.0, 0.0, 0.0, 0.0, 0.010, 0.00 , 0.00 , 0.00 },
    { 0.0,-2.0, 1.0, 0.0, 0.0, 0.010, 0.00 , 0.00 , 0.00 },
    { 0.0, 0.0,-1.0, 0.0, 0.0, 0.020, 0.010, 0.0  , 0.0  },
    { 0.0, 0.0, 1.0, 0.0, 0.0, 0.060, 0.020, 0.000, 0.010},
    { 0.0, 0.0, 1.0, 1.0, 0.0, 0.010, 0.0  , 0.0 ,  0.0  },
    { 0.0, 2.0,-1.0, 0.0, 0.0, 0.010, 0.00 , 0.00 , 0.00 },
    { 1.0,-3.0, 0.0, 0.0, 1.0,-0.060, 0.00 , 0.00 , 0.00 },
    { 1.0,-2.0, 0.0,-1.0, 0.0, 0.010, 0.00 , 0.00 , 0.00 },
    { 1.0,-2.0, 0.0, 0.0, 0.0,-1.230,-0.050, 0.060,-0.060},
    { 1.0,-1.0, 0.0, 0.0,-1.0, 0.020, 0.00 , 0.00 , 0.00 },
    { 1.0,-1.0, 0.0, 0.0, 1.0, 0.040, 0.00 , 0.00 , 0.00 },
    { 1.0, 0.0, 0.0,-1.0, 0.0,-0.220, 0.010, 0.010, 0.00 },
    { 1.0, 0.0, 0.0, 0.0, 0.0,12.020,-0.450,-0.660, 0.170},
    { 1.0, 0.0, 0.0, 1.0, 0.0, 1.730,-0.070,-0.100, 0.020},
    { 1.0, 0.0, 0.0, 2.0, 0.0,-0.040, 0.00 , 0.00 , 0.00 },
    { 1.0, 1.0, 0.0, 0.0,-1.0,-0.500, 0.00 , 0.030, 0.00 },
    { 1.0, 1.0, 0.0, 0.0, 1.0, 0.010, 0.00 , 0.00 , 0.00 },
    { 0.0, 1.0, 0.0, 1.0,-1.0,-0.010, 0.00 , 0.00 , 0.00 },
    { 1.0, 2.0,-2.0, 0.0, 0.0,-0.010, 0.00 , 0.00 , 0.00 },
    { 1.0, 2.0, 0.0, 0.0, 0.0,-0.120, 0.010, 0.010, 0.00 },
    { 2.0,-2.0, 1.0, 0.0, 0.0,-0.010, 0.00 , 0.00 , 0.00 },
    { 2.0, 0.0,-1.0, 0.0, 0.0,-0.020, 0.020, 0.00 , 0.010},
    { 3.0, 0.0, 0.0, 0.0, 0.0, 0.00 , 0.010, 0.00 , 0.010},
    { 3.0, 0.0, 0.0, 1.0, 0.0, 0.00 , 0.010, 0.00 , 0.00 } };
    double TT   = T*T;
    double TTT  = T*T*T;
    double TTTT = T*T*T*T;
    double S    = 218.316645630+481267.881940*T-0.00146638890*TT+0.000001851390*TTT;
    double TAU  = FHR*15.0+280.46061840+36000.77005360*T+0.000387930*TT-0.00000002580*TTT-S;
    double PR   = 1.3969712780*T+0.0003088890*TT+0.0000000210*TTT+0.0000000070*TTTT;
    S +=PR;
    double H    = 280.466450+36000.76974890*T+0.000303222220*TT+0.0000000200*TTT-0.000000006540*TTTT;
    double P    = 83.353243120+4069.013635250*T-0.010321722220*TT-0.00001249910*TTT+0.000000052630*TTTT;
    double ZNS  = 234.955444990 +1934.136261970*T-0.002075611110*TT-0.000002139440*TTT+0.000000016500*TTTT;
    double PS   = 282.937340980+1.719457666670*T+0.000456888890*TT-0.000000017780*TTT-0.000000003340*TTTT;
    //REDUCE ANGLES TO BETWEEN 0 AND 360
    S   = fmod(fmod(S,360.0)+360.0,360.0);
    TAU = fmod(fmod(TAU,360.0)+360.0,360.0);
    H   = fmod(fmod(H,360.0)+360.0,360.0);
    P   = fmod(fmod(P,360.0)+360.0,360.0);
    ZNS = fmod(fmod(ZNS,360.0)+360.0,360.0);
    PS  = fmod(fmod(PS,360.0)+360.0,360.0);

    double RSTA =XSTA.norm();
    double SINPHI = XSTA[2]/RSTA;
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
    double SINLA = XSTA[1]/COSPHI/RSTA;
    double COSLA = XSTA[0]/COSPHI/RSTA;
    double ZLA = atan2(XSTA[1],XSTA[0]);
    for(int i=0;i<3;i++)
    {
        XCORSTA[i]=0.0;
    }
    double THETAF,DR,DN,DE;
    double COSSIN_2=COSPHI*COSPHI-SINPHI*SINPHI;
    const double D2R = 0.0174532925;
    for (int J=0;J<31;J++)
    {
        THETAF = (TAU+DATDI[J][0]*S+DATDI[J][1]*H+DATDI[J][2]*P+DATDI[J][3]*ZNS+DATDI[J][4]*PS)*D2R;
        DR = DATDI[J][5]*2.0*SINPHI*COSPHI*sin(THETAF+ZLA)+DATDI[J][6]*2.0*SINPHI*COSPHI*cos(THETAF+ZLA);
        DN = DATDI[J][7]*COSSIN_2*sin(THETAF+ZLA)+DATDI[J][8]*COSSIN_2*cos(THETAF+ZLA);
        DE = DATDI[J][7]*SINPHI*cos(THETAF+ZLA)+DATDI[J][8]*SINPHI*sin(THETAF+ZLA);
        XCORSTA[0]= XCORSTA[0]+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
        XCORSTA[1]= XCORSTA[1]+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
        XCORSTA[2]= XCORSTA[2]+DR*SINPHI+DN*COSPHI;
    }
    XCORSTA[0] /= 1000.0;
    XCORSTA[1] /= 1000.0;
    XCORSTA[2] /= 1000.0;
    return XCORSTA;
}

Vector3d tide_correction::ST2LON(const Vector3d &XSTA, double T)
{
    Vector3d XCOSTA;
    double DATDI[5][9] = {
        {0, 0, 0, 1, 0,  0.470, 0.230, 0.160, 0.070},
        {0, 2, 0, 0, 0, -0.200,-0.120,-0.110,-0.050},
        {1, 0,-1, 0, 0, -0.110,-0.080,-0.090,-0.040},
        {2, 0, 0, 0, 0, -0.130,-0.110,-0.150,-0.070},
        {2, 0, 0, 1, 0, -0.050,-0.050,-0.060,-0.030} };
    double TT   = T*T;
    double TTT  = T*T*T;
    double TTTT = T*T*T*T;
    double S = 218.316645630+481267.881940*T-0.00146638890*TT+0.000001851390*TTT;
    double PR = 1.3969712780*T+0.0003088890*TT+0.0000000210*TTT+0.0000000070*TTTT;
    S += PR;
    double H = 280.466450+36000.76974890*T+0.000303222220*TT+0.0000000200*TTT-0.000000006540*TTTT;
    double P = 83.353243120+4069.013635250*T-0.010321722220*TT-0.00001249910*TTT+0.000000052630*TTTT;
    double ZNS = 234.955444990 +1934.136261970*T-0.002075611110*TT-0.000002139440*TTT+0.000000016500*TTTT;
    double PS = 282.937340980+1.719457666670*T+0.000456888890*TT-0.000000017780*TTT-0.000000003340*TTTT;
    //REDUCE ANGLES TO BETWEEN 0 AND 360
    S = fmod(fmod(S,360.0)+360.0,360.0);
    H = fmod(fmod(H,360.0)+360.0,360.0);
    P = fmod(fmod(P,360.0)+360.0,360.0);
    ZNS = fmod(fmod(ZNS,360.0)+360,360);
    PS = fmod(fmod(PS,360.0)+360.0,360.0);
    double RSTA = XSTA.norm();
    double SINPHI = XSTA[2]/RSTA;
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
    double SINLA = XSTA[1]/COSPHI/RSTA;
    double COSLA = XSTA[0]/COSPHI/RSTA;
    XCOSTA[0] = 0.0;
    XCOSTA[1] = 0.0;
    XCOSTA[2] = 0.0;
    double THETAF,DR,DN,DE = 0.0;
    double SINPHI2 = SINPHI*SINPHI;
    const double D2R = 0.0174532925;
    for (int J=0;J<5;J++)
    {
        THETAF = (DATDI[J][0]*S+DATDI[J][1]*H+DATDI[J][2]*P+DATDI[J][3]*ZNS+DATDI[J][4]*PS)*D2R;
        DR = DATDI[J][5]*(3.0*SINPHI2-1.0)/2.0*cos(THETAF)+DATDI[J][7]*(3.0*SINPHI2-1.0)/2.0*sin(THETAF);
        DN = DATDI[J][6]*(COSPHI*SINPHI*2.0)*cos(THETAF)+DATDI[J][8]*(COSPHI*SINPHI*2.0)*sin(THETAF);

        XCOSTA[0] = XCOSTA[0]+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
        XCOSTA[1] = XCOSTA[1]+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
        XCOSTA[2] = XCOSTA[2]+DR*SINPHI+DN*COSPHI;
    }
    XCOSTA[0] /= 1000.0;
    XCOSTA[1] /= 1000.0;
    XCOSTA[2] /= 1000.0;
    return XCOSTA;
}

Vector3d tide_correction::ST2DIU(const Vector3d &XSTA, double FHR, double T)
{
    Vector3d XCOSTA;
    double DATDI[11][7] = {
        {-2.0, 0.0, 1.0, 0.0, 0.0,-0.090, 0.000},
        {-1.0, 0.0, 0.0,-1.0, 0.0,-0.100, 0.000},
        {-1.0, 0.0, 0.0, 0.0, 0.0,-0.530, 0.020},
        { 0.0, 0.0, 1.0, 0.0, 0.0, 0.060, 0.000},
        { 1.0,-3.0, 0.0, 0.0, 1.0,-0.050, 0.000},
        { 1.0,-2.0, 0.0, 0.0, 0.0,-1.230, 0.070},
        { 1.0, 0.0, 0.0,-1.0, 0.0,-0.220, 0.010},
        { 1.0, 0.0, 0.0, 0.0, 0.0,12.040,-0.720},
        { 1.0, 0.0, 0.0, 1.0, 0.0, 1.740,-0.100},
        { 1.0, 1.0, 0.0, 0.0,-1.0,-0.500, 0.030},
        { 1.0, 2.0, 0.0, 0.0, 0.0,-0.110, 0.010}  };
    double TT   = T*T;
    double TTT  = T*T*T;
    double TTTT = T*T*T*T;
    double S = 218.316645630+481267.881940*T-0.00146638890*TT+0.000001851390*TTT;
    double TAU = FHR*15.0+280.46061840+36000.77005360*T+0.000387930*TT-0.00000002580*TTT-S;
    double PR = 1.3969712780*T+0.0003088890*TT+0.0000000210*TTT+0.0000000070*TTTT;
    S += PR;

    double H = 280.466450+36000.76974890*T+0.000303222220*TT+0.0000000200*TTT-0.000000006540*TTTT;
    double P = 83.353243120+4069.013635250*T-0.010321722220*TT-0.00001249910*TTT+0.000000052630*TTTT;
    double ZNS = 234.955444990 +1934.136261970*T-0.002075611110*TT-0.000002139440*TTT+0.000000016500*TTTT;
    double PS = 282.937340980+1.719457666670*T+0.000456888890*TT-0.000000017780*TTT-0.000000003340*TTTT;

    S = fmod(S,360.0);
    TAU = fmod(TAU,360.0);
    H   = fmod(H,360.0);
    P   = fmod(P,360.0);
    ZNS = fmod(ZNS,360.0);
    PS  = fmod(PS,360.0);

    double RSTA = XSTA.norm();
    double SINPHI = XSTA[2]/RSTA;
    double COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
    double SINLA = XSTA[1]/COSPHI/RSTA;
    double COSLA = XSTA[0]/COSPHI/RSTA;
    double ZLA = atan2(XSTA[1],XSTA[0]);
    XCOSTA[0] = 0.0;
    XCOSTA[1] = 0.0;
    XCOSTA[2] = 0.0;
    double THETAF,DR,DN,DE;
    double COS_SIN = COSPHI*COSPHI-SINPHI*SINPHI;
    for(int J=0;J<11;J++)
    {
        THETAF = (TAU+DATDI[J][0]*S+DATDI[J][1]*H+DATDI[J][2]*P+DATDI[J][3]*ZNS+DATDI[J][4]*PS)*DR;
        DR = DATDI[J][5]*2.0*SINPHI*COSPHI*sin(THETAF+ZLA);
        DN = DATDI[J][6]*(COS_SIN)*sin(THETAF+ZLA);
        DE = DATDI[J][6]*SINPHI*cos(THETAF+ZLA);
        XCOSTA[0] = XCOSTA[0]+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
        XCOSTA[1] = XCOSTA[1]+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
        XCOSTA[2] = XCOSTA[2]+DR*SINPHI+DN*COSPHI;
    }

    XCOSTA /= 1000.0;
    return XCOSTA;
}

MatrixXd tide_correction::getArg(GPS_Time t)
{
    double DEG_TO_RAD=0.017453;
    int num_harmonics = 11;
    MatrixXd sig(1,num_harmonics);
    sig(0,0)  = 1.40519e-4;
    sig(0,1)  = 1.45444e-4;
    sig(0,2)  = 1.37880e-4;
    sig(0,3)  = 1.45842e-4;
    sig(0,4)  = 0.72921e-4;
    sig(0,5)  = 0.67598e-4;
    sig(0,6)  = 0.72523e-4;
    sig(0,7)  = 0.64959e-4;
    sig(0,8)  = 0.053234e-4;
    sig(0,9)  = 0.026392e-4;
    sig(0,10) = 0.003982e-4;
    MatrixXd angfac(4,11);
    angfac<< 2.0,-2.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0,
             2.0,-3.0, 1.0, 0.0,
             2.0, 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0, 0.25,
             1.0,-2.0, 0.0,-0.25,
            -1.0, 0.0, 0.0,-0.25,
             1.0,-3.0, 1.0,-0.25,
             0.0, 2.0, 0.0, 0.0,
             0.0, 1.0,-1.0, 0.0,
             2.0, 0.0, 0.0, 0.0;
    MatrixXd arguments(1,num_harmonics);
    dayofYear time;
    gregorianTime Time;
    timeTransform::GPSTToGregorian(t,Time);
    timeTransform::gregorianToDOY(Time,time);
    double year = time.year;
    double fday = time.sod;
    double d  = time.DOY+365.0*(year-1975.0)+floor((year-1973.0)/4.0);
    double tt = (27392.500528+1.000000035*d)/36525.0;
    double H0 = (279.69668+(36000.768930485+3.03e-4*tt)*tt)*DEG_TO_RAD;
    double S0 = (((1.9e-6*tt-0.001133)*tt+481267.88314137)*tt+270.434358)*DEG_TO_RAD;
    double P0 = (((-1.2e-5*tt-0.010325)*tt +4069.0340329577)*tt+334.329653)*DEG_TO_RAD;

    for (int k=0;k<num_harmonics;k++)
    {
        double temp = sig(0,k)*fday + angfac(0,k)*H0 + angfac(1,k)*S0 + angfac(2,k)*P0 + angfac(3,k)*2*pi;
        arguments(0,k) = fmod(temp,2*pi);
        if(arguments(0,k) < 0.0)
        {
            arguments(0,k) = arguments(0,k)+2*pi;
        }
    }
    return arguments;
}
