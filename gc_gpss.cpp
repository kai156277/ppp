#include "gc_gpss.h"

double GC_GPSS::JDGPSS = 2444244.5;

GC_GPSS::GC_GPSS()
    :JD(0),
     Year(0),Month(0),Day(0),
     Hour(0),Minute(0),Second(0),
     GPSW(0),GPSS(0)
{

}

void GC_GPSS::setJD(double jd)
{
    JD = jd;
}

void GC_GPSS::setGPS(double gpss, int gpsw)
{
    GPSW = gpsw;
    GPSS = gpss;
}

void GC_GPSS::setGC(int year, int month, int day,
                               double hour, double minute, double second)
{
    Hour = hour;
    Minute = minute;
    Second = second;
    Year = year;
    Month = month;
    Day = day;
}

double GC_GPSS::getJD()
{
    return JD;
}

double GC_GPSS::getGPSS()
{
    return GPSS;
}

void GC_GPSS::clear()
{
    JD = 0;
    Year = 0;Month = 0;Day = 0;
    Hour = 0;Minute = 0;Second = 0;
    GPSW = 0;GPSS = 0;
}

int GC_GPSS::getGPSW()
{
    return GPSW;
}

void GC_GPSS::GCtoGPS()
{
    if( Month > 2)
        {
             JD=static_cast<int>(365.25*Year)+static_cast<int>(30.6001*(Month+1))+Day+(Hour*3600+Minute*60+Second)/(24*3600)+1720981.5;
        }
    else
        {
            JD=static_cast<int>(365.25*(Year-1))+static_cast<int>(30.6001*(Month+13))+Day+(Hour*3600+Minute*60+Second)/(24*3600)+1720981.5;
        }
        GPSW=static_cast<int>((JD-JDGPSS)/7);
        GPSS=((JD-JDGPSS) - GPSW * 7 ) * 24 * 3600;
}

void GC_GPSS::GPStoGC()
{
    int a,b,c,d,e;
    JD=GPSS/(24*3600)+GPSW*7+JDGPSS;
    a=static_cast<int>(JD+0.5);
    b=a+1537;
    c=static_cast<int>((b-122.1)/365.25);
    d=static_cast<int>(365.25*c);
    e=static_cast<int>((b-d)/30.6001);
    Day=b-d-static_cast<int>(30.6001*e)+(JD + 0.5 - static_cast<int>(JD+0.5));
    Month=e-1-12*static_cast<int>(e/14);
    Year=c-4715-static_cast<int>((7+Month)/10);
}
GC_GPSS::~GC_GPSS()
{

}

