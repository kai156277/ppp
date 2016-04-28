#include "time_system.h"

using namespace std;

const double GPS_Time::JDGPSS = 2444244.5;

time_Transform::time_Transform()
{

}

time_Transform::~time_Transform()
{

}

julian_Time time_Transform::GCtoJulian(const GC_Time &GC)
{
    double Second_of_Day = GC.hour*3600+GC.minute*60+GC.second;
    double julian = 0;
    if( GC.month > 2)
    {
        julian = (int)(365.25 * GC.year) + (int)(30.6001 * (GC.month + 1)) + GC.day +
                Second_of_Day / (24*3600) + 1720981.5;
    }
    else
    {
        julian = (int)(365.25 * (GC.year-1)) + (int)(30.6001 * (GC.month+13)) + GC.day +
                         Second_of_Day / (24*3600) + 1720981.5;
    }
    julian_Time Julian(julian);
    return Julian;
}

GC_Time time_Transform::JuliantoGC(const julian_Time &Julian)
{
    GC_Time GC;
    double julian = Julian.julian;
    int a = (int)(julian + 0.5);
    int b = a + 1537;
    int c = (int)((b - 122.1) / 365.25);
    int d = (int)(365.25 * c);
    int e = (int)((b - d) / 30.6001);
    int hour   = 24 * (julian + 0.5 - (int)(julian + 0.5));
    int minute = (int)(60 * (24 * (julian + 0.5 - (int)(julian + 0.5)) - hour));
    double second = 60 * (60 * (24 * (julian + 0.5 - (int)(julian + 0.5)) - hour) - minute);
    int day    = b - d -  (int)(30.6001 * e);
    int month  = e - 1 - 12 * (int)(e / 14);
    int year   = c - 4715 - (int)((7 + month) / 10);
    GC.setGC_Time(year,month,day,hour,minute,second);
    return GC;
}

GC_Time time_Transform::GPSTtoGC(const GPS_Time &GPST)
{
    return time_Transform::JuliantoGC(time_Transform::GPSTtoJulian(GPST));
}

GC_Time time_Transform::DOYtoGC(const DOY_Time &DOY)
{
    return time_Transform::JuliantoGC(time_Transform::DOYtoJulian(DOY));
}

julian_Time time_Transform::GPSTtoJulian(const GPS_Time &GPST)
{
    double julian = GPST.GPSW * 7 + GPST.GPSS / (60 * 60 * 24) + 2444244.5;
    julian_Time Julian(julian);
    return Julian;
}

julian_Time time_Transform::DOYtoJulian(const DOY_Time &DOY)
{
    GC_Time GC(DOY.year,1,1,0,0,0);
    julian_Time julian = time_Transform::GCtoJulian(GC);
    double doy = julian.julian + DOY.doy + DOY.sod;
    julian_Time DOY_julian(doy);
    return DOY_julian;
}

GPS_Time time_Transform::JuliantoGPST(const julian_Time &Julian)
{
    double julian = Julian.julian;
    double GPSW = (int)((julian - GPS_Time::JDGPSS)/7);
    double GPSS = ((julian - GPS_Time::JDGPSS) - GPSW * 7 ) * 24 * 3600;
    GPS_Time gps(GPSW,GPSS);
    return gps;
}

GPS_Time time_Transform::GCtoGPST(const GC_Time &GC)
{
    return time_Transform::JuliantoGPST( time_Transform::GCtoJulian(GC) );
}

GPS_Time time_Transform::DOYtoGPST(const DOY_Time &DOY)
{
    return time_Transform::JuliantoGPST( time_Transform::DOYtoJulian(DOY));
}

DOY_Time time_Transform::GCtoDOY(const GC_Time &GC)
{
    return time_Transform::JuliantoDOY( time_Transform::GCtoJulian(GC) );
}

DOY_Time time_Transform::GPSTtoDOY(const GPS_Time &GPST)
{
    return time_Transform::JuliantoDOY( time_Transform::GPSTtoJulian(GPST) );
}

DOY_Time time_Transform::JuliantoDOY(const julian_Time &Julian)
{
    GC_Time GC = time_Transform::JuliantoGC(Julian);
    int year = GC.year;
    int hour = GC.hour;
    int minute = GC.minute;
    double second = GC.second;
    double sod = hour * 3600 + minute * 60 + second;
    GC_Time doy_GC(year,1,1,hour,minute,second);
    int doy = GC - doy_GC + 1;
    DOY_Time DOY(year,doy,sod);
    return DOY;
}


GC_Time::GC_Time()
    :year(0),month(0),day(0),
    hour(0),minute(0),second(0.0)
{

}

GC_Time::GC_Time(int year, int month, int day, int hour, double minute, double second)
{
    setGC_Time(year,month,day,hour,minute,second);
}

GC_Time::~GC_Time()
{

}

void GC_Time::setGC_Time(int year, int month, int day, int hour, double minute, double second)
{
    this->year   = year;
    this->month  = month;
    this->day    = day;
    this->hour   = hour;
    this->minute = minute;
    this->second = second;
}

double GC_Time::operator -(const GC_Time &right) const
{
    julian_Time this_julian =  time_Transform::GCtoJulian(*this);
    julian_Time right_julian = time_Transform::GCtoJulian(right);
    return  this_julian.julian - right_julian.julian ;
}

bool GC_Time::operator ==(const GC_Time &right) const
{
    if((*this - right) < 1e-6)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool GC_Time::operator >(const GC_Time &right) const
{
    julian_Time this_julian =  time_Transform::GCtoJulian(*this);
    julian_Time right_julian = time_Transform::GCtoJulian(right);
    if(this_julian.julian > right_julian.julian)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool GC_Time::operator <(const GC_Time &right) const
{
    julian_Time this_julian =  time_Transform::GCtoJulian(*this);
    julian_Time right_julian = time_Transform::GCtoJulian(right);
    if(this_julian.julian > right_julian.julian)
    {
        return true;
    }
    else
    {
        return false;
    }
}

julian_Time::julian_Time()
    :julian(0.0)
{

}

julian_Time::julian_Time(double value)
{
    julian = value;
}

julian_Time::~julian_Time()
{

}


GPS_Time::GPS_Time()
{

}

GPS_Time::GPS_Time(int GPSW, double GPSS)
{
    setGPS_Time(GPSW,GPSS);
}

GPS_Time::~GPS_Time()
{

}

double GPS_Time::operator -(const GPS_Time &right) const
{
    julian_Time this_julian = time_Transform::GPSTtoJulian(*this);
    julian_Time right_julian = time_Transform::GPSTtoJulian(right);
    return  this_julian.julian - right_julian.julian;
}

bool GPS_Time::operator ==(const GPS_Time &right) const
{
    if((*this - right) < 1e-6)
    {
        return true;
    }
    else
    {
        return false;
    }

}

bool GPS_Time::operator >(const GPS_Time &right) const
{
    julian_Time this_julian = time_Transform::GPSTtoJulian(*this);
    julian_Time right_julian = time_Transform::GPSTtoJulian(right);
    if(this_julian.julian > right_julian.julian)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool GPS_Time::operator <(const GPS_Time &right) const
{
    julian_Time this_julian = time_Transform::GPSTtoJulian(*this);
    julian_Time right_julian = time_Transform::GPSTtoJulian(right);
    if(this_julian.julian < right_julian.julian)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void GPS_Time::setGPS_Time(int GPSW, double GPSS)
{
    this->GPSW = GPSW;
    this->GPSS = GPSS;
}



DOY_Time::DOY_Time()
    :year(0),doy(0),sod(0)
{

}

DOY_Time::DOY_Time(int year, int doy, double sod)
{
    setDOY(year,doy,sod);
}

DOY_Time::~DOY_Time()
{

}

void DOY_Time::setDOY(int year, int doy, double sod)
{
    this->year = year;
    this->doy = doy;
    this->sod = sod;
}


