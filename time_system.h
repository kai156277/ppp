#ifndef TIME_SYSTEM_H
#define TIME_SYSTEM_H

#include <iostream>

class julian_Time
{
public:
    julian_Time();
    julian_Time(double value);
    ~julian_Time();

    double julian;
};

class GC_Time
{
public:
    GC_Time();
    GC_Time(int year, int month, int day, int hour, double minute, double second);
    ~GC_Time();

    void setGC_Time(int year, int month, int day, int hour, double minute, double second);

    double operator -(const GC_Time &right) const;
    bool operator ==(const GC_Time &right) const;
    bool operator >(const GC_Time &right) const;
    bool operator <(const GC_Time &right) const;

    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
};

class GPS_Time
{
public:
    GPS_Time();
    GPS_Time(int GPSW,double GPSS);
    ~GPS_Time();

    double operator -(const GPS_Time &right) const;
    bool operator ==(const GPS_Time &right) const;
    bool operator >(const GPS_Time &right) const;
    bool operator <(const GPS_Time &right) const;
    void setGPS_Time(int GPSW,double GPSS);


    const static double JDGPSS;
    int GPSW;
    double GPSS;

};

class DOY_Time
{
public:
    DOY_Time();
    DOY_Time(int year,int doy,double sod);
    ~DOY_Time();
    void setDOY(int year,int doy,double sod);
    double operator -(const DOY_Time &right) const;
    bool operator ==(const DOY_Time &right) const;
    bool operator >(const DOY_Time &right) const;
    bool operator <(const DOY_Time &right) const;

    int year;
    int doy;
    double sod;
};

class time_Transform
{
public:
    time_Transform();
    ~time_Transform();
    static julian_Time GCtoJulian(const GC_Time &GC);
    static julian_Time GPSTtoJulian(const GPS_Time &GPST);
    static julian_Time DOYtoJulian(const DOY_Time &DOY);

    static GC_Time JuliantoGC(const julian_Time &Julian);
    static GC_Time GPSTtoGC(const GPS_Time &GPST);
    static GC_Time DOYtoGC(const DOY_Time &DOY);

    static GPS_Time JuliantoGPST(const julian_Time &Julian);
    static GPS_Time GCtoGPST(const GC_Time &GC);
    static GPS_Time DOYtoGPST(const DOY_Time &DOY);

    static DOY_Time GCtoDOY(const GC_Time &GC);
    static DOY_Time GPSTtoDOY(const GPS_Time &GPST);
    static DOY_Time JuliantoDOY(const julian_Time &Julian);

};

#endif // TIME_SYSTEM_H
