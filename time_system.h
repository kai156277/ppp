#ifndef TIME_SYSTEM_H
#define TIME_SYSTEM_H

#include <ostream>

class julian_Time
{
public:
    julian_Time();
    julian_Time(double value);
    ~julian_Time();
    double getJulian() const;
    void setJulian(double value);

private:
    double julian;
};

class GC_Time
{
public:
    GC_Time();
    GC_Time(int year, int month, int day, int hour, double minute, double second);
    ~GC_Time();

    void setGC_Time(int year, int month, int day, int hour, double minute, double second);
    void getGC_Time(int &year, int &month, int &day, int &hour, double &minute, double &second) const;

    double operator -(const GC_Time &right) const;
    bool operator ==(const GC_Time &right) const;
    bool operator >(const GC_Time &right) const;
    bool operator <(const GC_Time &right) const;

    int getYear() const;
    void setYear(int value);
    int getMonth() const;
    void setMonth(int value);
    int getDay() const;
    void setDay(int value);
    int getHour() const;
    void setHour(int value);
    int getMinute() const;
    void setMinute(int value);
    double getSecond() const;
    void setSecond(double value);

private:
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
    int getGPSW() const;
    double getGPSS() const;

private:
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
private:
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
    static GC_Time JuliantoGC(const julian_Time &Julian);
    static GC_Time GPSTtoGC(const GPS_Time &GPST);
    static GPS_Time JuliantoGPST(const julian_Time &Julian);
    static GPS_Time GCtoGPST(const GC_Time &GC);

private:
    static const double JDGPSS;
};

#endif // TIME_SYSTEM_H
