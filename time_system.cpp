#include "time_system.h"

const double GPS_Time::JDGPSS = 2444244.5;

time_Transform::time_Transform()
{

}

time_Transform::~time_Transform()
{

}

julian_Time time_Transform::GCtoJulian(const GC_Time &GC)
{
    double Second_of_Day = GC.getHour()*3600+GC.getMinute()*60+GC.getSecond();
    double julian = 0;
    if( GC.getMonth() > 2)
    {
        julian = (int)(365.25 * GC.getYear()) + (int)(30.6001 * (GC.getMonth() + 1)) + GC.getDay() +
                Second_of_Day / (24*3600) + 1720981.5;
    }
    else
    {
        julian = (int)(365.25 * (GC.getYear()-1)) + (int)(30.6001 * (GC.getMonth()+13)) + GC.getDay() +
                         Second_of_Day / (24*3600) + 1720981.5;
    }
    julian_Time Julian(julian);
    return Julian;
}

GC_Time time_Transform::JuliantoGC(const julian_Time &Julian)
{
    GC_Time GC;
    double julian = Julian.getJulian();
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

julian_Time time_Transform::GPSTtoJulian(const GPS_Time &GPST)
{
    double julian = GPST.getGPSW() * 7 + GPST.getGPSS() / (60 * 60 * 24) + 2444244.5;
    julian_Time Julian(julian);
    return Julian;
}

GPS_Time time_Transform::JuliantoGPST(const julian_Time &Julian)
{
    double julian = Julian.getJulian();
    double GPSW = (int)((julian-JDGPSS)/7);
    double GPSS = ((julian-JDGPSS) - GPSW * 7 ) * 24 * 3600;
    GPS_Time gps(GPSW,GPSS);
    return gps;
}

GPS_Time time_Transform::GCtoGPST(const GC_Time &GC)
{
    return time_Transform::JuliantoGPST( time_Transform::GCtoJulian(GC) );
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

void GC_Time::getGC_Time(int &year, int &month, int &day, int &hour, double &minute, double &second) const
{
    year    = getYear();
    month   = getMonth();
    day     = getDay();
    hour    = getHour();
    minute  = getMinute();
    second  = getSecond();
}

double GC_Time::operator -(const GC_Time &right) const
{
    return time_Transform::GCtoJulian(*this) - time_Transform::GCtoJulian(right);
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
    if(time_Transform::GCtoJulian(*this) > time_Transform::GCtoJulian(right))
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
    if(time_Transform::GCtoJulian(*this) < time_Transform::GCtoJulian(right))
    {
        return true;
    }
    else
    {
        return false;
    }
}
int GC_Time::getYear() const
{
    return year;
}

void GC_Time::setYear(int value)
{
    if(value <= 0)
    {
        std::cerr << "You want set a error year for GC time! " << endl
                  << "That value is " << value << endl;
        exit(1);
    }
    else
    {
        year = value;
    }

}
int GC_Time::getMonth() const
{
    return month;
}

void GC_Time::setMonth(int value)
{
    if(value <= 0 || value > 12)
    {
        std::cerr << "You want set a error month for GC time! " << endl
                  << "That value is " << value << endl;
        exit(1);
    }
    else
    {
        month = value;
    }
}
int GC_Time::getDay() const
{
    return day;
}

void GC_Time::setDay(int value)
{
    bool leap_year = false;
    /* -------------------------------------------------------------------------
     * If month is 1,3,5,7,8,10,12 month_flag = true
     * If month is 4,6,9,11 month_flag = false
     * -------------------------------------------------------------------------*/
    bool month_flag = false;
    if(month <= 0 || year <= 0)
    {
        std::cerr << "please put year and month in GC time!" << endl;
        exit(1);
    }
    if(year%400 == 0)
    {
        leap_year = true;
    }
    else
    {
        if((year%4 == 0) && (year%100!=0))
        {
            leap_year = true;
        }
    }
    switch (month) {
        case 1:  month_flag = true; break;
        case 3:  month_flag = true; break;
        case 5:  month_flag = true; break;
        case 7:  month_flag = true; break;
        case 8:  month_flag = true; break;
        case 10: month_flag = true; break;
        case 12: month_flag = true; break;
        case 4:  month_flag = false; break;
        case 6:  month_flag = false; break;
        case 9:  month_flag = false; break;
        case 11: month_flag = false; break;
        case 2:
        {
            if(leap_year == true && value > 29)
            {
                std::cerr << "You want set a error day  for GC time! " << endl
                          << "That value is "
                          << "year = " << year <<"leap_year = true" << endl
                          << "month = " << month << endl;
                          << "day = " << value << endl;
                exit(1);
            }
            else if(leap_year == false && value > 28)
            {
                std::cerr << "You want set a error day  for GC time! " << endl
                          << "That value is "
                          << "year = " << year <<"leap_year = false" << endl
                          << "month = " << month << endl;
                          << "day = " << value << endl;
                exit(1);
            }
            else
            {
                day = value;
            }
        }
    }
    if(month_flag == true && day > 31)
    {
        std::cerr << "You want set a error day  for GC time! " << endl
                  << "month = " << month << endl;
                  << "day = " << value << endl;
        exit(1);
    }
    else if(month_flag == false && day > 30)
    {
        std::cerr << "You want set a error day  for GC time! " << endl
                  << "month = " << month << endl;
                  << "day = " << value << endl;
        exit(1);
    }
    else
    {
        day = value;
    }
}
int GC_Time::getHour() const
{
    return hour;
}

void GC_Time::setHour(int value)
{
    if(value < 0 || value > 60)
    {
        std::cerr << "You want set a error hour for GC time! " << endl
                  << "That value is " << value << endl;
        exit(1);
    }
    else
    {
        hour = value;
    }
}
int GC_Time::getMinute() const
{
    return minute;
}

void GC_Time::setMinute(int value)
{
    if(value < 0 || value > 60)
    {
        std::cerr << "You want set a error minute for GC time! " << endl
                  << "That value is " << value << endl;
        exit(1);
    }
    else
    {
        minute = value;
    }
}
double GC_Time::getSecond() const
{
    return second;
}

void GC_Time::setSecond(double value)
{
    if(value < 0 || value > 60.0)
    {
        std::cerr << "You want set a error second for GC time! " << endl
                  << "That value is " << value << endl;
        exit(1);
    }
    else
    {
        second = value;
    }
}








julian_Time::julian_Time()
    :julian(0.0)
{

}

julian_Time::julian_Time(double value)
{
    setJulian(value);
}

double julian_Time::getJulian() const
{
    return julian;
}

void julian_Time::setJulian(double value)
{
    if(value < 0)
    {
        std::cerr << "You want set a error Julian time! " << endl
                  << "That value is " << value << endl;
        exit(1);
    }
    else
    {
        julian = value;
    }

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
    return time_Transform::GPSTtoJulian(*this) - time_Transform::GPSTtoJulian(right);
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
    if(time_Transform::GPSTtoJulian(*this) > time_Transform::GPSTtoJulian(right))
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
    if(time_Transform::GPSTtoJulian(*this) < time_Transform::GPSTtoJulian(right))
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
    if(GPSW < 0)
    {
        std::cerr << "You want set a error GPS time! " << endl;
                  << "That value is "
                  << "GPSW = " << GPSW << endl;
                  << "GPSS = " << GPSS << endl;
        exit(1);
    }
    else
    {
        if(GPSS < 0)
        {
            std::cerr << "You want set a error GPS time! " << endl;
                      << "That value is "
                      << "GPSW = " << GPSW << endl;
                      << "GPSS = " << GPSS << endl;
            exit(1);
        }
        else
        {
            this->GPSW = GPSW;
            this->GPSS = GPSS;
        }
    }

}

int GPS_Time::getGPSW() const
{
    return GPSW;
}

double GPS_Time::getGPSS() const
{
    return GPSS;
}



void DOY_Time::setDOY(int year, int doy, double sod)
{

}
