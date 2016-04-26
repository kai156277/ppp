#ifndef GC_GPSS_H
#define GC_GPSS_H


class GC_GPSS
{
public:
    GC_GPSS();
    ~GC_GPSS();
    double JD;
    double Hour;
    double Minute;
    double Second;
    int Year;
    int Month;
    int Day;
    int GPSW;
    double GPSS;
    int DOY;
    static double JDGPSS;
    void setGC(int year, int month, int day,
               double hour, double minute, double second);
    void setGPS(double gpss, int gpsw);
    void setJD(double jd);
    void GCtoGPS();
    void GPStoGC();
    void GCtoDOY();
    void clear();
};

#endif // GC_GPSS_H
