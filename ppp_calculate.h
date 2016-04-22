#ifndef PPP_CALCULATE_H
#define PPP_CALCULATE_H

#include<QVector>
#include<QString>

#include"o_date.h"
#include"sp3_date.h"
#include"clock_date.h"
#include"ppp_date.h"
#include"snx_date.h"

class ppp_calculate
{
public:
    ppp_calculate();
    void ppp_coordinate_clock(const o_file_date &ofile,const sp3_file &sp3file,const clock_file &clockfile,ppp_file &ppp);
    void set_station_coordinate(const snx_date &snx);
private:
    void sate_angle(ppp_sate &date);
    void sate_sagnac(ppp_sate &date);
    void sate_relativity(ppp_sate &date);
    static double c;
    static double w; //sagnac
    static double u; //relativity
    double station_x;
    double station_y;
    double station_z;
};

#endif // PPP_CALCULATE_H
