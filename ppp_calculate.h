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
    void sate_troposphere(ppp_sate &date,int doy);
    double station_x;
    double station_y;
    double station_z;
    double station_B;
    double station_L;
    double station_H;
    const static double c;
    const static double w; //sagnac
    const static double u; //relativity
    const static double Pi;

    /*UNB3参数-------------------------------------------------------------------*/
    const static double UNB3_average[5][6];
    const static double UNB3_amplitude[5][6];
    const static double Neil_wet_average[5][4];
    const static double R;
    const static double k1;
    const static double k2;
    const static double k3;
    /*--------------------------------------------------------------------------*/
    const static double Neil_hyd_average[5][4];
    const static double Neil_hyd_amplitude[5][4];
    /*----------------------------------------------------------------------------*/
};

#endif // PPP_CALCULATE_H
