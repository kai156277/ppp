#include <algorithm>
#include <iostream>
#include <cmath>

#include "ppp_calculate.h"
#include "gc_gpss.h"
#include "math_function.h"
#include "coordinate.h"

double ppp_calculate::c = 299792458;
double ppp_calculate::w = 7.292115e-5;
double ppp_calculate::u = 3.986004418e+14;


ppp_calculate::ppp_calculate()
{

}

void ppp_calculate::ppp_coordinate_clock(const o_file_date &ofile,const sp3_file &sp3file,const clock_file &clockfile, ppp_file &ppp)
{
    //钟差
    double interval = clockfile.file[1].GPSS - clockfile.file[0].GPSS;
    double clock_interval = interval * 7 / 2.0;
    //坐标
    double sp3_interval = sp3file.heard.interval * 7 / 2.0;
    for(int i = 0; i<ofile.satellite_file.size(); i++)
    {
        ppp_epoch epoch;
        epoch.year   = ofile.satellite_file[i].year;
        epoch.month  = ofile.satellite_file[i].month;
        epoch.day    = ofile.satellite_file[i].day;
        epoch.hour   = ofile.satellite_file[i].hour;
        epoch.minute = ofile.satellite_file[i].minute;
        epoch.second = ofile.satellite_file[i].second;
        epoch.sate_number = ofile.satellite_file[i].number_of_satellite;
        epoch.GPSS = ofile.satellite_file[i].GPSS;

        //寻找符合历元的clock lagrange的插值起点
        double clock_find_time = epoch.GPSS - clock_interval;
        QVector<clock_epoch>::const_iterator clock_find_time_item = NULL;
        clock_find_time_item = std::find_if(clockfile.file.begin(),clockfile.file.end(),clock_GPSS_finder(clock_find_time));
        int clock_item = clock_find_time_item - clockfile.file.begin(); //lagrange的插值起点

        //寻找符合历元的坐标lagrange的插值起点
        double sp3_find_time = epoch.GPSS - sp3_interval;
        QVector<sp3_epoch>::const_iterator sp3_find_time_item = NULL;
        sp3_find_time_item = std::find_if(sp3file.file.begin(),sp3file.file.end(),sp3_GPSS_finder(sp3_find_time));
        int sp3_item = sp3_find_time_item - sp3file.file.begin();

        //寻找每一个对应卫星的lagrange的插值
        for(int j = 0; j<ofile.satellite_file[i].satellite_epoch.size(); j++)
        {
            ppp_sate sate;
            //每一个卫星的7个插值
            double sate_clock[14] = {0};    //前七个数值为GPS时间，后七个数值为钟差值
            double sate_x[14] = {0};    //前七个数值为GPS时间，后七个数值为坐标值
            double sate_y[14] = {0};
            double sate_z[14] = {0};
            double velocity_x[14] = {0};
            double velocity_y[14] = {0};
            double velocity_z[14] = {0};

            sate.P1 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[0];
            sate.P2 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[1];
            sate.P3 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[2];
            sate.L1 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[3];
            sate.L2 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[4];
            sate.L3 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[5];


            //计算钟差
            for(int k = 0; k<7; k++)
            {
                QVector<clock_info>::const_iterator sate_find = NULL;
                sate_find = std::find(clockfile.file[clock_item+k].GPS_epoch.begin(),
                                      clockfile.file[clock_item+k].GPS_epoch.end(),
                                      ofile.satellite_file[i].satellite_epoch[j]);
                sate_clock[k] = clockfile.file[clock_item+k].GPSS;
                sate_clock[k+7] = sate_find->clock_bias;
            }
            double clock = math_function::lagrange(epoch.GPSS,sate_clock,14);
            double clock_error = clock + sate.P2 / ppp_calculate::c;
            //计算坐标
            for(int k = 0; k<7; k++)
            {
                QVector<sp3_sate>::const_iterator sate_find = NULL;
                sate_find = std::find(sp3file.file[sp3_item+k].GPS_epoch.begin(),
                                      sp3file.file[sp3_item+k].GPS_epoch.end(),
                                      ofile.satellite_file[i].satellite_epoch[j]);
                sate_x[k] = sp3file.file[sp3_item+k].GPSS + clock_error ;
                sate_y[k] = sp3file.file[sp3_item+k].GPSS + clock_error ;
                sate_z[k] = sp3file.file[sp3_item+k].GPSS + clock_error ;
                sate_x[k+7] = sate_find->x*1000;
                sate_y[k+7] = sate_find->y*1000;
                sate_z[k+7] = sate_find->z*1000;
            }
            sate.position_x = math_function::lagrange(epoch.GPSS,sate_x,14);
            sate.position_y = math_function::lagrange(epoch.GPSS,sate_y,14);
            sate.position_z = math_function::lagrange(epoch.GPSS,sate_z,14);

            //计算速度
            for(int k = 0; k<7; k++)
            {
                QVector<sp3_sate>::const_iterator sate_find = NULL;
                sate_find = std::find(sp3file.file[sp3_item+k].GPS_epoch.begin(),
                                      sp3file.file[sp3_item+k].GPS_epoch.end(),
                                      ofile.satellite_file[i].satellite_epoch[j]);
                velocity_x[k] = sp3file.file[sp3_item+k].GPSS + clock_error - 0.5;
                velocity_y[k] = sp3file.file[sp3_item+k].GPSS + clock_error - 0.5;
                velocity_z[k] = sp3file.file[sp3_item+k].GPSS + clock_error - 0.5;
                velocity_x[k+7] = sate_find->x*1000;
                velocity_y[k+7] = sate_find->y*1000;
                velocity_z[k+7] = sate_find->z*1000;
            }
            sate.velocity_x = (sate.position_x - math_function::lagrange(epoch.GPSS,velocity_x,14)) / 0.5;
            sate.velocity_y = (sate.position_y - math_function::lagrange(epoch.GPSS,velocity_y,14)) / 0.5;
            sate.velocity_z = (sate.position_z - math_function::lagrange(epoch.GPSS,velocity_z,14)) / 0.5;


            sate.PRN = ofile.satellite_file[i].satellite_epoch[j].satellite_infomation;
            sate.clock = clock * ppp_calculate::c;

            sate_angle(sate);       //计算卫星角度信息
            sate_sagnac(sate);      //计算地球自转效应
            sate_relativity(sate);  //计算相对论效应
            epoch.sate_info.push_back(sate);
        }
        ppp.file.push_back(epoch);
    }
}

void ppp_calculate::set_station_coordinate(const snx_date &snx)
{
    station_x = snx.station_x;
    station_y = snx.station_y;
    station_z = snx.station_z;
}

void ppp_calculate::sate_angle(ppp_sate &date)
{
    Coordinate angle;
    angle.setENU(date.position_x, date.position_y, date.position_z ,
                 this->station_x, this->station_y, this->station_z).ENUparameter();
    date.azimuth   = angle.azimuth;
    date.elevation = angle.elevation;
}

void ppp_calculate::sate_sagnac(ppp_sate &date)
{
    date.sagnac = -( ppp_calculate::w/ppp_calculate::c )
                * ((date.position_x - this->station_x)*date.position_y -
                   (date.position_y - this->station_y)*date.position_x);
}

void ppp_calculate::sate_relativity(ppp_sate &date)
{
    double Relativity = 0, Relativity1 = 0;
    double dj = 0,di = 0,dji = 0;
    dj = sqrt( pow( date.position_x,2 ) + pow( date.position_y,2 ) + pow( date.position_z,2 ));
    di = sqrt( pow( this->station_x,2 ) + pow( this->station_y,2 ) +pow( this->station_z,2));
    dji = sqrt( pow(( date.position_x - this->station_x ),2 ) +
              pow(( date.position_y - this->station_y ),2 ) +
              pow(( date.position_z - this->station_z ),2 ));

    Relativity = ( -2.0 / ppp_calculate::c ) * ((date.position_x * date.velocity_x) +
                                     (date.position_y * date.velocity_y) +
                                     (date.position_z * date.velocity_z));
    Relativity1 = ( 2*ppp_calculate::u / (ppp_calculate::c * ppp_calculate::c) ) * log( (dj+di+dji) / (dj+di-dji));
    date.relativity -= Relativity + Relativity1;
}



