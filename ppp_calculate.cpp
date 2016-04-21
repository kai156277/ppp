#include <algorithm>
#include <iostream>

#include "ppp_calculate.h"
#include "gc_gpss.h"
#include "math_function.h"

double ppp_calculate::c = 299792458;

ppp_calculate::ppp_calculate()
{

}

void ppp_calculate::ppp_coordinate(const o_file_date &ofile, const sp3_file &sp3file, ppp_file &ppp)
{
    double sp3_interval = sp3file.heard.interval * 7 / 2.0;
    for(int i = 0; i<ofile.satellite_file.size(); i++)
    {
        ppp_epoch epoch;
        //寻找符合历元的坐标lagrange的插值起点
        double sp3_find_time = epoch.GPSS - sp3_interval;
        QVector<sp3_epoch>::const_iterator sp3_find_time_item = NULL;
        sp3_find_time_item = std::find_if(sp3file.file.begin(),sp3file.file.end(),sp3_GPSS_finder(sp3_find_time));
        int sp3_item = sp3_find_time_item - sp3file.file.begin(); //lagrange的插值起点
        //寻找每一个对应卫星的lagrange的插值
        for(int j = 0; j<ofile.satellite_file[i].satellite_epoch.size(); j++)
        {
            ppp_sate sate;
            //每一个卫星的7个插值
            double sate_x[14] = {0};    //前七个数值为GPS时间，后七个数值为坐标值
            double sate_y[14] = {0};
            double sate_z[14] = {0};
            for(int k = 0; k<7; k++)
            {
                QVector<sp3_sate>::const_iterator sate_find = NULL;
                sate_find = std::find(sp3file.file[sp3_item+k].GPS_epoch.begin(),
                                      sp3file.file[sp3_item+k].GPS_epoch.end(),
                                      ofile.satellite_file[i].satellite_epoch[j]);
                sate_x[k] = sp3file.file[sp3_item+k].GPSS;
                sate_y[k] = sp3file.file[sp3_item+k].GPSS;
                sate_z[k] = sp3file.file[sp3_item+k].GPSS;
                sate_x[k+7] = sate_find->x*1000;
                sate_y[k+7] = sate_find->y*1000;
                sate_z[k+7] = sate_find->z*1000;
            }
            //计算x,y,z坐标
            sate.PRN = ofile.satellite_file[i].satellite_epoch[j].satellite_infomation;
            sate.position_x = math_function::lagrange(epoch.GPSS,sate_x,14);
            sate.position_y = math_function::lagrange(epoch.GPSS,sate_y,14);
            sate.position_z = math_function::lagrange(epoch.GPSS,sate_z,14);
            epoch.sate_info.push_back(sate);
        }
        ppp.file.push_back(epoch);
    }
}

void ppp_calculate::ppp_clock(const o_file_date &ofile,const sp3_file &sp3file,const clock_file &clockfile, ppp_file &ppp)
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
        std::cout << clock_item << std::endl;

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
            sate.PRN = ofile.satellite_file[i].satellite_epoch[j].satellite_infomation;
            sate.clock = clock * ppp_calculate::c;
            epoch.sate_info.push_back(sate);
        }
        ppp.file.push_back(epoch);
    }
}



