#include <algorithm>
#include <iostream>

#include "ppp_calculate.h"
#include "gc_gpss.h"
#include "math_function.h"

ppp_calculate::ppp_calculate()
{

}

void ppp_calculate::ppp_coordinate(const o_file_date &ofile, const sp3_file &sp3file, ppp_file &ppp)
{
    double half_time_interval = sp3file.heard.interval * 7 / 2.0;
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

        //寻找符合历元的lagrange的插值起点
        double find_time = epoch.GPSS - half_time_interval;
        QVector<sp3_epoch_date>::const_iterator sp3_find_time_item = NULL;
        sp3_find_time_item = std::find_if(sp3file.file.begin(),sp3file.file.end(),sp3_GPSS_finder(find_time));
        int item = sp3_find_time_item - sp3file.file.begin(); //lagrange的插值起点
        std::cout << item << std::endl;
        //寻找每一个对应卫星的lagrange的插值
        for(int j = 0; j<ofile.satellite_file[i].GPS_satellite_epoch.size(); j++)
        {
            ppp_sate sate;
            //每一个卫星的7个插值
            double sate_x[14] = {0};    //前七个数值为GPS时间，后七个数值为坐标值
            double sate_y[14] = {0};
            double sate_z[14] = {0};
            for(int k = 0; k<7; k++)
            {
                QVector<sp3_sate_date>::const_iterator sate_find = NULL;
                sate_find = std::find(sp3file.file[item+k].GPS_epoch.begin(),
                                      sp3file.file[item+k].GPS_epoch.end(),
                                      ofile.satellite_file[i].GPS_satellite_epoch[j]);
                sate_x[k] = sp3file.file[item+k].GPSS;
                sate_y[k] = sp3file.file[item+k].GPSS;
                sate_z[k] = sp3file.file[item+k].GPSS;
                sate_x[k+7] = sate_find->x*1000;
                sate_y[k+7] = sate_find->y*1000;
                sate_z[k+7] = sate_find->z*1000;
            }
            //计算x,y,z坐标
            sate.PRN = ofile.satellite_file[i].GPS_satellite_epoch[j].satellite_infomation;
            sate.position_x = math_function::lagrange(epoch.GPSS,sate_x,14);
            sate.position_y = math_function::lagrange(epoch.GPSS,sate_y,14);
            sate.position_z = math_function::lagrange(epoch.GPSS,sate_z,14);
            epoch.sate_info.push_back(sate);
        }
        ppp.file.push_back(epoch);
    }
}



