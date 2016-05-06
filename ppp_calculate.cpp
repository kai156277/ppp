#include <algorithm>
#include <iostream>
#include <cmath>
#include <QDebug>
#include <Eigen/Eigen>

#include "ppp_calculate.h"
#include "time_system.h"
#include "math_function.h"
#include "coordinate_system.h"
#include "sun_moon_position.h"
#include "tide_correction.h"

using namespace Eigen;
const double ppp_calculate::c = 299792458;
const double ppp_calculate::w = 7.292115e-5;
const double ppp_calculate::u = 3.986004418e+14;
const double ppp_calculate::Pi = 3.141592653;

/*UNB3参数-----------------------------------------------------------------------*/
/*------------------------------------------------纬度    大气压强    温度 水蒸气压力 气温变化率 水气压变化率*/
/*------------------------------------------------lat       P         T     WVP     Beta        Lambda-*/
const double ppp_calculate::UNB3_average[5][6] = {15.0 , 1013.25 , 299.65 , 75.00 , 6.30e-3 , 2.77,
                                                  30.0 , 1017.25 , 294.15 , 80.00 , 6.05e-3 , 3.15,
                                                  45.0 , 1015.75 , 283.15 , 76.00 , 5.58e-3 , 2.57,
                                                  60.0 , 1011.75 , 272.15 , 77.50 , 5.39e-3 , 1.81,
                                                  75.0 , 1013.00 , 263.65 , 82.50 , 4.53e-3 , 1.55};
const double ppp_calculate::UNB3_amplitude[5][6] = {15.0 ,  0.00 ,  0.00 ,  0.00 , 0.00e-3 , 0.00,
                                                    30.0 , -3.75 ,  7.00 ,  0.00 , 0.25e-3 , 0.33,
                                                    45.0 , -2.25 , 11.00 , -1.00 , 0.32e-3 , 0.46,
                                                    60.0 , -1.75 , 15.00 , -2.50 , 0.81e-3 , 0.74,
                                                    75.0 , -0.50 , 14.50 ,  2.50 , 0.62e-3 , 0.30};
/*----------------------------------------------------lat------a_hyd---------b_hyd---------c_hyd-------*/
const double ppp_calculate::Neil_hyd_average[5][4] = {15.0 ,1.2769934e-3 ,2.9153695e-3 ,62.610505e-3 ,
                                                      30.0 ,1.2683230e-3 ,2.9152299e-3 ,62.837393e-3 ,
                                                      45.0 ,1.2465397e-3 ,2.9288445e-3 ,63.721774e-3 ,
                                                      60.0 ,1.2196049e-3 ,2.9022565e-3 ,63.824265e-3 ,
                                                      75.0 ,1.2045996e-3 ,2.9024912e-3 ,64.258455e-3};
const double ppp_calculate::Neil_hyd_amplitude[5][4] = {15.0 ,         0.0 ,         0.0 ,         0.0,
                                                        30.0 ,1.2709626e-5 ,2.1414979e-5 ,9.0128400e-5,
                                                        45.0 ,2.6523662e-5 ,3.0160779e-5 ,4.3497037e-5,
                                                        60.0 ,3.4000452e-5 ,7.2562722e-5 ,84.795348e-5,
                                                        75.0 ,4.1202191e-5 ,11.723375e-5 ,170.37206e-5};
/*-------------------------------------------------------------------------------------------------------*/
const double ppp_calculate::Neil_wet_average[5][4] = {15.0 ,5.8021897e-4 ,1.4275268e-3 ,4.3472961e-2,
                                                      30.0 ,5.6794847e-4 ,1.5138625e-3 ,4.6729510e-2,
                                                      45.0 ,5.8118019e-4 ,1.4572752e-3 ,4.3908931e-2,
                                                      60.0 ,5.9727542e-4 ,1.5007428e-3 ,4.4626982e-2,
                                                      75.0 ,6.1641693e-4 ,1.7599082e-3 ,5.4736038e-2};

const double ppp_calculate::R = 287.054;
const double ppp_calculate::k1 = 77.604;
const double ppp_calculate::k2 = 16.6;
const double ppp_calculate::k3 = 377600;

/*GPS载波频率--------------------------------------------------------------------*/
const double ppp_calculate::f1 = 154 * 10.23 * 1e6;
const double ppp_calculate::f2 = 120 * 10.23 * 1e6;
const double ppp_calculate::f5 = 115 * 10.23 * 1e6;
/*GPS载波波长--------------------------------------------------------------------*/
const double ppp_calculate::lambda1 = ppp_calculate::c / ppp_calculate::f1;
const double ppp_calculate::lambda2 = ppp_calculate::c / ppp_calculate::f2;
/*天线相位缠绕预报----------------------------------------------------------------*/
double ppp_calculate::satellite_phase[40] = {0,0,0,0,0,0,0,0,0,0,
                                             0,0,0,0,0,0,0,0,0,0,
                                             0,0,0,0,0,0,0,0,0,0,
                                             0,0,0,0,0,0,0,0,0,0};
double ppp_calculate::station_phase[40] = {0,0,0,0,0,0,0,0,0,0,
                                           0,0,0,0,0,0,0,0,0,0,
                                           0,0,0,0,0,0,0,0,0,0,
                                           0,0,0,0,0,0,0,0,0,0};


ppp_calculate::ppp_calculate()
    :station_x(0),station_y(0),station_z(0),
      station_B(0),station_L(0),station_H(0),
      antenna_E(0),antenna_N(0),antenna_H(0),station_ant()
{

}

void ppp_calculate::ppp_spp(const o_file &ofile,const sp3_file &sp3file,const clock_file &clockfile,const antmod_file &ant,const erp_file &erp_data, ppp_file &ppp)
{
    ppp.receiver_x = ofile.heard.position_X;
    ppp.receiver_y = ofile.heard.position_Y;
    ppp.receiver_z = ofile.heard.position_Z;
    bool cmp(const ppp_sate &a,const ppp_sate &b);
    //观测间隔
    double interval = clockfile.file[1].GPSS - clockfile.file[0].GPSS;
    double clock_interval = interval * 7 / 2.0;
    //坐标
    double sp3_interval = sp3file.heard.interval * 7 / 2.0;
    XYZ_coordinate receiver(station_x,station_y,station_z);
    for(int i = 0; i<ofile.satellite_file.size(); i++)
    {
        ppp_epoch epoch;
        double DOY = 0; //day of year;
        epoch.year   = ofile.satellite_file[i].year;
        epoch.month  = ofile.satellite_file[i].month;
        epoch.day    = ofile.satellite_file[i].day;
        epoch.hour   = ofile.satellite_file[i].hour;
        epoch.minute = ofile.satellite_file[i].minute;
        epoch.second = ofile.satellite_file[i].second;
        epoch.GPSW = ofile.satellite_file[i].GPSW;
        epoch.GPSS = ofile.satellite_file[i].GPSS;
        if(epoch.hour == 2 && epoch.minute == 34)
        {
            int j  = 0;
        }

        GC_Time station_GC;
        station_GC.setGC_Time(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second);
        DOY_Time station_DOY = time_Transform::GCtoDOY(station_GC);
        julian_Time Julian = time_Transform::GCtoJulian(station_GC);
        GPS_Time station_GPST(epoch.GPSW,epoch.GPSS);
        DOY = station_DOY.doy;
        XYZ_coordinate sunPostion = sun_moon_position::sunPostion(Julian);
        XYZ_coordinate moonPostion = sun_moon_position::moonPostion(Julian);
        erp erp_record;


        /*--------------------------三潮改正基本值--------------------------------*/
        Vector3d solidEffect = tide_correction::solidTide(receiver,sunPostion,moonPostion,station_GPST);
        Vector3d oceanEffect = tide_correction::oceanTide(station_ocean,station_GPST);
        for(int k = 0; k <erp_data.record.size(); k++ )
        {
            if (fabs(erp_data.record[k].MJD - Julian.julian) < 1)
            {
                erp_record = erp_data.record[k];
            }
        }
        Vector3d poleEffect = tide_correction::poleTide(erp_record,receiver,station_GPST);
        Vector3d tide = poleEffect + solidEffect + oceanEffect;


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
            satellite_antmod sate_ant;

            //每一个卫星的7个插值
            double sate_clock[14] = {0};    //前七个数值为GPS时间，后七个数值为钟差值
            double sate_x[14] = {0};    //前七个数值为GPS时间，后七个数值为坐标值
            double sate_y[14] = {0};
            double sate_z[14] = {0};
            double velocity_x[14] = {0};
            double velocity_y[14] = {0};
            double velocity_z[14] = {0};
            double clock = 0;
            double clock_error = 0;
            bool find_flag = false;

            sate.PRN = ofile.satellite_file[i].satellite_epoch[j].satellite_infomation;
            sate_ant.sate_info = sate.PRN;
            sate_ant.start_year   = epoch.year;
            sate_ant.start_month  = epoch.month;
            sate_ant.start_day    = epoch.day;
            sate_ant.start_hour   = epoch.hour;
            sate_ant.start_minute = epoch.minute;
            sate_ant.start_second = epoch.second;
            sate.P1 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[0];
            sate.P2 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[1];
            sate.P3 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[2];
            sate.L1 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[3];
            sate.L2 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[4];
            sate.L3 = ofile.satellite_file[i].satellite_epoch[j].satellite_observation_value[5];

            if(sate.P1 == 0||sate.P2 == 0)
            {
                continue;
            }

            //计算钟差
            for(int k = 0; k<7; k++)
            {
                QVector<clock_info>::const_iterator sate_find = NULL;
                sate_find = std::find(clockfile.file[clock_item+k].GPS_epoch.begin(),
                        clockfile.file[clock_item+k].GPS_epoch.end(),
                        ofile.satellite_file[i].satellite_epoch[j]);
                if(sate_find != clockfile.file[clock_item+k].GPS_epoch.end())
                {
                    find_flag = true;
                    sate_clock[k] = clockfile.file[clock_item+k].GPSS;
                    sate_clock[k+7] = sate_find->clock_bias;
                }
            }
            if(find_flag == true)
            {
                clock = math_function::lagrange(epoch.GPSS,sate_clock,14);
                clock_error = clock + sate.P2 / ppp_calculate::c;
            }
            //计算坐标
            for(int k = 0; k<7; k++)
            {
                QVector<sp3_sate>::const_iterator sate_find = NULL;
                sate_find = std::find(sp3file.file[sp3_item+k].GPS_epoch.begin(),
                        sp3file.file[sp3_item+k].GPS_epoch.end(),
                        ofile.satellite_file[i].satellite_epoch[j]);
                if(sate_find != sp3file.file[sp3_item+k].GPS_epoch.end())
                {
                    find_flag = true;
                    sate_x[k] = sp3file.file[sp3_item+k].GPSS + clock_error ;
                    sate_y[k] = sp3file.file[sp3_item+k].GPSS + clock_error ;
                    sate_z[k] = sp3file.file[sp3_item+k].GPSS + clock_error ;
                    sate_x[k+7] = sate_find->x*1000;
                    sate_y[k+7] = sate_find->y*1000;
                    sate_z[k+7] = sate_find->z*1000;
                }
            }
            if(find_flag == true)
            {
                sate.position_x = math_function::lagrange(epoch.GPSS,sate_x,14);
                sate.position_y = math_function::lagrange(epoch.GPSS,sate_y,14);
                sate.position_z = math_function::lagrange(epoch.GPSS,sate_z,14);

            }

            //计算速度
            for(int k = 0; k<7; k++)
            {
                QVector<sp3_sate>::const_iterator sate_find = NULL;
                sate_find = std::find(sp3file.file[sp3_item+k].GPS_epoch.begin(),
                        sp3file.file[sp3_item+k].GPS_epoch.end(),
                        ofile.satellite_file[i].satellite_epoch[j]);
                if(sate_find != sp3file.file[sp3_item+k].GPS_epoch.end())
                {
                    velocity_x[k] = sp3file.file[sp3_item+k].GPSS + clock_error - 0.5;
                    velocity_y[k] = sp3file.file[sp3_item+k].GPSS + clock_error - 0.5;
                    velocity_z[k] = sp3file.file[sp3_item+k].GPSS + clock_error - 0.5;
                    velocity_x[k+7] = sate_find->x*1000;
                    velocity_y[k+7] = sate_find->y*1000;
                    velocity_z[k+7] = sate_find->z*1000;
                }

            }


            if(find_flag == true)
            {
                sate.velocity_x = (sate.position_x - math_function::lagrange(epoch.GPSS,velocity_x,14)) / 0.5;
                sate.velocity_y = (sate.position_y - math_function::lagrange(epoch.GPSS,velocity_y,14)) / 0.5;
                sate.velocity_z = (sate.position_z - math_function::lagrange(epoch.GPSS,velocity_z,14)) / 0.5;
                sate.clock = clock * ppp_calculate::c;


                sate_angle(sate);       //计算卫星角度信息
                if(sate.elevation < 15)
                {
                    continue;
                }
                epoch.PRN.push_back(sate.PRN);
                sate_sagnac(sate);      //计算地球自转效应
                sate_relativity(sate);  //计算相对论效应
                sate_troposphere(sate,DOY); //计算对流层效应
                receiver_antenna(sate);
                satellite_antenna_info(sate_ant,ant);
                satellite_antenna(sate,sate_ant,sunPostion);
                satellite_windup(sate,sate_ant,sunPostion);
                satellite_tide(sate,tide);
                satellite_ionized(sate);
                epoch.sate_info.push_back(sate);
            }
        }
        epoch.sate_number = epoch.sate_info.size();
        if(epoch.sate_number >=6 )
        {
            std::sort(epoch.PRN.begin(),epoch.PRN.end());
            std::sort(epoch.sate_info.begin(),epoch.sate_info.end(),cmp);
            ppp.file.push_back(epoch);
        }

    }
}

bool cmp(const ppp_sate &a, const ppp_sate &b)
{
    return a < b;
}

void ppp_calculate::ppp_pretreatment(const o_file &ofile, const antmod_file &ant, const ocean_file &ocean_data)
{
    /*---------------------------测站的坐标信息-----------------------------------*/
    station_x = ofile.heard.position_X;
    station_y = ofile.heard.position_Y;
    station_z = ofile.heard.position_Z;
    XYZ_coordinate XYZ_station;
    XYZ_station.setXYZ(ppp_calculate::station_x,ppp_calculate::station_y,ppp_calculate::station_z);
    BLH_coordinate BLH_station = coordinate_transform::XYZtoBLH(XYZ_station);
    station_B = BLH_station.get_latitude_angle();
    station_L = BLH_station.get_longitude_angle();
    station_H = BLH_station.getHeight();



    /*----------------------------测站的天线信息----------------------------------*/
    antenna_E = ofile.heard.antenna_E;
    antenna_N = ofile.heard.antenna_N;
    antenna_H = ofile.heard.antenna_H;
    bool ant_flag = false;
    for(int i = 0; i<ant.station.size(); i++)
    {
        if(ant.station[i].antenna_type == ofile.heard.antenna_type)
        {
            station_ant = ant.station[i];
            ant_flag = true;
            break;
        }
    }
    if(ant_flag == false)
    {
        station_ant.L1_NOAZI = Eigen::MatrixXd::Zero(2,2);
        station_ant.L2_NOAZI = Eigen::MatrixXd::Zero(2,2);
    }

    /*----------------------------测站的海潮信息----------------------------------*/
    for(int i = 0; i<ocean_data.record.size(); i++)
    {
        if(ofile.heard.marker_name == ocean_data.record[i].stationName)
        {
            station_ocean = ocean_data.record[i];
        }
    }

}

void ppp_calculate::sate_angle(ppp_sate &date)
{
    XYZ_coordinate station_XYZ(this->station_x, this->station_y, this->station_z);
    XYZ_coordinate satellite_XYZ(date.position_x, date.position_y, date.position_z );
    ENU_coordiante satellite_ENU = coordinate_transform::sateENU(station_XYZ,satellite_XYZ);
    date.azimuth   = satellite_ENU.getAzimuth();
    date.elevation = satellite_ENU.getElevation();
    date.distance  = satellite_ENU.getDistance();
    date.position_e = satellite_ENU.getSate_e();
    date.position_n = satellite_ENU.getSate_n();
    date.position_u = satellite_ENU.getSate_u();
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

void ppp_calculate::sate_troposphere(ppp_sate &date, int doy)
{
    double computed_avg[6] = {0};    //平均值的计算值
    double computed_amp[6] = {0};    //振幅的计算值
    double Neil_wet_avg[4] = {0};

    /*计算Avg,Amp,---------------------------------------------------------------*/
    if(station_B <= 15)
    {
        for(int i = 0; i<6; i++)
        {
            computed_avg[i] = UNB3_average[0][i];
            computed_amp[i] = UNB3_amplitude[0][i];
        }
    }
    else if(station_B >= 75)
    {
        for(int i = 0; i<6; i++)
        {
            computed_avg[i] = UNB3_average[4][i];
            computed_amp[i] = UNB3_amplitude[4][i];
        }
    }
    else
    {
        int index = ((int)station_B / 15)-1;
        for(int i = 0; i<6; i++)
        {
            computed_avg[i] = UNB3_average[index][i] +
                    (UNB3_average[index+1][i] - UNB3_average[index][i]) / 15.0 * (station_B - UNB3_average[index][0]);

            computed_amp[i] = UNB3_amplitude[index][i] +
                    (UNB3_amplitude[index+1][i] - UNB3_amplitude[index][i]) / 15.0 * (station_B - UNB3_amplitude[index][0]);
        }
    }

    /*计算干延迟投影函数Neil_hyd_Avg,Neil_hyd_Amp---------------------------------*/
    double Neil_hyd_avg[4] = {0};
    double Neil_hyd_amp[4] = {0};
    if(station_B <= 15)
    {
        for(int i = 0; i<4; i++)
        {
            Neil_hyd_avg[i] = Neil_hyd_average[0][i];
            Neil_hyd_amp[i] = Neil_hyd_amplitude[0][i];
        }
    }
    else if(station_B >= 75)
    {
        for(int i = 0; i<4; i++)
        {
            Neil_hyd_avg[i] = Neil_hyd_average[4][i];
            Neil_hyd_amp[i] = Neil_hyd_amplitude[4][i];
        }
    }
    else
    {
        int index = ((int)station_B / 15)-1;
        for(int i = 0; i<4; i++)
        {
            Neil_hyd_avg[i] = Neil_hyd_average[index][i] +
                    (Neil_hyd_average[index+1][i] - Neil_hyd_average[index][i]) / 15.0 * (station_B - Neil_hyd_average[index][0]);

            Neil_hyd_amp[i] = Neil_hyd_amplitude[index][i] +
                    (Neil_hyd_amplitude[index+1][i] - Neil_hyd_amplitude[index][i]) / 15.0 * (station_B - Neil_hyd_amplitude[index][0]);
        }
    }
    double  a_hyd_doy =     Neil_hyd_avg[1] -     Neil_hyd_amp[1] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);;
    double  b_hyd_doy =     Neil_hyd_avg[2] -     Neil_hyd_amp[2] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);;
    double  c_hyd_doy =     Neil_hyd_avg[3] -     Neil_hyd_amp[3] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);;

    /*-----------------------------------------------------------------------------*/

    /*计算湿延迟投影函数Neil_wet_Avg----------------------------------------------*/
    if(station_B <= 15)
    {
        for(int i = 0; i<4; i++)
        {
            Neil_wet_avg[i] = Neil_wet_average[0][i];
        }
    }
    else if(station_B >= 75)
    {
        for(int i = 0; i<4; i++)
        {
            Neil_wet_avg[i] = Neil_wet_average[4][i];
        }
    }
    else
    {
        int index = ((int)station_B / 15)-1;
        for(int i = 0; i<4; i++)
        {
            Neil_wet_avg[i] = Neil_wet_average[index][i] +
                    (Neil_wet_average[index+1][i] - Neil_wet_average[index][i]) / 15.0 * (station_B - Neil_wet_average[index][0]);
        }
    }

    /*当天的气象条件--------------------------------------------------------------*/
    double      P_doy = computed_avg[1] - computed_amp[1] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double      T_doy = computed_avg[2] - computed_amp[2] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double    WVP_doy = computed_avg[3] - computed_amp[3] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double   Beta_doy = computed_avg[4] - computed_amp[4] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double Lambda_doy = computed_avg[5] - computed_amp[5] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);

    double  a_wet_doy = Neil_wet_avg[1];
    double  b_wet_doy = Neil_wet_avg[2];
    double  c_wet_doy = Neil_wet_avg[3];

    /*计算湿延迟分量---------------------------------------------------------------*/

    double B = station_B * Pi / 180;
    double gm = 9.784 * (1 - 2.66e-3 * cos(2 * B) - 2.8e-7 * ppp_calculate::station_H);
    double g = 9.78049* (1 + 5.29e-3 * pow(sin(B),2));
    double T = T_doy - Beta_doy * ppp_calculate::station_H;
    double Tm = T * (1 - (Beta_doy * ppp_calculate::R) / (gm * (Lambda_doy + 1)));

    double dzh1 = (1e-6 * ppp_calculate::k1 * ppp_calculate::R) / gm;
    double dzh2 = P_doy;
    double dzh3 = 1 - (Beta_doy * ppp_calculate::station_H) / T_doy;
    double dzh4 = g / (ppp_calculate::R * Beta_doy);
    double dzh = dzh1 * dzh2 * pow(dzh3,dzh4);
    /*计算干延迟分量--------------------------------------------------------------*/

    double dzw1 = 1e-6 * (Tm * ppp_calculate::k2 + ppp_calculate::k3) * ppp_calculate::R;
    double dzw2 = gm * (Lambda_doy + 1) - Beta_doy * ppp_calculate::R;
    double dzw3 = WVP_doy / T_doy;
    double dzw4 = dzh3;
    double dzw5 = (((Lambda_doy + 1) * g) / (ppp_calculate::R * Beta_doy)) - 1;
    double dzw = (dzw1 / dzw2) * dzw3 * pow(dzw4,dzw5);

    /*Neil 湿延迟投影函数---------------------------------------------------------*/
    double E = sin(date.elevation * Pi / 180);
    double m_wet1 = 1 + (a_wet_doy / (1 + b_wet_doy / (1 + c_wet_doy)));
    double m_wet2 = E + (a_wet_doy / (E + b_wet_doy / (E + c_wet_doy)));
    double m_wet = m_wet1 / m_wet2;

    /*Neil 干延迟投影函数---------------------------------------------------------*/
    double a_hgt = 2.53e-5;
    double b_hgt = 5.49e-3;
    double c_hgt = 1.14e-3;

    double m_hyd1 = 1 + (a_hyd_doy / (1 + b_hyd_doy / (1 + c_hyd_doy)));
    double m_hyd2 = E + (a_hyd_doy / (E + b_hyd_doy / (E + c_hyd_doy)));
    double m_hyd3 = 1 + (a_hgt / (1 + b_hgt / (1 + c_hgt)));
    double m_hyd4 = E + (a_hgt / (E + b_hgt / (E + c_hgt)));
    double m_hyd = (m_hyd1 / m_hyd2) + (1/E - m_hyd3 / m_hyd4) * (ppp_calculate::station_H / 1000);
    /*----------------------------------------------------------------------------*/



    date.trop_map = m_wet;
    date.trop_delay = dzh * m_hyd;
}

void ppp_calculate::receiver_antenna(ppp_sate &date)
{
    /*当计算接收机的表格时使用天顶角------------------------------------------------*/
    double ei = date.position_e / date.distance;
    double ni = date.position_n / date.distance;
    double ui = date.position_u / date.distance;
    double station_zenith = 90 - date.elevation;
    date.stat_ant_height = antenna_E * ei + antenna_N * ni +  antenna_H * ui;
    double PCO1 = station_ant.L1_APC_E * ei + station_ant.L1_APC_N * ni + station_ant.L1_APC_U * ui;
    double PCO2 = station_ant.L2_APC_E * ei + station_ant.L2_APC_N * ni + station_ant.L2_APC_U * ui;
    int e = station_zenith / (int)station_ant.DAZI;
    int a = date.azimuth / (int)station_ant.DZEN;
    if(e == 18)
    {
        e -= 1;
    }
    if(a == 72)
    {
        a -= 1;
    }
    double offsetL1 = station_ant.L1_NOAZI(a,  e) + station_ant.L1_NOAZI(a+1,e+1) +
               station_ant.L1_NOAZI(a+1,e) + station_ant.L1_NOAZI(a,  e+1);
    double offsetL2 = station_ant.L2_NOAZI(a,  e) + station_ant.L2_NOAZI(a+1,e+1) +
               station_ant.L2_NOAZI(a+1,e) + station_ant.L2_NOAZI(a,  e+1);
    date.offsetL1 = (PCO1  - offsetL1 / 4.0) / 1000.0 / ppp_calculate::lambda1;
    date.offsetL2 = (PCO2  - offsetL2 / 4.0) / 1000.0 / ppp_calculate::lambda2;
}

void ppp_calculate::satellite_antenna(ppp_sate &date, const satellite_antmod &sate_ant, const XYZ_coordinate &sunPostion)
{
    Vector3d satPos(date.position_x,date.position_y,date.position_z);
    Vector3d sunPos(sunPostion.getX(),sunPostion.getY(),sunPostion.getZ());
    Vector3d sat_sun = satPos.cross(sunPos);
    Vector3d dsat_sun = satPos.cross(sat_sun);
    Vector3d ez(3);
    Vector3d ey(3);
    Vector3d ex(3);

    /*1.计算星固系在惯性系中的单位向量 [ez,ey,ex]-----------------------------------*/
    ez = -1.0 * satPos / satPos.norm();
    ey = -1.0 * sat_sun / sat_sun.norm();
    ex = -1.0 * dsat_sun / dsat_sun.norm();

    /*2.计算卫星到接收机的单位向量 r-----------------------------------------------*/
    Vector3d PCO(sate_ant.APC_x,sate_ant.APC_y,sate_ant.APC_z);
    Vector3d P_rec_sat(station_x - satPos(0),station_y - satPos(1),station_z - satPos(2));
    Vector3d r = P_rec_sat / P_rec_sat.norm();

    /*3.计算ez与r之间的夹角-------------------------------------------------------*/
    double s = ez.dot(r);
    (s < -1.0) ? (s = -1.0) : 0 ;
    (s > 1.0)  ? (s =  1.0) : 0 ;
    double angle = acos(s) * 180 / Pi;
    if( angle < sate_ant.ZEN1 || angle >sate_ant.ZEN2 )
    {
        qDebug() << "This error is in the satellite antenna function! "
                 << "satellite " << date.PRN  << " angle is out of ZEN1~ZEN2."
                 << "satellite enu coordinate is:" <<endl
                 << "elevation :" << date.elevation
                 << "azimuth :" << date.azimuth << endl;
        exit(1);
    }
    if( sate_ant.DZEN == 0)
    {
        qDebug() << "This error is in the satellite antenna function! "
                 << "satellite " << date.PRN  << " type is wrong!."
                 << "satellite enu coordinate is:" <<endl
                 << "elevation :" << date.elevation
                 << "azimuth :" << date.azimuth << endl;
        exit(1);
    }
    /*4.根据夹角内插PCV*/
    int c = ceil(angle) ;
    int f = floor(angle);
    double PCV_U = (sate_ant.NOAZI[c] + sate_ant.NOAZI[f]) / 2;
    Vector3d PCV(0,0,PCV_U);
    /*5.将星固系下的PCC转换到地固系下的PCC------------------------------------------*/
    Vector3d PCC = PCO - PCV;
    Vector3d PCC_ECEF = ex * PCC(0) + ey * PCC(1) + ez * PCC(2);

    date.sate_antenna = r.dot(PCC_ECEF) / 1000;
}

void ppp_calculate::satellite_windup(ppp_sate &date, const satellite_antmod &sate_ant, const XYZ_coordinate &sunPostion)
{
    /*--------------------------------星固系-------------------------------------
     * 1. y轴平行与太阳能帆板，垂直与太阳矢量
     * 2. z轴平行与天线指向，始终指向地球
     * 3. x轴由y轴和z轴的矢量积确定指向
     * -------------------------------------------------------------------------*/
    Vector3d sun_position(sunPostion.getX(),sunPostion.getY(),sunPostion.getZ());
    Vector3d sate_position(date.position_x,date.position_y,date.position_z);
    Vector3d stat_position(station_x,station_y,station_z);
    double station_latitude = station_B * Pi / 180.0;
    double station_longitude = station_L * Pi / 180.0;


    /*----------------------获得当前卫星的星固系坐标-------------------------------*/

    /*1.获得在地心地固系下卫星质量中心到太阳质量中心的单位向量 sate_sun_i--------------*/
    Vector3d sate_sun_i = sun_position - sate_position;
    sate_sun_i = sate_sun_i / sate_sun_i.norm();

    /*2.获得在地心地固系下卫星质量中心到地球质量中心的单位向量 sate_earth_i------------*/
    Vector3d sate_z_i = -1 * sate_position / sate_position.norm();

    /*3.获得星固系坐标的y轴，即sate_sun_i和sate_earth_i矢量积的单位向量 sate_y_i-----*/
    Vector3d sate_y_i = sate_z_i.cross(sate_sun_i);
    sate_y_i = sate_y_i / sate_y_i.norm();

    /*4. 获得星固系坐标的x轴------------------------------------------------------*/
    Vector3d sate_x_i = sate_y_i.cross(sate_z_i);
    sate_x_i = sate_x_i / sate_x_i.norm();



    /*-------------------------计算卫星的绕转角-----------------------------------*/

    /*1.计算卫星到测站的单位矢量 sate_stat_i---------------------------------------*/
    Vector3d sate_stat_i = stat_position - sate_position;
    sate_stat_i = sate_stat_i / sate_stat_i.norm();

    /*2.将sate_z_i投影到sate_stat_i方向上-----------------------------------------*/
    double zk = sate_stat_i.dot(sate_z_i);

    /*3.获得一个没有sate_z_i分量的向量，即位于sate_x_i和sate_y_i的平面上-------------*/
    Vector3d dpp = sate_stat_i - zk * sate_z_i;

    /*4.将dpp投影到sate_x_i和sate_y_i方向上---------------------------------------*/
    double xk = dpp.dot(sate_x_i);
    double yk = dpp.dot(sate_y_i);

    /*5.计算卫星的绕转角(弧度)----------------------------------------------------*/
    double alpha1 = atan2(yk,xk);



    /*--------------------------获得接收机的旋转角--------------------------------*/

    /*1.计算接收机到地心的单位向量-------------------------------------------------*/
    Vector3d stat_earth_i = -1.0 * stat_position / stat_position.norm();

    /*2.在UEN坐标系中定义一个仅包含北方向的向量-------------------------------------*/
    Vector3d delta(0,0,1);

    /*3.定义ENU到XYZ的旋转矩阵----------------------------------------------------*/
    Matrix3d rot = Matrix3d::Zero();
    rot << cos(station_latitude),0,-sin(station_latitude),
                               0,1,0,
           sin(station_latitude),0,cos(station_latitude);
    delta = rot * delta;
    rot << cos(-station_longitude),sin(-station_longitude),0,
          -sin(-station_longitude),cos(-station_longitude),0,
           0,0,1;
    delta = rot * delta;
    delta = delta.transpose();

    /*4.计算测站参考系y轴的单位向量------------------------------------------------*/
    Vector3d stat_y_i = stat_earth_i.cross(delta);
    stat_y_i = stat_y_i / stat_y_i.norm();

    /*5.计算测站参考系x轴的单位向量------------------------------------------------*/
    Vector3d stat_x_i = stat_y_i.cross(stat_earth_i);
    stat_x_i = stat_x_i / stat_x_i.norm();

    /*6.将将stat_earth_i投影到sate_stat_i方向上-----------------------------------*/
    double rzk = sate_stat_i.dot(stat_earth_i);

    /*7.获得一个没有stat_earth_i分量的向量，即位于stat_x_i和stat_y_i的平面上---------*/
    Vector3d rdpp = sate_stat_i - rzk * stat_earth_i;

    /*8.将rdpp投影到stat_x_i和stat_y_i方向上--------------------------------------*/
    double rxk = rdpp.dot(stat_x_i);
    double ryk = rdpp.dot(stat_y_i);

    /*9.计算测站的绕转角(弧度)----------------------------------------------------*/
    double alpha2 = atan2(ryk,rxk);



    /*--------------------------计算天线相位缠绕的影响-----------------------------*/

    /*1.查找卫星天线类型是否含有"IIR",这种类型的天线有180的相位-----------------------*/
    double wind_up = 0.0;
    if(sate_ant.antenna_type.indexOf("IIR")>=0)
    {
        wind_up = Pi;
    }
    alpha1 = alpha1 + wind_up;

    int PRN = date.PRN.mid(1,2).toInt();
    double phase_satellite = satellite_phase[PRN];
    double da1 = alpha1 - phase_satellite;
    double phase_station   = station_phase[PRN];
    double da2 = alpha2 - phase_station;
    phase_satellite = phase_satellite + atan2(sin(da1),cos(da1));
    phase_station = phase_station + atan2(sin(da2),cos(da2));

    wind_up = phase_satellite - phase_station;
    date.windup = wind_up / (2 *Pi);

    satellite_phase[PRN] = phase_satellite;
    station_phase[PRN] = phase_station;
}

void ppp_calculate::satellite_tide(ppp_sate &data, const Vector3d &tide)
{
    double ele = data.elevation / 180 * Pi;
    double azi = data.azimuth / 180 *Pi;
    double cosel = cos(ele);
    double sinel = sin(ele);
    double cosaz = cos(azi);
    double sinaz = sin(azi);
    data.tide_effect = tide[0]*cosel*cosaz + tide[1]*cosel*sinaz + tide[2]*sinel;
}



int ppp_calculate::satellite_antenna_info(satellite_antmod &sate_ant, const antmod_file &ant)
{
    julian_Time JD;
    GC_Time time;
    time.setGC_Time(sate_ant.start_year  ,sate_ant.start_month,
               sate_ant.start_day   ,sate_ant.start_hour,
               sate_ant.start_minute,sate_ant.start_second);
    JD = time_Transform::GCtoJulian(time);
    double sate_ant_JD = JD.julian;


    QVector<satellite_antmod>::const_iterator sate_ant_find = ant.satellite.begin();
    for(int i = 0; i<ant.satellite.size(); i++)
    {
        double start_JD = 0;
        double end_JD = 0;
        sate_ant_find = std::find(sate_ant_find,ant.satellite.end(),sate_ant);
        if(sate_ant_find != ant.satellite.end())
        {

            time.setGC_Time(sate_ant.start_year  ,sate_ant.start_month,
                       sate_ant.start_day   ,sate_ant.start_hour,
                       sate_ant.start_minute,sate_ant.start_second);
            JD = time_Transform::GCtoJulian(time);
            start_JD = JD.julian;

            if(sate_ant_find->end_time !=  "")
            {
                time.setGC_Time(sate_ant_find->end_year  ,sate_ant_find->end_month,
                           sate_ant_find->end_day   ,sate_ant_find->end_hour,
                           sate_ant_find->end_minute,sate_ant_find->end_second);
                JD = time_Transform::GCtoJulian(time);
                end_JD = JD.julian;
                if(sate_ant_JD>=start_JD && sate_ant_JD <= end_JD)
                {
                    sate_ant = ant.satellite.at(sate_ant_find-ant.satellite.begin());
                    return 1;
                }
            }
            else
            {
                if(sate_ant_JD>=start_JD )
                {
                    sate_ant = ant.satellite.at(sate_ant_find-ant.satellite.begin());
                    return 1;
                }
            }
            sate_ant_find += 1;
        }
        else
        {
            return 0;
        }
    }
    return 0;
}

void ppp_calculate::satellite_ionized(ppp_sate &date)
{
    double pseudo1 = date.P1 + lambda1*date.offsetL1;
    double pseudo2 = date.P2 + lambda2*date.offsetL2;
    double carrier1  = lambda1 * (date.L1+date.offsetL1-date.windup);
    double carrier2  = lambda2 * (date.L2+date.offsetL2-date.windup);
    double line1 =  pow(f1,2) / (pow(f1,2)-pow(f2,2));
    double line2 = -pow(f2,2) / (pow(f1,2)-pow(f2,2));
    date.ionized_pseudo = line1 * pseudo1 + line2 * pseudo2;
    date.ionized_carrier = line1 * carrier1 + line2 * carrier2;
}
