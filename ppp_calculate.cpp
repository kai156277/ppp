#include <algorithm>
#include <iostream>
#include <cmath>
#include <QDebug>
#include <Eigen/Eigen>

#include "ppp_calculate.h"
#include "gc_gpss.h"
#include "math_function.h"
#include "coordinate.h"

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

ppp_calculate::ppp_calculate()
    :station_x(0),station_y(0),station_z(0),
      station_B(0),station_L(0),station_H(0),
      antenna_E(0),antenna_N(0),antenna_H(0),station_ant()
{

}

void ppp_calculate::ppp_spp(const o_file &ofile,const sp3_file &sp3file,const clock_file &clockfile,const antmod_file &ant, ppp_file &ppp)
{
    //观测间隔
    double interval = clockfile.file[1].GPSS - clockfile.file[0].GPSS;
    double clock_interval = interval * 7 / 2.0;
    //坐标
    double sp3_interval = sp3file.heard.interval * 7 / 2.0;
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
        epoch.sate_number = ofile.satellite_file[i].number_of_satellite;
        epoch.GPSS = ofile.satellite_file[i].GPSS;
        GC_GPSS station;
        station.setGC(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second);
        station.GCtoDOY();
        DOY = station.DOY;
        double sunPos[3] = {0};
        sunPosition(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second,sunPos);//计算太阳位置

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
            sate_troposphere(sate,DOY); //计算对流层效应
            satellite_antenna(sate,sunPos);
            satellite_antenna_info(sate_ant,ant);
            //receiver_antenna(sate);
            epoch.sate_info.push_back(sate);
        }
        ppp.file.push_back(epoch);
    }
}

void ppp_calculate::ppp_pretreatment(const o_file &ofile, const antmod_file &ant)
{
    station_x = ofile.heard.position_X;
    station_y = ofile.heard.position_Y;
    station_z = ofile.heard.position_Z;
    Coordinate station;
    station.setXYZ(ppp_calculate::station_x,ppp_calculate::station_y,ppp_calculate::station_z).XYZtoBLH();
    station_B = fabs(station.B) * 180 / ppp_calculate::Pi;
    station_L = station.L * 180 / ppp_calculate::Pi;
    station_H = station.H;
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
}

void ppp_calculate::sate_angle(ppp_sate &date)
{
    Coordinate angle;
    angle.setENU(date.position_x, date.position_y, date.position_z ,
                 this->station_x, this->station_y, this->station_z).ENUparameter();
    date.azimuth   = angle.azimuth;
    date.elevation = angle.elevation;
    date.distance  = angle.distance;
    date.position_u =date.distance * sin(date.elevation * Pi / 180);
    double line = date.distance * cos(date.elevation * Pi / 180);
    date.position_n = line * cos(date.azimuth * Pi / 180);
    date.position_e = line * sin(date.azimuth * Pi / 180);
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
    double lat_hyd_doy=     Neil_hyd_avg[0] -     Neil_hyd_amp[0] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);;
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
    double    lat_doy = computed_avg[0] - computed_amp[0] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double      P_doy = computed_avg[1] - computed_amp[1] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double      T_doy = computed_avg[2] - computed_amp[2] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double    WVP_doy = computed_avg[3] - computed_amp[3] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double   Beta_doy = computed_avg[4] - computed_amp[4] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);
    double Lambda_doy = computed_avg[5] - computed_amp[5] * cos((doy - 28) * 2 * ppp_calculate::Pi / 365.25);

    double lat_wet_doy= Neil_wet_avg[0];
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
    double ei = date.position_e / date.distance;
    double ni = date.position_n / date.distance;
    double ui = date.position_u / date.distance;
    date.receiver_antenna_height = 0 * ei + 0 * ni +  antenna_H * ui;
    double OffsetL1 = station_ant.L1_APC_E * ei + station_ant.L1_APC_N * ni + station_ant.L1_APC_U * ui;
    double OffsetL2 = station_ant.L2_APC_E * ei + station_ant.L2_APC_N * ni + station_ant.L2_APC_U * ui;
    int e = date.elevation / (int)station_ant.DAZI;
    int a = date.azimuth / (int)station_ant.DZEN;
    if(e == 0)
    {
        e += 1;
    }
    if(a == 0)
    {
        a += 1;
    }
    double offsetL1 = station_ant.L1_NOAZI(a,  e) + station_ant.L1_NOAZI(a-1,e-1) +
               station_ant.L1_NOAZI(a-1,e) + station_ant.L1_NOAZI(a,  e-1);
    double offsetL2 = station_ant.L2_NOAZI(a,  e) + station_ant.L2_NOAZI(a-1,e-1) +
               station_ant.L2_NOAZI(a-1,e) + station_ant.L2_NOAZI(a,  e-1);
    date.offsetL1 = OffsetL1 / 1000.0 + offsetL1 / 4.0;
    date.offsetL2 = OffsetL2 / 1000.0 + offsetL2 / 4.0;


}

void ppp_calculate::satellite_antenna(ppp_sate &date,const double *posCTS)
{
    /*
    Vector3d satPos(date.position_x,date.position_y,date.position_z);
    Vector3d sunPos(posCTS[0],posCTS[1],posCTS[2]);
    Vector3d sat_sun = satPos.corss(sunPos);
    Vector3d dsat_sun = satPos.corss(sat_sun);
    double ez[3] = {0};
    double ey[3] = {0};
    double ex[3] = {0};
    double sat_distance = sqrt( pow(satPos(0),2) + pow(satPos(1),2) +pow(satPos(2),2));
    double sat_sun_distance = sqrt( pow(sat_sun(0),2) + pow(sat_sun(1),2) +pow(sat_sun(2),2));
    double dsat_sun_distance = sqrt( pow(dsat_sun(0),2) + pow(dsat_sun(1),2) + pow(dsat_sun(2),2));
    for(int i = 0; i < 3; i++)
    {
        ez[i] = satPos(i) / sat_distance;
        ey[i] = sat_sun(i) / sat_sun_distance;
        ex[i] = dsat_sun(i) / dsat_sun_distance;
    }

    Matrix3d e;
    for(int i = 0; i<3; i++)
    {
        e(i,0) = ex[i];
        e(i,1) = ey[i];
        e(i,2) = ez[i];
    }
    Vector3d PCO();
    //Vector3d PCO_ecef = e * PCO;

    double r[3] = {0};
    double r_dx = station_x - satPos(0);
    double r_dy = station_y - satPos(1);
    double r_dz = station_z - satPos(2);
    double r_distance = sqrt( pow(r_dx,2) + pow(r_dy,2) + pow(r_dz,2) );
    r[0] = r_dx / r_distance;
    r[1] = r_dy / r_distance;
    r[2] = r_dz / r_distance;

    double s =
    */
}

void ppp_calculate::sunPosition(int year, int month, int day, int hour, int minute, double second, double *posCTS)
{
    GC_GPSS sid_time;
    sid_time.setGC(year,month,day,hour,minute,second);
    sid_time.GCtoDOY();
    int doy = sid_time.DOY;
    double SOD = hour*3600+minute*60+second;
    double jd = sid_time.JD;
    double posCIS[3] ={0};
    double TWO_PI = 6.2831853071796 ;
    /* Mean Earth-Moon barycenter (EMB) distance (AU).--------------------------*/
    double MeanEarthMoonBary=3.12e-5;
    /* Astronomical Unit value (AU), in meters.---------------------------------*/
    double AU_CONST=1.49597870e11;
    /* Compute the years, and fraction of year, pased since J1900.0-------------*/
    /* Current year-------------------------------------------------------------*/
    /* Day of current year------------------------------------------------------*/
    double fd=SOD/86400.0;

    /* Fraction of day----------------------------------------------------------*/
    int years=floor(year-1900) ;
    /* Integer number of years since J1900.0------------------------------------*/
    int iy4=floor(((year%4)+4)%4);

    /* Compute fraction of year-------------------------------------------------*/
    double yearfrac=((4*(doy-1/(iy4+1))-iy4-2)+4.0*fd)/1461.0;

    double times=years+yearfrac;

    /* Compute the geometric mean longitude of the Sun--------------------------*/
    double elm1 = 4.881628+TWO_PI*yearfrac+0.0001342*times;
    double elm = elm1 - ((int)(elm1 / TWO_PI)) * TWO_PI;

    /* Mean longitude of perihelion---------------------------------------------*/
    double gamma=4.90823+0.00030005*times;

    /* Mean anomaly-------------------------------------------------------------*/
    double em=elm-gamma;

    /* Mean obliquity-----------------------------------------------------------*/
    double eps0=0.40931975-2.27e-6*times;

    /* Eccentricity-------------------------------------------------------------*/
    double e=0.016751-4.2e-7*times;
    double esq=e*e;

    /* True anomaly-------------------------------------------------------------*/
    double v=em+2.0*e*sin(em)+1.25*esq*sin(2.0*em);

    /* True ecliptic longitude--------------------------------------------------*/
    double elt=v+gamma;

    /* True distance------------------------------------------------------------*/
    double r=(1.0-esq)/(1.0+e*cos(v));

    /* Moon's mean longitude----------------------------------------------------*/
    double elmm1 = 4.72 + 83.9971*times;
    double elmm= elmm1 - ((int)(elmm1/TWO_PI)) * TWO_PI;

    /* Useful definitions-------------------------------------------------------*/
    double coselt=cos(elt);
    double sineps=sin(eps0);
    double coseps=cos(eps0);
    double w1=-r*sin(elt);
    double selmm=sin(elmm);
    double celmm=cos(elmm);

    /* Sun position is the opposite of Earth position---------------------------*/
    posCIS[0] = (r*coselt+MeanEarthMoonBary*celmm)*AU_CONST;
    posCIS[1] = (MeanEarthMoonBary*selmm-w1)*coseps*AU_CONST;
    posCIS[2] = (-w1*sineps)*AU_CONST;

    /*SID恒星时------------------------------------------------------------------*/

    /*Hours of day(decimal)-----------------------------------------------------*/
    double sid_h = SOD / 3600.0;
    double tt = (jd - 2451545.0)/36525.0;
    double sid = 24110.54841 + tt*((8640184.812866)+tt*((0.093104)-(6.2e-6*tt)));
    sid = sid / 3600.0 +sid_h;
    sid = sid - ((int)(sid / 24.0))*24.0;
    if(sid < 0)
    {
        sid = sid + 24.0;
    }

    double ts = sid * TWO_PI / 24.0;

    posCTS[0] = cos(ts)*posCIS[0]+sin(ts)*posCIS[1];
    posCTS[1] = -sin(ts)*posCIS[0]+cos(ts)*posCIS[1];
    posCTS[2] = posCIS[2];

}

void ppp_calculate::satellite_antenna_info(satellite_antmod &sate_ant, const antmod_file &ant)
{

}



