#ifndef ANTMOD_DATE_H
#define ANTMOD_DATE_H

#include<QString>
#include<QVector>
#include <Eigen/Eigen>

class satellite_antmod
{
public:
    satellite_antmod();
    bool operator ==(const satellite_antmod &right) const;
    QString antenna_type;
    QString sate_info;
    double DAZI;
    double ZEN1;
    double ZEN2;
    double DZEN;
    int num_of_frequencies;
    QString start_time;     //有效期开始的时间
    int start_year;
    int start_month;
    int start_day;
    int start_hour;
    int start_minute;
    double start_second;
    QString end_time;
    int end_year;
    int end_month;
    int end_day;
    int end_hour;
    int end_minute;
    double end_second;
    QString frequency_type;

    QString frequency;
    double APC_x;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double APC_y;
    double APC_z;
    QVector<double> NOAZI;  //PCV 天线相位中心偏差

};

class station_antmod
{
public:
    station_antmod();
    QString antenna_type;
    QString receiver_info;
    double DAZI;
    double ZEN1;
    double ZEN2;
    double DZEN;
    int num_of_frequencies; //频率的个数
    QString start_time;     //有效期开始的时间
    QString end_time;
    QString frequency_type;

    QString L1_frequency;
    double L1_APC_N;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double L1_APC_E;
    double L1_APC_U;
    Eigen::MatrixXd L1_NOAZI;  //PCV 天线相位中心偏差

    QString L2_frequency;
    double L2_APC_N;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double L2_APC_E;
    double L2_APC_U;
    Eigen::MatrixXd L2_NOAZI;  //PCV 天线相位中心偏差，
};

class antmod_file
{
public:
    antmod_file();
    QVector<satellite_antmod> satellite;
    QVector<station_antmod> station;
};

#endif // ANTMOD_DATE_H
