#ifndef ANTMOD_DATE_H
#define ANTMOD_DATE_H

#include<QString>
#include<QVector>
#include <Eigen/Eigen>

class satellite_antmod_date
{
public:
    satellite_antmod_date();
    QString antenna_type;
    QString sate_info;
    double DAZI;
    double ZEN1;
    double ZEN2;
    double DZEN;
    int num_of_frequencies;
    QString start_time;     //有效期开始的时间
    QString end_time;
    QString frequency_type;

    const QString L1_frequency;
    double L1_APC_x;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double L1_APC_y;
    double L1_APC_z;
    QVector<double> L1_NOAZI;  //PCV 天线相位中心偏差

    const QString L2_frequency;
    double L2_APC_x;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double L2_APC_y;
    double L2_APC_z;
    QVector<double> L2_NOAZI;  //PCV 天线相位中心偏差，

};

class station_antmod_date
{
public:
    station_antmod_date();
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

    const QString L1_frequency;
    double L1_APC_N;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double L1_APC_E;
    double L1_APC_U;
    Eigen::MatrixXd L1_NOAZI;  //PCV 天线相位中心偏差

    const QString L2_frequency;
    double L2_APC_N;       //x方向上卫星天线相位中心与卫星中心相对偏差
    double L2_APC_E;
    double L2_APC_U;
    Eigen::MatrixXd L2_NOAZI;  //PCV 天线相位中心偏差，
};

class antmod_file
{
public:
    antmod_file();
    QVector<satellite_antmod_date> sate;
    QVector<station_antmod_date> stat;
};

#endif // ANTMOD_DATE_H
