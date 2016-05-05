#ifndef PPP_DATE_H
#define PPP_DATE_H

#include<QString>
#include<QVector>

class ppp_sate
{
public:
    ppp_sate();
    bool operator<(const ppp_sate &right) const;


    QString PRN;
    double position_x;
    double position_y;
    double position_z;
    double position_n;
    double position_e;
    double position_u;
    double velocity_x;
    double velocity_y;
    double velocity_z;
    double clock;
    double elevation;
    double azimuth;
    double distance;
    double P1;
    double P2;
    double L1;
    double L2;
    double ionized_pseudo;
    double ionized_carrier;
    double trop_delay;
    double trop_map;
    double relativity;
    double sagnac;
    double tide_effect;
    double stat_ant_height;
    double sate_antenna;
    double offsetL1;
    double offsetL2;
    double windup;
    double P3;
    double L3;

};

class ppp_epoch
{
public:
    ppp_epoch();
    int sate_number;
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    double zhd;
    double GPSS;
    double GPSW;
    QStringList PRN;
    QVector<ppp_sate> sate_info;
};

class ppp_file
{
public:
    ppp_file();
    double receiver_x;
    double receiver_y;
    double receiver_z;
    QVector<ppp_epoch> file;
};

class ppp_X
{
public:
    ppp_X();

    double dx;
    double dy;
    double dz;
    double dt;
    double trop;
    QVector<double> N;
};

class result_file
{
public:
    result_file();
    QVector<ppp_X> file;
};

#endif // PPP_DATE_H
