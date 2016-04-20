#ifndef CLOCK_DATE_H
#define CLOCK_DATE_H

#include<QString>
#include<QStringList>
#include<QVector>

class station_infomation
{
public:
    station_infomation();
    QString station_name;
    QString station_id;
    double station_x;
    double station_y;
    double station_z;
};

class clock_date_heard
{
public:
    clock_date_heard();

    /*RINEX VERSION / TYPE*/
    QString rinex_format_version;
    QString rinex_file_type;
    QString satellite_system;

    /*PGM / RUN BY / DATE*/
    QString creating_program_name;
    QString creating_agency_name;
    QString creation_time;

    /*# / TYPES OF DATA*/
    int clock_type_num;
    QStringList clock_type;

    /**/

    /*# OF SOLN STA / TRF*/
    int reference_num;
    QString reference_frame;

    /*SOLN STA NAME / NUM*/
    QVector<station_infomation> satation_info;

    /*# OF SOLN SATS*/
    int satellite_num;

    /*PRN LIST*/
    QStringList satellite_list;
    int BDS_satellites;
    int Galileo_satellites;
    int GPS_satellites;
    int QZSS_satellites;
    int GLONASS_satellites;
    int SBAS_satellites;
};

class clock_info
{
public:
    clock_info();
    QString sate_name;           //Receiver or Satellite Name
    QVector<double> record;
};

class clock_epoch
{
public:
    clock_epoch();
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    double GPSS;
    int GPSW;
    int number_of_data;
    QVector<clock_info> GPS_epoch;
    QVector<clock_info> BDS_epoch;
    QVector<clock_info> GLONASS_epoch;
};

class clock_file
{
public:
    clock_file();
    clock_date_heard heard;
    QVector<clock_epoch> file;
};

#endif // CLOCK_DATE_H
