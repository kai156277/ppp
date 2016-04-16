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
    int satation_num;
    QString satation_reference_frame;

    /*SOLN STA NAME / NUM*/
    QVector<station_infomation> satation_info;

    /*# OF SOLN SATS*/
    int satellite_num;

    /*PRN LIST*/
    QStringList satellite_list;
};

class clock_info
{
public:
    clock_info();
    QString clock_type;
    QString R_S_name;           //Receiver or Satellite Name
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    int number_of_data;
    QVector<double> info;
};

class clock_date
{
public:
    clock_date();
    QVector<clock_info> file;
};

#endif // CLOCK_DATE_H
