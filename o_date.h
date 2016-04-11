#ifndef O_DATE_H
#define O_DATE_H

#include<QString>
#include<QVector>
#include<QList>
#include<QStringList>


class system_record
{
public:
    system_record();
    QString system_type;
    double observation_name;
    QStringList observation_descriptor;
};

class o_heard_date
{
public:
    o_heard_date();

    /*RINEX VERSION / TYPE*/
    QString format_version;
    QString file_type;
    QString satellite_system;

    /*PGM / RUN BY / DATE */
    QString creating_program_name;
    QString creating_agency_name;
    QString creation_time;

    /*MARKER NAME*/
    QString marker_name;

    /*MARKER NUMBER*/
    QString marker_number;

    /*MARKER TYPE*/
    QString marke_type;

    /*OBSERVER / AGENCY*/
    QString observer_name;
    QString agency_name;

    /*REC # / TYPE / VERS*/
    QString receiver_number;
    QString receiver_type;
    QString receiver_version;

    /*ANT # / TYPE*/
    QString antenna_number;
    QString antenna_type;

    /*APPROX POSITION XYZ*/
    double position_X;
    double position_Y;
    double position_Z;

    /*ANTENNA: DELTA H/E/N*/
    double antenna_H;
    double antenna_E;
    double antenna_N;

    /*SYS / # / OBS TYPES*/
    QVector<system_record> system_record;

    /*INTERVAL*/
    double interval;

    /*TIME OF FIRST OBS*/
    QString first_time;
    QString time_system;

    /*SYS / PHASE SHIFT*/


    /*SIGNAL STRENGTH UNIT*/
    QString dbhz;
};

class sate_date
{
public:
    sate_date();
    QString satellite_infomation;
    /*vector中容器的数据顺序为P1,P2,P5,L1,L2,L5。缺省值为0 */
    QVector<double> satellite_observation_value;
    /*容器最大为6*/
};

class epoch_date
{
public:
    epoch_date();
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    int epoch_flag;
    int number_of_satellite;
    double clock_offset;
    QVector<sate_date> satellite_epoch;
};

class file_date
{
public:
    file_date();
    QVector<epoch_date> satellite_file;
};

/*用于选择P1,P2,P5,L1,L2,L5*/
class obs_signal
{
public:
    obs_signal();
    QList<double> P1;   //P1
    QList<double> P2;   //P2
    QList<double> P5;   //P5
    QList<double> L1;   //L1
    QList<double> L2;   //L2
    QList<double> L5;   //L5

    //适配顺序检测表
    static const char *P1_list;
    static const char *P2_list;
    static const char *P5_list;
    static const char *L1_list;
    static const char *L2_list;
    static const char *L5_list;
};

#endif // O_DATE_H


























