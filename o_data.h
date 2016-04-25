#ifndef O_DATE_H
#define O_DATE_H

#include<QString>
#include<QVector>
#include<QList>
#include<QStringList>

class phase_shift
{
public:
    phase_shift();
    QString sate_type;
    QString cpoc;       //Carrier phase observation code
    double correction;  //Correction applied
    int satellite_number;
    QStringList satellite_list;
};

class o_heard
{
public:
    o_heard();

    /*RINEX VERSION / TYPE*/
    QString format_version;
    QString file_type;

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
    //QVector<sys_record> system_record;
    int GPS_number;
    QStringList observation_descriptor;

    /*INTERVAL*/
    double interval;

    /*TIME OF FIRST OBS*/
    QString first_time;
    QString time_system;

    /*SYS / PHASE SHIFT*/
    QVector<phase_shift> phase_shift_record;

    /*SIGNAL STRENGTH UNIT*/
    QString dbhz;
};

class o_sate
{
public:
    o_sate();
    QString satellite_infomation;
    /*
     * GPS--------P1,P2,P5,L1,L2,L5
     * 缺省值为0
     */
    QVector<double> satellite_observation_value;
    QVector<int> satellite_LLI;
    QVector<int> satellite_signal_strength;
};

class o_epoch
{
public:
    o_epoch();
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    int epoch_flag;
    int number_of_satellite;
    double clock_offset;
    int GPSW;
    double GPSS;
    QVector<o_sate> satellite_epoch;
};

class o_file
{
public:
    o_file();
    QVector<o_epoch> satellite_file;
};

/*用于选择P1,P2,P5,L1,L2,L5*/
class system_signal
{
public:
    system_signal();
    /*GPS_signal        G*/
    int GPS_P1,GPS_P2,GPS_P5,
        GPS_L1,GPS_L2,GPS_L5;
    QStringList GPS_P1_list,GPS_P2_list,GPS_P5_list,
                GPS_L1_list,GPS_L2_list,GPS_L5_list;
};

#endif // O_DATE_H


























