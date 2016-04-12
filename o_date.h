#ifndef O_DATE_H
#define O_DATE_H

#include<QString>
#include<QVector>
#include<QList>
#include<QStringList>


class sys_record
{
public:
    sys_record();
    QString system_type;
    int observation_number;
    QStringList observation_descriptor;
};

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

class o_heard_date
{
public:
    o_heard_date();

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
    QVector<sys_record> system_record;

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

class o_sate_date
{
public:
    o_sate_date();
    QString satellite_infomation;
    /*
     * GPS--------P1,P2,P5,L1,L2,L5
     * GLONASS----G1,G2,G3,L1,L2,L3
     * Galileo----E1,E5a,E5b,E5,E6,L1,L5a,L5b,L5,L6
     * SBAS-------P1,P5,L1,L5
     * QZSS-------P1,P2,P5,P6,L1,L2,L5,L6
     * BDS--------B1,B2,B3,L1,L2,L3
     * 缺省值为0
     *
     */
    QVector<double> satellite_observation_value;
    QVector<int> satellite_LLI;
    QVector<int> satellite_signal_strength;
};

class o_epoch_date
{
public:
    o_epoch_date();
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    int epoch_flag;
    int number_of_satellite;
    double clock_offset;
    QVector<o_sate_date> satellite_epoch;
};

class o_file_date
{
public:
    o_file_date();
    QVector<o_epoch_date> satellite_file;
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
    /*GLONASS_signal    R*/
    int GLONASS_G1,GLONASS_G2,GLONASS_G3,
        GLONASS_L1,GLONASS_L2,GLONASS_L3;
    QStringList GLONASS_G1_list,GLONASS_G2_list,GLONASS_G3_list,
                GLONASS_L1_list,GLONASS_L2_list,GLONASS_L3_list;
    /*Galileo_signal    E*/
    int Galileo_E1,Galileo_E5a,Galileo_E5b,Galileo_E5,Galileo_E6,
        Galileo_L1,Galileo_L5a,Galileo_L5b,Galileo_L5,Galileo_L6;
    QStringList Galileo_E1_list,Galileo_E5a_list,Galileo_E5b_list,Galileo_E5_list,Galileo_E6_list,
                Galileo_L1_list,Galileo_L5a_list,Galileo_L5b_list,Galileo_L5_list,Galileo_L6_list;
    /*SBAS_signal       S*/
    int SBAS_P1,SBAS_P5,
        SBAS_L1,SBAS_L5;
    QStringList SBAS_P1_list,SBAS_P5_list,
                SBAS_L1_list,SBAS_L5_list;
    /*QZSS_signal       J*/
    int QZSS_P1,QZSS_P2,QZSS_P5,QZSS_P6,
        QZSS_L1,QZSS_L2,QZSS_L5,QZSS_L6;
    QStringList QZSS_P1_list,QZSS_P2_list,QZSS_P5_list,QZSS_P6_list,
                QZSS_L1_list,QZSS_L2_list,QZSS_L5_list,QZSS_L6_list;
    /*BDS_signal        B*/
    int BDS_B1,BDS_B2,BDS_B3,
        BDS_L1,BDS_L2,BDS_L3;
    QStringList BDS_B1_list,BDS_B2_list,BDS_B3_list,
                BDS_L1_list,BDS_L2_list,BDS_L3_list;
};

#endif // O_DATE_H


























