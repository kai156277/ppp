#include<QString>
#include<QVector>
#include<QStringList>
#include<QFile>
#include<QTextStream>
#include<qdebug.h>
#include<cmath>
#include<Eigen/Eigen>

#include"file_read.h"
#include"o_data.h"
#include"time_system.h"

using namespace std;
using namespace Eigen;

file_read::file_read()
{

}

void file_read::ppp_o_read(const QString &file_path, o_file &obs)
{
    QFile ppp_o_file( file_path );
    if(!ppp_o_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit(EXIT_FAILURE);
    }
    QTextStream read( &ppp_o_file );
    QString readString;

    /*read_heard*/

    do
    {
        readString = read.readLine();
        if(readString.indexOf("RINEX VERSION / TYPE") >= 0)
        {
            obs.heard.format_version        = readString.mid(0,9).trimmed();
            obs.heard.file_type             = readString.mid(40,1);
        }
        else if(readString.indexOf("PGM / RUN BY / DATE") >= 0)
        {
            obs.heard.creating_program_name = readString.mid(0,20).trimmed();
            obs.heard.creating_agency_name  = readString.mid(20,20).trimmed();
            obs.heard.creation_time         = readString.mid(40,20).trimmed();
        }
        else if(readString.indexOf("MARKER NAME") >= 0)
        {
            obs.heard.marker_name = readString.mid(0,60).trimmed();
        }
        else if(readString.indexOf("MARKER NUMBER") >= 0)
        {
            obs.heard.marker_number = readString.mid(0,60).trimmed();
        }
        else if(readString.indexOf("MARKER TYPE") >= 0)
        {
            obs.heard.marke_type = readString.mid(0,60).trimmed();
        }
        else if(readString.indexOf("OBSERVER / AGENCY") >= 0)
        {
            obs.heard.observer_name = readString.mid(0,20);
            obs.heard.agency_name   = readString.mid(20,40).trimmed();
        }
        else if(readString.indexOf("REC # / TYPE / VERS") >= 0)
        {
            obs.heard.receiver_number  = readString.mid(0,20).trimmed();
            obs.heard.receiver_type    = readString.mid(20,20).trimmed();
            obs.heard.receiver_version = readString.mid(40,20).trimmed();
        }
        else if(readString.indexOf("ANT # / TYPE") >= 0)
        {
            obs.heard.antenna_number = readString.mid(0,20).trimmed();
            obs.heard.antenna_type   = readString.mid(20,20).trimmed();
        }
        else if(readString.indexOf("APPROX POSITION XYZ") >= 0)
        {
            obs.heard.position_X = readString.mid(0,14).toDouble();
            obs.heard.position_Y = readString.mid(14,14).toDouble();
            obs.heard.position_Z = readString.mid(28,14).toDouble();
        }
        else if(readString.indexOf("ANTENNA: DELTA H/E/N") >= 0)
        {
            obs.heard.antenna_H = readString.mid(0,14).toDouble();
            obs.heard.antenna_E = readString.mid(0,14).toDouble();
            obs.heard.antenna_N = readString.mid(0,14).toDouble();
        }
        else if(readString.indexOf("SYS / # / OBS TYPES") >= 0)
        {
            if(readString.mid(0,1) == "G")
            {
                obs.heard.GPS_number = readString.mid(1,5).toInt();
                int i = ceil(obs.heard.GPS_number / 13.0); //有几行同样卫星的数据
                QString sate_info;
                for(; i>0; i--)
                {
                    sate_info += readString.mid(7,53);
                    if( i!=1 )
                    {
                        readString = read.readLine();
                    }
                }
                sate_info = sate_info.trimmed();
                obs.heard.observation_descriptor = sate_info.split(' ');
            }

        }
        else if(readString.indexOf("INTERVAL") >= 0)
        {
            obs.heard.interval = readString.mid(0,10).toDouble();
        }
        else if(readString.indexOf("TIME OF FIRST OBS") >= 0)
        {
            obs.heard.first_time = readString.mid(0,43);
            obs.heard.time_system = readString.mid(43,8).trimmed();
        }
        else if(readString.indexOf("SIGNAL STRENGTH UNIT") >= 0)
        {
            obs.heard.dbhz = readString.mid(0,20).trimmed();
        }
    }while(readString.indexOf("END OF HEADER")<=0);

    /*phase matching*/
    system_signal sys_signal;
    phase_matching(obs.heard.observation_descriptor ,sys_signal);

    /*date*/

    while(!read.atEnd())
    {
        o_epoch epoch;
        readString = read.readLine();
        epoch.year   = readString.mid(2,4).toInt();
        epoch.month  = readString.mid(6,3).toInt();
        epoch.day    = readString.mid(9,3).toInt();
        epoch.hour   = readString.mid(12,3).toInt();
        epoch.minute = readString.mid(15,3).toInt();
        epoch.second = readString.mid(18,11).toDouble();
        epoch.epoch_flag = readString.mid(29,3).toInt();
        epoch.number_of_satellite = readString.mid(32,3).toInt();
        epoch.clock_offset = readString.mid(41,15).toDouble();

        GC_Time GC;
        GC.setGC_Time(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second);
        GPS_Time GPST = time_Transform::GCtoGPST(GC);
        epoch.GPSW = GPST.GPSW;
        epoch.GPSS = GPST.GPSS;

        for(int i = 0; i<epoch.number_of_satellite; i++)
        {
            readString = read.readLine();
            o_sate satellite_data;
            satellite_data.satellite_infomation = readString.mid(0,3);

            if(satellite_data.satellite_infomation.mid(0,1) == "G")
            {
                QVector<QString> sate;
                sate.reserve(obs.heard.GPS_number);
                for(int j = 0;j<obs.heard.GPS_number ;j++)
                {
                    sate.push_back(readString.mid(3+j*16,16));
                }
                //添加一个用于指示没有匹配的观测信号的状况
                sate.push_back("0000000000000000");

                satellite_data.satellite_observation_value[0] = sate[sys_signal.GPS_P1].mid(0,14).toDouble();
                satellite_data.satellite_LLI[0] = sate[sys_signal.GPS_P1].mid(14,1).toInt();
                satellite_data.satellite_signal_strength[0] = sate[sys_signal.GPS_P1].mid(15,1).toInt();

                satellite_data.satellite_observation_value[1] = sate[sys_signal.GPS_P2].mid(0,14).toDouble();
                satellite_data.satellite_LLI[1] = sate[sys_signal.GPS_P2].mid(14,1).toInt();
                satellite_data.satellite_signal_strength[1] = sate[sys_signal.GPS_P2].mid(15,1).toInt();

                satellite_data.satellite_observation_value[2] = sate[sys_signal.GPS_P5].mid(0,14).toDouble();
                satellite_data.satellite_LLI[2] = sate[sys_signal.GPS_P5].mid(14,1).toInt();
                satellite_data.satellite_signal_strength[2] = sate[sys_signal.GPS_P5].mid(15,1).toInt();

                satellite_data.satellite_observation_value[3] = sate[sys_signal.GPS_L1].mid(0,14).toDouble();
                satellite_data.satellite_LLI[3] = sate[sys_signal.GPS_L1].mid(14,1).toInt();
                satellite_data.satellite_signal_strength[3] = sate[sys_signal.GPS_L1].mid(15,1).toInt();

                satellite_data.satellite_observation_value[4] = sate[sys_signal.GPS_L2].mid(0,14).toDouble();
                satellite_data.satellite_LLI[4] = sate[sys_signal.GPS_L2].mid(14,1).toInt();
                satellite_data.satellite_signal_strength[4] = sate[sys_signal.GPS_L2].mid(15,1).toInt();

                satellite_data.satellite_observation_value[5] = sate[sys_signal.GPS_L5].mid(0,14).toDouble();
                satellite_data.satellite_LLI[5] = sate[sys_signal.GPS_L5].mid(14,1).toInt();
                satellite_data.satellite_signal_strength[5] = sate[sys_signal.GPS_L5].mid(15,1).toInt();

                epoch.satellite_epoch.push_back(satellite_data);
            }
        }
        epoch.number_of_satellite = epoch.satellite_epoch.size();
        obs.satellite_file.push_back(epoch);
    }
    ppp_o_file.close();
}

void file_read::ppp_sp3_read(const QString &file_path, sp3_file &sp3)
{
    QFile ppp_sp3_file( file_path );
    if(!ppp_sp3_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit( EXIT_FAILURE );
    }
    QTextStream read( &ppp_sp3_file );
    QString readString;

    /*read_heard*/

    /*#C*/
    readString = read.readLine();
    sp3.heard.mode_flag  = readString.mid(2,1);
    sp3.heard.year       = readString.mid(3,4).toInt();
    sp3.heard.month      = readString.mid(7,3).toInt();
    sp3.heard.day        = readString.mid(10,3).toInt();
    sp3.heard.hour       = readString.mid(13,3).toInt();
    sp3.heard.minute     = readString.mid(16,3).toInt();
    sp3.heard.second     = readString.mid(19,12).toDouble();
    sp3.heard.number_of_epoch = readString.mid(31,8).toInt();
    sp3.heard.data_use   = readString.mid(39,6).trimmed();
    sp3.heard.coordinate_system = readString.mid(45,6).trimmed();
    sp3.heard.orbit_type = readString.mid(51,4).trimmed();
    sp3.heard.agency     = readString.mid(55,5).trimmed();

    /*##*/
    readString = read.readLine();
    sp3.heard.GPS_week_number  = readString.mid(3,4).toInt();
    sp3.heard.second_of_week   = readString.mid(7,16).toDouble();
    sp3.heard.interval         = readString.mid(23,15).toDouble();
    sp3.heard.modified_JD      = readString.mid(38,6).toInt();
    sp3.heard.fractional_day   = readString.mid(44,16).toDouble();

    /*+*/
    sp3.heard.satellites.reserve(85);
    for(int i = 0; i<5; i++)
    {
        readString = read.readLine();
        sp3.heard.number_of_satellites += readString.mid(4,2).toInt();
        for(int j = 0; j<17; j++)
        {
            sp3.heard.satellites.push_back(readString.mid(9+j*3,3));
        }
    }
    for(int k = 85; k>sp3.heard.number_of_satellites; k--)
    {
        sp3.heard.satellites.pop_back();
    }
    QStringList num ;


    num = sp3.heard.satellites.filter("C");
    sp3.heard.before_GPS += num.size();
    num = sp3.heard.satellites.filter("E");
    sp3.heard.before_GPS += num.size();
    num = sp3.heard.satellites.filter("G");
    sp3.heard.GPS_satellites = num.size();
    num = sp3.heard.satellites.filter("J");
    sp3.heard.after_GPS += num.size();
    num = sp3.heard.satellites.filter("R");
    sp3.heard.after_GPS += num.size();
    num = sp3.heard.satellites.filter("S");
    sp3.heard.after_GPS += num.size();
    /*++*/
    sp3.heard.accuracy.reserve(85);
    for(int i=0; i<5; i++)
    {
        readString = read.readLine();
        for(int j = 0; j<17; j++)
        {
            sp3.heard.accuracy.push_back(readString.mid(9+j*3,3).toInt());
        }
    }
    for(int k = 85; k>sp3.heard.number_of_satellites; k--)
    {
        sp3.heard.accuracy.pop_back();
    }

    /*%c*/
    readString = read.readLine();
    sp3.heard.file_type = readString.mid(3,1);
    sp3.heard.time_type = readString.mid(9,3);
    readString = read.readLine();

    /*%f*/
    readString = read.readLine();
    sp3.heard.p_v_base = readString.mid(3 ,10).toDouble();
    sp3.heard.c_c_base = readString.mid(14,12).toDouble();
    for(int i = 0; i<7; i++)
    {
        read.readLine();
    }
    /*date*/

    for(int i = 0; i<sp3.heard.number_of_epoch; i++)
    {
        sp3_epoch epoch;
        readString = read.readLine();
        epoch.year   = readString.mid(3,4).toInt();
        epoch.month  = readString.mid(7,3).toInt();
        epoch.day    = readString.mid(10,3).toInt();
        epoch.hour   = readString.mid(13,3).toInt();
        epoch.minute = readString.mid(16,3).toInt();
        epoch.second = readString.mid(19,12).toDouble();

        GC_Time GC;
        GC.setGC_Time(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second);
        GPS_Time GPST = time_Transform::GCtoGPST(GC);
        epoch.GPSW = GPST.GPSW;
        epoch.GPSS = GPST.GPSS;

        //读取GPS卫星之前的信息
        for(int j = 0; j<sp3.heard.before_GPS; j++)
        {
            readString = read.readLine();
        }
        //G
        for(int j = 0; j<sp3.heard.GPS_satellites; j++)
        {
            sp3_sate sate_data;
            readString = read.readLine();
            sate_data.flag  = readString.mid(0,1);
            sate_data.sate_info = readString.mid(1,3);
            sate_data.x     = readString.mid(4,14).toDouble();
            sate_data.y     = readString.mid(18,14).toDouble();
            sate_data.z     = readString.mid(32,14).toDouble();
            sate_data.clock = readString.mid(46,14).toDouble();
            sate_data.x_SD  = readString.mid(60,3).toInt();
            sate_data.y_SD  = readString.mid(63,3).toInt();
            sate_data.z_SD  = readString.mid(66,3).toInt();
            sate_data.clock_SD = readString.mid(69,4).toInt();
            epoch.GPS_epoch.push_back(sate_data);
        }
        //读取GPS卫星之后的信息
        for(int j = 0; j<sp3.heard.after_GPS; j++)
        {
            readString = read.readLine();
        }
        sp3.file.push_back(epoch);
    }

    ppp_sp3_file.close();
}

void file_read::ppp_clock_read(const QString &file_path, clock_file &clock)
{
    QFile ppp_clock_file( file_path );
    if(!ppp_clock_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit( EXIT_FAILURE );
    }
    QTextStream read( &ppp_clock_file );
    QString readString;

    /*read heard*/

    do
    {
        readString = read.readLine();
        if(readString.indexOf("RINEX VERSION / TYPE")>=0)
        {
            clock.heard.rinex_format_version = readString.mid(0,9);
            clock.heard.rinex_file_type      = readString.mid(20,1);
            clock.heard.satellite_system     = readString.mid(40,20).trimmed();
        }
        else if(readString.indexOf("PGM / RUN BY / DATE")>=0)
        {
            clock.heard.creating_program_name = readString.mid(0,20).trimmed();
            clock.heard.creating_agency_name  = readString.mid(20,20).trimmed();
            clock.heard.creation_time         = readString.mid(40,20).trimmed();
        }
        else if(readString.indexOf("# / TYPES OF DATA")>=0)
        {
            clock.heard.clock_type_num = readString.mid(0,6).toInt();
            for(int i = 0; i<clock.heard.clock_type_num; i++)
            {
                clock.heard.clock_type.push_back(readString.mid(6+i*6+6).trimmed());
            }
        }
        else if(readString.indexOf("# OF SOLN STA / TRF")>=0)
        {
            clock.heard.reference_num   = readString.mid(0,6).toInt();
            clock.heard.satation_info.reserve(clock.heard.reference_num);
            clock.heard.reference_frame = readString.mid(10,50).trimmed();
        }
        else if(readString.indexOf("SOLN STA NAME / NUM")>=0)
        {
            station_infomation info;
            info.station_name = readString.mid(0,5);
            info.station_id   = readString.mid(5,20);
            info.station_x    = readString.mid(25,12).toDouble();
            info.station_y    = readString.mid(37,12).toDouble();
            info.station_z    = readString.mid(49,11).toDouble();
            clock.heard.satation_info.push_back(info);
        }
        else if(readString.indexOf("# OF SOLN SATS")>=0)
        {
            clock.heard.satellite_num = readString.mid(0,6).toInt();
        }
        else if(readString.indexOf("PRN LIST")>=0)
        {
            int num = ceil(clock.heard.satellite_num / 15.0);
            for(int i = 0; i<num; i++)
            {
                for(int j = 0; j<15; j++)
                {
                    clock.heard.satellite_list.push_back(readString.mid(0+j*4,3));
                }
                readString = read.readLine();
            }
        }
    }while(readString.indexOf("END OF HEADER")<=0);

    /*获取每种卫星数量*/
    QStringList num ;
    num = clock.heard.satellite_list.filter("C");
    clock.heard.BDS_satellites = num.size();
    num = clock.heard.satellite_list.filter("G");
    clock.heard.GPS_satellites = num.size();
    num = clock.heard.satellite_list.filter("R");
    clock.heard.GLONASS_satellites = num.size();
    num = clock.heard.satellite_list.filter("J");
    clock.heard.QZSS_satellites = num.size();
    num = clock.heard.satellite_list.filter("E");
    clock.heard.Galileo_satellites = num.size();
    num = clock.heard.satellite_list.filter("S");
    clock.heard.SBAS_satellites = num.size();

    /*read data*/
    bool readFlag = true;
    do
    {
        if(readFlag == true)
        {
            readString = read.readLine();
        }
        QString clock_type = readString.mid(0,2);
        if(clock_type == "AS")
        {
            readFlag = false;
            clock_epoch epoch;
            epoch.year       = readString.mid(8,4).toInt();
            epoch.month      = readString.mid(12,3).toInt();
            epoch.day        = readString.mid(15,3).toInt();
            epoch.hour       = readString.mid(18,3).toInt();
            epoch.minute     = readString.mid(21,3).toInt();
            epoch.second     = readString.mid(24,10).toDouble();
            epoch.number_of_data = readString.mid(34,6).toInt();

            GC_Time GC;
            GC.setGC_Time(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second);
            GPS_Time GPST = time_Transform::GCtoGPST(GC);
            epoch.GPSW = GPST.GPSW;
            epoch.GPSS = GPST.GPSS;

            for(int i = 0; i<6; i++)
            {
                if(readString.mid(3,1) == "G")
                {
                    for(int j = 0; j<clock.heard.GPS_satellites; j++)
                    {

                        clock_info sate_date;
                        readString = read.readLine();
                        sate_date.sate_name = readString.mid(3,3);
                        sate_date.clock_bias = readString.mid(40,20).toDouble();
                        sate_date.clock_bias_sigma = readString.mid(60,20).toDouble();
                        epoch.GPS_epoch.push_back(sate_date);
                    }
                }
                else if(readString.mid(3,1) == "R")
                {
                    for(int j = 0; j<clock.heard.GLONASS_satellites; j++)
                    {
                        readString = read.readLine();
                    }
                }
                else if(readString.mid(3,1) == "E")
                {
                    for(int k = 0; k<clock.heard.Galileo_satellites; k++)
                    {
                        readString = read.readLine();
                    }
                }
                else if(readString.mid(3,1) == "C")
                {
                    for(int j = 0; j<clock.heard.BDS_satellites; j++)
                    {
                        readString = read.readLine();
                    }
                }
                else if(readString.mid(3,1) == "S")
                {
                    for(int k = 0; k<clock.heard.SBAS_satellites; k++)
                    {
                        readString = read.readLine();
                    }
                }
                else if(readString.mid(3,1) == "J")
                {
                    for(int k = 0; k<clock.heard.QZSS_satellites; k++)
                    {
                        readString = read.readLine();
                    }
                }
            }
            clock.file.push_back(epoch);
        }
        else
        {
            readFlag = true;
        }
    }while(!read.atEnd());
    ppp_clock_file.close();
}

void file_read::ppp_snx_read(const QString &file_path, snx_data &snx, QString mark_name)
{
    QFile ppp_snx_file( file_path );
    if(!ppp_snx_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit( EXIT_FAILURE );
    }
    QTextStream read( &ppp_snx_file );
    QString readString;

    do
    {
        readString = read.readLine();
    }while(readString.indexOf("+SOLUTION/ESTIMATE")<0);

    do
    {
        readString = read.readLine();
        if((readString.indexOf(mark_name)>=0)&&(readString.indexOf("STAX")>=0))
        {
            snx.station_x = readString.mid(47,21).toDouble();
        }
        else if((readString.indexOf(mark_name)>=0)&&(readString.indexOf("STAY")>=0))
        {
            snx.station_y = readString.mid(47,21).toDouble();
        }
        else if((readString.indexOf(mark_name)>=0)&&(readString.indexOf("STAZ")>=0))
        {
            snx.station_z = readString.mid(47,21).toDouble();
        }
    }while(readString.indexOf("-SOLUTION/ESTIMATE")<0);
    ppp_snx_file.close();
}

void file_read::ppp_ant_read(const QString &file_path, antmod_file &ant)
{
    /*Open antmond file and Create input stream---------------------------------*/
    QFile ppp_ant_file( file_path );
    if(!ppp_ant_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit( EXIT_FAILURE );
    }
    QTextStream read( &ppp_ant_file );
    QString readString;

    /*read header of antmod file------------------------------------------------*/
    while(read.readLine().indexOf("END OF HEADER")<0)
    {}

    /*read antmond file data----------------------------------------------------*/
    do
    {
        readString = read.readLine();

        if(readString.indexOf("START OF ANTENNA")>=0)
        {
    /*use DAZI value to judge antenna type--------------------------------------*/
            QString antenna_type,name;
            QString DAZI = "";
            do
            {
                readString = read.readLine();
                if(readString.indexOf("TYPE / SERIAL NO")>=0)
                {
                    antenna_type = readString.mid(0,20).trimmed();
                    name = readString.mid(20,20).trimmed();
                }
                else if(readString.indexOf("DAZI")>=0)
                {
                    DAZI = readString.mid(0,8).trimmed();
                }
            }while(readString.indexOf("DAZI")<0);
    /*if antenna type is satellite----------------------------------------------*/
            if(DAZI == "0.0")
            {
                satellite_antmod sate;
                sate.antenna_type = antenna_type;
                sate.sate_info = name;
                sate.DAZI = DAZI.toDouble();
                do
                {
                    readString = read.readLine();
                    if(readString.indexOf("ZEN1 / ZEN2 / DZEN")>=0)
                    {
                        sate.ZEN1 = readString.mid(0,8).toDouble();
                        sate.ZEN2 = readString.mid(8,6).toDouble();
                        sate.DZEN = readString.mid(14,6).toDouble();
                    }
                    else if(readString.indexOf("# OF FREQUENCIES")>=0)
                    {
                        sate.num_of_frequencies = readString.mid(0,6).toInt();
                    }
                    else if(readString.indexOf("VALID FROM")>=0)
                    {
                        sate.start_time   = readString.mid(2,41);
                        sate.start_year   = readString.mid(2,4).toInt();
                        sate.start_month  = readString.mid(6,6).toInt();
                        sate.start_day    = readString.mid(12,6).toInt();
                        sate.start_hour   = readString.mid(18,6).toInt();
                        sate.start_minute = readString.mid(24,6).toInt();
                        sate.start_second = readString.mid(30,12).toDouble();

                    }
                    else if(readString.indexOf("VALID UNTIL")>=0)
                    {
                        sate.end_time   = readString.mid(2,41);
                        sate.end_year   = readString.mid(2,4).toInt();
                        sate.end_month  = readString.mid(6,6).toInt();
                        sate.end_day    = readString.mid(12,6).toInt();
                        sate.end_hour   = readString.mid(18,6).toInt();
                        sate.end_minute = readString.mid(24,6).toInt();
                        sate.end_second = readString.mid(30,12).toDouble();
                    }
                    else if(readString.indexOf("SINEX CODE")>=0)
                    {
                        sate.frequency_type = readString.mid(0,60).trimmed();
                    }
                    else if(readString.indexOf("START OF FREQUENCY")>=0)
                    {
    /*L1 frequency--------------------------------------------------------------*/
                        if(readString.indexOf("G01")>=0)
                        {
                            readString = read.readLine();
                            sate.APC_x = readString.mid(0,10).toDouble();
                            sate.APC_y = readString.mid(10,10).toDouble();
                            sate.APC_z = readString.mid(20,10).toDouble();
                            readString = read.readLine().simplified();
                            QStringList NOAZI = readString.split(" ");
                            for(int i = 1; i<NOAZI.size(); i++)
                            {
                                sate.NOAZI.push_back(NOAZI[i].toDouble());
                            }
                        }
    /*L2 frequency--------------------------------------------------------------*/
                    }
                }while(readString.indexOf("END OF ANTENNA")<0);
                ant.satellite.push_back(sate);
            }
    /*if antenna type is station------------------------------------------------*/
            else if(DAZI == "5.0")
            {
                station_antmod stat;
                stat.antenna_type = antenna_type;
                stat.receiver_info = name;
                stat.DAZI = DAZI.toDouble();
                do
                {
                    readString = read.readLine();
                    if(readString.indexOf("ZEN1 / ZEN2 / DZEN")>=0)
                    {
                        stat.ZEN1 = readString.mid(0,8).toDouble();
                        stat.ZEN2 = readString.mid(8,6).toDouble();
                        stat.DZEN = readString.mid(14,6).toDouble();
                    }
                    else if(readString.indexOf("# OF FREQUENCIES")>=0)
                    {
                        stat.num_of_frequencies = readString.mid(0,6).toInt();
                    }
                    else if(readString.indexOf("VALID FROM")>=0)
                    {
                        stat.start_time = readString.mid(2,41);
                    }
                    else if(readString.indexOf("VALID UNTIL")>=0)
                    {
                        stat.end_time = readString.mid(2,41);
                    }
                    else if(readString.indexOf("SINEX CODE")>=0)
                    {
                        stat.frequency_type = readString.mid(0,60).trimmed();
                    }
                    else if(readString.indexOf("START OF FREQUENCY")>=0)
                    {
    /*L1 frequency--------------------------------------------------------------*/
                        int line = 73;
                        int list = stat.ZEN2 / stat.DZEN + 1 ;
                        if(readString.indexOf("G01")>=0)
                        {
                            readString = read.readLine();
                            stat.L1_APC_N = readString.mid(0,10).toDouble();
                            stat.L1_APC_E = readString.mid(10,10).toDouble();
                            stat.L1_APC_U = readString.mid(20,10).toDouble();
                            read.readLine();
                            stat.L1_NOAZI =  MatrixXd::Random(line,list);
                            for(int i = 0; i<73; i++)
                            {
                                readString = read.readLine().simplified();
                                QStringList NOAZI = readString.split(" ");
                                for(int j = 1; j<NOAZI.size(); j++)
                                {
                                    stat.L1_NOAZI(i,j-1) = NOAZI[j].toDouble();
                                }
                            }
                        }
    /*L2 frequency--------------------------------------------------------------*/
                        else if(readString.indexOf("G02")>=0)
                        {
                            readString = read.readLine();
                            stat.L2_APC_N = readString.mid(0,10).toDouble();
                            stat.L2_APC_E = readString.mid(10,10).toDouble();
                            stat.L2_APC_U = readString.mid(20,10).toDouble();
                            read.readLine();
                            stat.L2_NOAZI = MatrixXd::Random(line,list);
                            for(int i = 0; i<73; i++)
                            {
                                readString = read.readLine().simplified();
                                QStringList NOAZI = readString.split(" ");
                                for(int j = 1; j<NOAZI.size(); j++)
                                {
                                    stat.L2_NOAZI(i,j-1) = NOAZI[j].toDouble();
                                }
                            }
                        }
                    }
                }while(readString.indexOf("END OF ANTENNA")<0);
                ant.station.push_back(stat);
            }
        }
    }while(!read.atEnd());
    ppp_ant_file.close();
}

void file_read::ppp_ocean_read(const QString &file_path, ocean_file &ocean_data)
{
    /*Open ocean file and Create input stream---------------------------------*/
    QFile ppp_ocean_file( file_path );
    if(!ppp_ocean_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit( EXIT_FAILURE );
    }
    QTextStream read( &ppp_ocean_file );
    QString readString;

    /*read ocean file heard-----------------------------------------------------*/
    do
    {
        readString = read.readLine();
        if(readString.indexOf("COLUMN ORDER")>=0)
        {
            ocean_data.signer = readString.mid(16).simplified().split(" ");
            ocean_data.col = ocean_data.signer.size();
        }
    }while(readString.indexOf("END HEADER")<0);

    /*read ocean file data------------------------------------------------------*/
    QStringList month;
    month.reserve(12);
    month.push_back("JAN");month.push_back("FEB");month.push_back("MAR");
    month.push_back("APR");month.push_back("MAY");month.push_back("JUN");
    month.push_back("Jul");month.push_back("AUG");month.push_back("SEP");
    month.push_back("OCT");month.push_back("NOV");month.push_back("DEC");

    while(!read.atEnd())
    {
        readString = read.readLine();
        ocean temp;
        readString = read.readLine();
        temp.stationName = readString.mid(2,4);
        readString = read.readLine();
        temp.oceanModel = readString.mid(12,6);
        readString = read.readLine();
        temp.year = readString.mid(63,4).toInt();
        temp.month = month.indexOf(readString.mid(68,3)) + 1;
        temp.day = readString.mid(72,2).toInt();
        readString = read.readLine();
        temp.station_B = readString.mid(49,10).toDouble();
        temp.station_L = readString.mid(59,10).toDouble();
        temp.station_H = readString.mid(69,10).toDouble();
        for (int i=0;i<6;i++)
        {
            readString = read.readLine();
            for(int j=0;j<ocean_data.col;j++)
            {
                temp.data[i][j]=(readString.mid(1+7*j,7)).toDouble();
            }
        }
        ocean_data.record.push_back( temp );
    }
    ppp_ocean_file.close();
}

void file_read::ppp_erp_read(const QString &file_path, erp_file &erp_data)
{
    /*Open erp file and Create input stream---------------------------------*/
    QFile ppp_erp_file( file_path );
    if(!ppp_erp_file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Can`t open" << file_path << endl;
        exit( EXIT_FAILURE );
    }
    QTextStream read( &ppp_erp_file );
    QString readString = read.readLine();
    const double jdToMjd = 2400000.5;
    bool flag1 = false;
    bool flag2 = false;
    while(1)
    {
        readString = read.readLine();
        if(readString.indexOf("dpsi")>=0)
        {
            flag1 = true;
        }
        if(readString.indexOf("deps")>=0)
        {
            flag2 = true;
        }
        if(readString.indexOf("MJD")>=0)
        {
            break;
        }
    }
    readString = read.readLine();
    while(!read.atEnd())
    {
        readString = read.readLine();
        erp temp;
        temp.MJD       = readString.mid( 0, 8).toDouble() + jdToMjd;
        temp.Xpole     = readString.mid( 8, 9).toDouble()*1e-6;
        temp.Ypole     = readString.mid(17, 9).toDouble()*1e-6;
        temp.UT1_UTC   = readString.mid(26,10).toDouble();
        if(flag1 == true)
        {
            temp.dpsi = readString.mid(112,7).toDouble();
        }
        if(flag2 == true)
        {
            temp.deps = readString.mid(119,8).toDouble();
        }
        erp_data.record.push_back( temp );
    }

    ppp_erp_file.close();
}

void file_read::phase_matching(const QStringList &match_list, system_signal &sys_list)
{

    sys_list.GPS_P1 = each_phase_matching(sys_list.GPS_P1_list,match_list);
    sys_list.GPS_P2 = each_phase_matching(sys_list.GPS_P2_list,match_list);
    sys_list.GPS_P5 = each_phase_matching(sys_list.GPS_P5_list,match_list);
    sys_list.GPS_L1 = each_phase_matching(sys_list.GPS_L1_list,match_list);
    sys_list.GPS_L2 = each_phase_matching(sys_list.GPS_L2_list,match_list);
    sys_list.GPS_L5 = each_phase_matching(sys_list.GPS_L5_list,match_list);
}

int file_read::each_phase_matching(const QStringList &reference_list, const QStringList &obs_descriptor)
{
    int i = reference_list.size();
    for(;i>0;i--)//zuijinque
    {
        int j = obs_descriptor.indexOf(reference_list[i-1]);
        if( j != -1 )
        {
            return j;
        }
    }
    return obs_descriptor.size();; //当没有观测信号时我们将它指向提供的序列的末尾
}

