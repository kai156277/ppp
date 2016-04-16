#include<QString>
#include<QVector>
#include<QStringList>
#include<QFile>
#include<QTextStream>
#include<qdebug.h>
#include<cmath>
#include<ctype.h>

#include "file_read.h"
#include "o_date.h"

using namespace std;

file_read::file_read()
{

}

void file_read::ppp_o_read(const QString &file_path, o_file_date &obs)
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
    o_heard_date heard_date;
    do
    {
        readString = read.readLine();
        if(readString.indexOf("RINEX VERSION / TYPE") >= 0)
        {
            heard_date.format_version        = readString.mid(0,9).simplified();
            heard_date.file_type             = readString.mid(40,1);
        }
        else if(readString.indexOf("PGM / RUN BY / DATE") >= 0)
        {
            heard_date.creating_program_name = readString.mid(0,20).simplified();
            heard_date.creating_agency_name  = readString.mid(20,20).simplified();
            heard_date.creation_time         = readString.mid(40,20).simplified();
        }
        else if(readString.indexOf("MARKER NAME") >= 0)
        {
            heard_date.marker_name = readString.mid(0,60).simplified();
        }
        else if(readString.indexOf("MARKER NUMBER") >= 0)
        {
            heard_date.marker_number = readString.mid(0,60).simplified();
        }
        else if(readString.indexOf("MARKER TYPE") >= 0)
        {
            heard_date.marke_type = readString.mid(0,60).simplified();
        }
        else if(readString.indexOf("OBSERVER / AGENCY") >= 0)
        {
            heard_date.observer_name = readString.mid(0,20);
            heard_date.agency_name   = readString.mid(20,40).simplified();
        }
        else if(readString.indexOf("REC # / TYPE / VERS") >= 0)
        {
            heard_date.receiver_number  = readString.mid(0,20).simplified();
            heard_date.receiver_type    = readString.mid(20,20).simplified();
            heard_date.receiver_version = readString.mid(40,20).simplified();
        }
        else if(readString.indexOf("ANT # / TYPE") >= 0)
        {
            heard_date.antenna_number = readString.mid(0,20).simplified();
            heard_date.antenna_type   = readString.mid(20,20).simplified();
        }
        else if(readString.indexOf("APPROX POSITION XYZ") >= 0)
        {
            heard_date.position_X = readString.mid(0,14).toDouble();
            heard_date.position_Y = readString.mid(14,14).toDouble();
            heard_date.position_Z = readString.mid(28,14).toDouble();
        }
        else if(readString.indexOf("ANTENNA: DELTA H/E/N") >= 0)
        {
            heard_date.antenna_H = readString.mid(0,14).toDouble();
            heard_date.antenna_E = readString.mid(0,14).toDouble();
            heard_date.antenna_N = readString.mid(0,14).toDouble();
        }
        else if(readString.indexOf("SYS / # / OBS TYPES") >= 0)
        {
            sys_record sate_record;
            sate_record.system_type = readString.mid(0,1);
            sate_record.observation_number = readString.mid(1,5).toInt();
            if(sate_record.observation_number <= 13)
            {
                for(int i = 0; i<sate_record.observation_number; i++)
                {
                    sate_record.observation_descriptor.push_back(readString.mid(7+i*4,3));
                }
            }
            else
            {
                int i = ceil(sate_record.observation_number / 13.0); //有几行同样卫星的数据
                for(int j = 0; j<i; j++)
                {
                    for(int k = 0; k<13; k++)
                    {
                        sate_record.observation_descriptor.push_back(readString.mid(7+k*4,3));
                    }
                    if(j < i-1)
                    {
                        readString = read.readLine();
                    }
                }
                int j = sate_record.observation_descriptor.size()-sate_record.observation_number;
                for(int k = 0; k < j; k++)
                {
                    sate_record.observation_descriptor.pop_back();
                }
            }
            heard_date.system_record.push_back(sate_record);
        }
        else if(readString.indexOf("INTERVAL") >= 0)
        {
            heard_date.interval = readString.mid(0,10).toDouble();
        }
        else if(readString.indexOf("TIME OF FIRST OBS") >= 0)
        {
            heard_date.first_time = readString.mid(0,43);
            heard_date.time_system = readString.mid(43,8).simplified();
        }
        else if(readString.indexOf("SYS / PHASE SHIFT") >= 0)
        {
            phase_shift record;
            record.sate_type = readString.mid(0,1);
            record.cpoc = readString.mid(2,3);
            record.correction = readString.mid(6,8).toDouble();
            record.satellite_number = readString.mid(16,2).toInt();
            if(record.satellite_number == 0)
            {

            }
            else if(record.satellite_number < 10)
            {
                for(int i = 0; i<record.satellite_number; i++)
                {
                    record.satellite_list.push_back(readString.mid(18+i*4,3));
                }
            }
            else
            {
                int i = ceil(record.satellite_number / 10.0); //有几行同样卫星的数据
                for(int j = 0; j<i; j++)
                {
                    for(int k = 0; k<13; k++)
                    {
                        record.satellite_list.push_back(readString.mid(18+k*4,3));
                    }
                    if(j < i-1)
                    {
                        readString = read.readLine();
                    }
                }
                int j = record.satellite_list.size()-record.satellite_number;
                for(int k = 0; k < j; k++)
                {
                    record.satellite_list.pop_back();
                }
            }
            heard_date.phase_shift_record.push_back(record);
        }
        else if(readString.indexOf("SIGNAL STRENGTH UNIT") >= 0)
        {
            heard_date.dbhz = readString.mid(0,20).simplified();
        }
    }while(readString.indexOf("END OF HEADER")<=0);

    QVector<int> sate_index={-1,-1,-1,-1,-1,-1};  /*C++11 new standard*/
    QVector<int> sate_num_index={0,0,0,0,0,0};
    /* GPS GLONASS SBAS Galileo BDS QZSS*/

    for(int i = 0; i<heard_date.system_record.size(); i++)
    {
        if(heard_date.system_record[i].system_type == "G")
        {
            sate_index[0] = i;
            sate_num_index[0] = heard_date.system_record[i].observation_number;
        }
        else if(heard_date.system_record[i].system_type == "R")
        {
            sate_index[1] = i;
            sate_num_index[1] = heard_date.system_record[i].observation_number;
        }
        else if(heard_date.system_record[i].system_type == "S")
        {
            sate_index[2] = i;
            sate_num_index[2] = heard_date.system_record[i].observation_number;
        }
        else if(heard_date.system_record[i].system_type == "E")
        {
            sate_index[3] = i;
            sate_num_index[3] = heard_date.system_record[i].observation_number;
        }
        else if(heard_date.system_record[i].system_type == "C")
        {
            sate_index[4] = i;
            sate_num_index[4] = heard_date.system_record[i].observation_number;
        }
        else if(heard_date.system_record[i].system_type == "J")
        {
            sate_index[5] = i;
            sate_num_index[5] = heard_date.system_record[i].observation_number;
        }
    }

    /*phase matching*/
    system_signal sys_signal;
    phase_matching(heard_date.system_record,sys_signal);

    /*date*/

    while(!read.atEnd())
    {
        o_epoch_date epoch;
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

        for(int i = 0; i<epoch.number_of_satellite; i++)
        {
            readString = read.readLine();
            o_sate_date satellite_date;
            satellite_date.satellite_infomation = readString.mid(0,3);
            /*use for QString to char* to ascii*/
            char *sate_info = NULL;
            QByteArray ba = satellite_date.satellite_infomation.mid(0,1).toLatin1();
            sate_info=ba.data();
            switch (toascii(sate_info[0]))
            {
            case 'G':
            {
                satellite_date.satellite_LLI.resize(6);
                satellite_date.satellite_observation_value.resize(6);
                satellite_date.satellite_signal_strength.resize(6);
                int G_num = sate_index[0];
                for(int j = 0;j<sate_num_index[G_num] ;j++)
                {
                    QString value = readString.mid(3+j*16,16);
                    if(sys_signal.GPS_P1 == j+1)
                    {
                        satellite_date.satellite_observation_value[0] = value.mid(0,14).toDouble();
                        satellite_date.satellite_LLI[0] = value.mid(14,1).toInt();
                        satellite_date.satellite_signal_strength[0] = value.mid(15,1).toInt();
                    }
                    else if(sys_signal.GPS_P2 == j+1)
                    {
                        satellite_date.satellite_observation_value[1] = value.mid(0,14).toDouble();
                        satellite_date.satellite_LLI[1] = value.mid(14,1).toInt();
                        satellite_date.satellite_signal_strength[1] = value.mid(15,1).toInt();
                    }
                    else if(sys_signal.GPS_P5 == j+1)
                    {
                        satellite_date.satellite_observation_value[2] = value.mid(0,14).toDouble();
                        satellite_date.satellite_LLI[2] = value.mid(14,1).toInt();
                        satellite_date.satellite_signal_strength[2] = value.mid(15,1).toInt();
                    }
                    else if(sys_signal.GPS_L1 == j+1)
                    {
                        satellite_date.satellite_observation_value[3] = value.mid(0,14).toDouble();
                        satellite_date.satellite_LLI[3] = value.mid(14,1).toInt();
                        satellite_date.satellite_signal_strength[3] = value.mid(15,1).toInt();
                    }
                    else if(sys_signal.GPS_L2 == j+1)
                    {
                        satellite_date.satellite_observation_value[4] = value.mid(0,14).toDouble();
                        satellite_date.satellite_LLI[4] = value.mid(14,1).toInt();
                        satellite_date.satellite_signal_strength[4] = value.mid(15,1).toInt();
                    }
                    else if(sys_signal.GPS_L5 == j+1)
                    {
                        satellite_date.satellite_observation_value[5] = value.mid(0,14).toDouble();
                        satellite_date.satellite_LLI[5] = value.mid(14,1).toInt();
                        satellite_date.satellite_signal_strength[5] = value.mid(15,1).toInt();
                    }

                }
                epoch.satellite_epoch.push_back(satellite_date);
                break;//GPS
            }
            case 'R':
            {
                int R_num = sate_index[1];
                break;//GLONASS
            }
            case 'S':
            {
                int S_num = sate_index[2];
                break;//SBAS
            }
            case 'E':
            {
                int E_num = sate_index[3];
                break;//Galileo
            }
            case 'C':
            {
                int C_num = sate_index[4];
                break;//BDS
            }
            case 'J':
            {
                int J_num = sate_index[5];
                break;//QZSS
            }
            default:
                break;
            }

        }
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
    sp3_heard_date heard_date;
    /*#C*/
    readString = read.readLine();
    heard_date.mode_flag  = readString.mid(2,1);
    heard_date.year       = readString.mid(3,4).toInt();
    heard_date.month      = readString.mid(7,3).toInt();
    heard_date.day        = readString.mid(10,3).toInt();
    heard_date.hour       = readString.mid(13,3).toInt();
    heard_date.minute     = readString.mid(16,3).toInt();
    heard_date.second     = readString.mid(19,12).toDouble();
    heard_date.number_of_epoch = readString.mid(31,8).toInt();
    heard_date.data_use   = readString.mid(39,6).simplified();
    heard_date.coordinate_system = readString.mid(45,6).simplified();
    heard_date.orbit_type = readString.mid(51,4).simplified();
    heard_date.agency     = readString.mid(55,5).simplified();

    /*##*/
    readString = read.readLine();
    heard_date.GPS_week_number  = readString.mid(3,4).toInt();
    heard_date.second_of_week   = readString.mid(7,16).toDouble();
    heard_date.interval         = readString.mid(23,15).toDouble();
    heard_date.modified_JD      = readString.mid(38,6).toInt();
    heard_date.fractional_day   = readString.mid(44,16).toDouble();

    /*+*/
    heard_date.satellites.reserve(85);
    for(int i = 0; i<5; i++)
    {
        readString = read.readLine();
        heard_date.number_of_satellites += readString.mid(4,2).toInt();
        for(int j = 0; j<17; j++)
        {
            heard_date.satellites.push_back(readString.mid(9+j*3,3));
        }
    }
    for(int k = 85; k>heard_date.number_of_satellites; k--)
    {
        heard_date.satellites.pop_back();
    }

    /*++*/
    heard_date.accuracy.reserve(85);
    for(int i=0; i<5; i++)
    {
        readString = read.readLine();
        for(int j = 0; j<17; j++)
        {
            heard_date.accuracy.push_back(readString.mid(9+j*3,3).toInt());
        }
    }
    for(int k = 85; k>heard_date.number_of_satellites; k--)
    {
        heard_date.accuracy.pop_back();
    }

    /*%c*/
    readString = read.readLine();
    heard_date.file_type = readString.mid(3,1);
    heard_date.time_type = readString.mid(9,3);
    readString = read.readLine();

    /*%f*/
    readString = read.readLine();
    heard_date.p_v_base = readString.mid(3 ,10).toDouble();
    heard_date.c_c_base = readString.mid(14,12).toDouble();
    for(int i = 0; i<7; i++)
    {
        read.readLine();
    }
    /*date*/

    for(int i = 0; i<heard_date.number_of_epoch; i++)
    {
        sp3_epoch_date epoch;
        readString = read.readLine();
        epoch.year   = readString.mid(3,4).toInt();
        epoch.month  = readString.mid(7,3).toInt();
        epoch.day    = readString.mid(10,3).toInt();
        epoch.hour   = readString.mid(13,3).toInt();
        epoch.minute = readString.mid(16,3).toInt();
        epoch.second = readString.mid(19,12).toDouble();
        for(int j = 0; j<heard_date.number_of_satellites; j++)
        {
            sp3_sate_date sate_date;
            readString = read.readLine();
            sate_date.flag  = readString.mid(0,1);
            sate_date.sate_info = readString.mid(1,3);
            sate_date.x     = readString.mid(4,14).toDouble();
            sate_date.y     = readString.mid(18,14).toDouble();
            sate_date.z     = readString.mid(32,14).toDouble();
            sate_date.clock = readString.mid(46,14).toDouble();
            sate_date.x_SD  = readString.mid(60,3).toInt();
            sate_date.y_SD  = readString.mid(63,3).toInt();
            sate_date.z_SD  = readString.mid(66,3).toInt();
            sate_date.clock_SD = readString.mid(69,4).toInt();
            epoch.epoch.push_back(sate_date);
        }
        sp3.file.push_back(epoch);
    }


}

void file_read::ppp_clock_read(const QString &file_path, clock_date &clock)
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

    clock_date_heard clock_heard;
    do
    {
        readString = read.readLine();
        if(readString.indexOf("RINEX VERSION / TYPE")>=0)
        {
            clock_heard.rinex_format_version = readString.mid(0,9);
            clock_heard.rinex_file_type      = readString.mid(20,1);
            clock_heard.satellite_system     = readString.mid(40,20).simplified();
        }
        else if(readString.indexOf("PGM / RUN BY / DATE")>=0)
        {
            clock_heard.creating_program_name = readString.mid(0,20).simplified();
            clock_heard.creating_agency_name  = readString.mid(20,20).simplified();
            clock_heard.creation_time         = readString.mid(40,20).simplified();
        }
        else if(readString.indexOf("# / TYPES OF DATA")>=0)
        {
            clock_heard.clock_type_num = readString.mid(0,6).toInt();
            for(int i = 0; i<clock_heard.clock_type_num; i++)
            {
                clock_heard.clock_type.push_back(readString.mid(6+i*6+6).simplified());
            }
        }
        else if(readString.indexOf("# OF SOLN STA / TRF")>=0)
        {
            clock_heard.reference_num   = readString.mid(0,6).toInt();
            clock_heard.satation_info.reserve(clock_heard.reference_num);
            clock_heard.reference_frame = readString.mid(10,50).simplified();
        }
        else if(readString.indexOf("SOLN STA NAME / NUM")>=0)
        {
            station_infomation info;
            info.station_name = readString.mid(0,5);
            info.station_id   = readString.mid(5,20);
            info.station_x    = readString.mid(25,12).toDouble();
            info.station_y    = readString.mid(37,12).toDouble();
            info.station_z    = readString.mid(49,11).toDouble();
            clock_heard.satation_info.push_back(info);
        }
        else if(readString.indexOf("# OF SOLN SATS")>=0)
        {
            clock_heard.satellite_num = readString.mid(0,6).toInt();
        }
        else if(readString.indexOf("PRN LIST")>=0)
        {
            int num = ceil(clock_heard.satellite_num / 15.0);
            for(int i = 0; i<num; i++)
            {
                for(int j = 0; j<15; j++)
                {
                    clock_heard.satellite_list.push_back(readString.mid(0+j*4,3));
                }
                readString = read.readLine();
            }
        }
    }while(readString.indexOf("END OF HEADER")<=0);


    /*read data*/
    do
    {
        readString = read.readLine();
        clock_info info;
        info.clock_type = readString.mid(0,2);
        if(info.clock_type != "AS")
        {
            continue;
        }
        info.R_S_name   = readString.mid(3,4);
        info.year       = readString.mid(8,4).toInt();
        info.month      = readString.mid(12,3).toInt();
        info.day        = readString.mid(15,3).toInt();
        info.hour       = readString.mid(18,3).toInt();
        info.minute     = readString.mid(21,3).toInt();
        info.second     = readString.mid(24,10).toDouble();
        info.number_of_data = readString.mid(34,6).toInt();
        if(info.number_of_data<=2)
        {
            for(int i = 0; i<info.number_of_data; i++)
            {
                info.record.push_back(readString.mid(40+20*i,20).toDouble());
            }
        }
        else
        {
            for(int i = 0; i<2; i++)
            {
                info.record.push_back(readString.mid(40+20*i,20).toDouble());
            }
            for(int i = 2; i<info.number_of_data; i++)
            {
                info.record.push_back(readString.mid(20*i,20).toDouble());
            }
        }
        clock.file.push_back(info);
    }while(!read.atEnd());
}

void file_read::phase_matching(const QVector<sys_record> &match_list, system_signal &sys_list)
{
    for(int i = 0; i<match_list.size(); i++)
    {
        if(match_list[i].system_type == "G")
        {
            sys_list.GPS_P1 = each_phase_matching(sys_list.GPS_P1_list,match_list[i].observation_descriptor);
            sys_list.GPS_P2 = each_phase_matching(sys_list.GPS_P2_list,match_list[i].observation_descriptor);
            sys_list.GPS_P5 = each_phase_matching(sys_list.GPS_P5_list,match_list[i].observation_descriptor);
            sys_list.GPS_L1 = each_phase_matching(sys_list.GPS_L1_list,match_list[i].observation_descriptor);
            sys_list.GPS_L2 = each_phase_matching(sys_list.GPS_L2_list,match_list[i].observation_descriptor);
            sys_list.GPS_L5 = each_phase_matching(sys_list.GPS_L5_list,match_list[i].observation_descriptor);
        }
        else if(match_list[i].system_type == "R")
        {
            sys_list.GLONASS_G1 = each_phase_matching(sys_list.GLONASS_G1_list,match_list[i].observation_descriptor);
            sys_list.GLONASS_G2 = each_phase_matching(sys_list.GLONASS_G2_list,match_list[i].observation_descriptor);
            sys_list.GLONASS_G3 = each_phase_matching(sys_list.GLONASS_G3_list,match_list[i].observation_descriptor);
            sys_list.GLONASS_L1 = each_phase_matching(sys_list.GLONASS_L1_list,match_list[i].observation_descriptor);
            sys_list.GLONASS_L2 = each_phase_matching(sys_list.GLONASS_L2_list,match_list[i].observation_descriptor);
            sys_list.GLONASS_L3 = each_phase_matching(sys_list.GLONASS_L3_list,match_list[i].observation_descriptor);
        }
        else if(match_list[i].system_type == "E")
        {
            sys_list.Galileo_E1  = each_phase_matching(sys_list.Galileo_E1_list, match_list[i].observation_descriptor);
            sys_list.Galileo_E5a = each_phase_matching(sys_list.Galileo_E5a_list,match_list[i].observation_descriptor);
            sys_list.Galileo_E5b = each_phase_matching(sys_list.Galileo_E5b_list,match_list[i].observation_descriptor);
            sys_list.Galileo_E5  = each_phase_matching(sys_list.Galileo_E5_list, match_list[i].observation_descriptor);
            sys_list.Galileo_E6  = each_phase_matching(sys_list.Galileo_E6_list, match_list[i].observation_descriptor);
            sys_list.Galileo_L1  = each_phase_matching(sys_list.Galileo_L1_list, match_list[i].observation_descriptor);
            sys_list.Galileo_L5a = each_phase_matching(sys_list.Galileo_L5a_list,match_list[i].observation_descriptor);
            sys_list.Galileo_L5b = each_phase_matching(sys_list.Galileo_L5b_list,match_list[i].observation_descriptor);
            sys_list.Galileo_L5  = each_phase_matching(sys_list.Galileo_L5_list, match_list[i].observation_descriptor);
            sys_list.Galileo_L6  = each_phase_matching(sys_list.Galileo_L6_list, match_list[i].observation_descriptor);
        }
        else if(match_list[i].system_type == "S")
        {
            sys_list.SBAS_P1 = each_phase_matching(sys_list.SBAS_P1_list,match_list[i].observation_descriptor);
            sys_list.SBAS_P5 = each_phase_matching(sys_list.SBAS_P5_list,match_list[i].observation_descriptor);
            sys_list.SBAS_L1 = each_phase_matching(sys_list.SBAS_L1_list,match_list[i].observation_descriptor);
            sys_list.SBAS_L5 = each_phase_matching(sys_list.SBAS_L5_list,match_list[i].observation_descriptor);
        }
        else if(match_list[i].system_type == "C")
        {
            sys_list.BDS_B1 = each_phase_matching(sys_list.BDS_B1_list,match_list[i].observation_descriptor);
            sys_list.BDS_B2 = each_phase_matching(sys_list.BDS_B2_list,match_list[i].observation_descriptor);
            sys_list.BDS_B3 = each_phase_matching(sys_list.BDS_B3_list,match_list[i].observation_descriptor);
            sys_list.BDS_L1 = each_phase_matching(sys_list.BDS_L1_list,match_list[i].observation_descriptor);
            sys_list.BDS_L2 = each_phase_matching(sys_list.BDS_L2_list,match_list[i].observation_descriptor);
            sys_list.BDS_L3 = each_phase_matching(sys_list.BDS_L3_list,match_list[i].observation_descriptor);
        }
        else if(match_list[i].system_type == "J")
        {
            sys_list.QZSS_P1 = each_phase_matching(sys_list.QZSS_P1_list,match_list[i].observation_descriptor);
            sys_list.QZSS_P2 = each_phase_matching(sys_list.QZSS_P2_list,match_list[i].observation_descriptor);
            sys_list.QZSS_P5 = each_phase_matching(sys_list.QZSS_P5_list,match_list[i].observation_descriptor);
            sys_list.QZSS_P6 = each_phase_matching(sys_list.QZSS_P6_list,match_list[i].observation_descriptor);
            sys_list.QZSS_L1 = each_phase_matching(sys_list.QZSS_L1_list,match_list[i].observation_descriptor);
            sys_list.QZSS_L2 = each_phase_matching(sys_list.QZSS_L2_list,match_list[i].observation_descriptor);
            sys_list.QZSS_L5 = each_phase_matching(sys_list.QZSS_L5_list,match_list[i].observation_descriptor);
            sys_list.QZSS_L6 = each_phase_matching(sys_list.QZSS_L6_list,match_list[i].observation_descriptor);
        }
    }
}

int file_read::each_phase_matching(const QStringList &reference_list, const QStringList &obs_descriptor)
{
    int i = reference_list.size();
    for(;i>0;i--)//zuijinque
    {
        int j = obs_descriptor.indexOf(reference_list[i-1]);
        if( j != -1 )
        {
            return j+1;
        }
    }
    return 0;
}

