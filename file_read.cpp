#include<QString>
#include<QVector>
#include<QStringList>
#include<QFile>
#include<QTextStream>
#include<qdebug.h>
#include<cmath>

#include "file_read.h"
#include "o_date.h"

using namespace std;

file_read::file_read()
{

}

void file_read::ppp_o_read(QString file_path, o_file_date &obs)
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

    /*phase matching*/
    system_signal sys_signal;
    phase_matching(heard_date.system_record,sys_signal);

    /*date*/
    o_epoch_date epoch;
    while(!read.atEnd())
    {

    }
    ppp_o_file.close();
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

