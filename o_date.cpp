#include "o_date.h"
#include<QString>
#include<QVector>
#include<QList>
#include<QStringList>

sys_record::sys_record()
    :system_type(""),observation_number(0)
{

}

o_heard_date::o_heard_date()
    :format_version(""),file_type(""),
    creating_program_name(""),creating_agency_name(""),creation_time(""),
    marker_name(""),marker_number(""),marke_type(""),
    observer_name(""),agency_name(""),
    receiver_number(""),receiver_type(""),receiver_version(""),
    antenna_number(""),antenna_type(""),
    position_X(0),position_Y(0),position_Z(0),
    antenna_H(0),antenna_E(0),antenna_N(0),
    interval(0),first_time(""),time_system(""),dbhz("")
{

}

o_sate_date::o_sate_date()
    :satellite_infomation("")
{
    satellite_observation_value.reserve(6);
    satellite_LLI.reserve(6);
    satellite_signal_strength.reserve(6);
}

o_epoch_date::o_epoch_date()
    :year(0),month(0),day(0),
    hour(0),minute(0),second(0),
    epoch_flag(0),number_of_satellite(0),clock_offset(0),
    GPSW(0),GPSS(0)
{

}

o_file_date::o_file_date()
{

}

phase_shift::phase_shift()
    :sate_type(""),cpoc(""),correction(0),satellite_number(0)
{

}


system_signal::system_signal()
    :GPS_P1(0),GPS_P2(0),GPS_P5(0),
     GPS_L1(0),GPS_L2(0),GPS_L5(0),
     GLONASS_G1(0),GLONASS_G2(0),GLONASS_G3(0),
     GLONASS_L1(0),GLONASS_L2(0),GLONASS_L3(0),
     BDS_B1(0),BDS_B2(0),BDS_B3(0),
     BDS_L1(0),BDS_L2(0),BDS_L3(0)
//     Galileo_E1(0),Galileo_E5a(0),Galileo_E5b(0),Galileo_E5(0),Galileo_E6(0),
//     Galileo_L1(0),Galileo_L5a(0),Galileo_L5b(0),Galileo_L5(0),Galileo_L6(0),
//     SBAS_P1(0),SBAS_P5(0),
//     SBAS_L1(0),SBAS_L5(0),
//     QZSS_P1(0),QZSS_P2(0),QZSS_P5(0),QZSS_P6(0),
//     QZSS_L1(0),QZSS_L2(0),QZSS_L5(0),QZSS_L6(0),

{
    GPS_P1_list << "C1C" << "C1S" << "C1L" << "C1X" << "C1P" << "C1W" << "C1Y" << "C1M";
    GPS_P2_list << "C2D" << "C2S" << "C2L" << "C2X" << "C2P" << "C2W" << "C2Y" << "C2M";
    GPS_P5_list << "C5I" << "C5Q" << "C5X";
    GPS_L1_list << "L1C" << "L1S" << "L1L" << "L1X" << "L1P" << "L1W" << "L1Y" << "L1M";
    GPS_L2_list << "L2D" << "L2S" << "L2L" << "L2X" << "L2P" << "L2W" << "L2Y" << "L2M";
    GPS_L5_list << "L5I" << "L5Q" << "L5X";
    GLONASS_G1_list << "C1C" << "C1P";
    GLONASS_G2_list << "C2C" << "C2P";
    GLONASS_G3_list << "C3I" << "C3Q" << "C3X";
    GLONASS_L1_list << "L1C" << "L1P";
    GLONASS_L2_list << "L2C" << "L2P";
    GLONASS_L3_list << "L3I" << "L3Q" << "C3X";
    BDS_B1_list << "C1I" << "C1Q" << "C1X";
    BDS_B2_list << "C7I" << "C7Q" << "C7X";
    BDS_B3_list << "C6I" << "C6Q" << "C6X";
    BDS_L1_list << "L1I" << "L1Q" << "L1X";
    BDS_L2_list << "L7I" << "L7Q" << "L7X";
    BDS_L3_list << "L6I" << "L6Q" << "L6X";
//    Galileo_E1_list  << "C1A" << "C1B" << "C1C" << "C1X" << "C1Z";
//    Galileo_E5a_list << "C5I" << "C5Q" << "C5X";
//    Galileo_E5b_list << "C7I" << "C7Q" << "C7X";
//    Galileo_E5_list  << "C8I" << "C8Q" << "C8X";
//    Galileo_E6_list  << "C6A" << "C6B" << "C6C" << "C6X" << "C6Z";
//    Galileo_L1_list  << "L1A" << "L1B" << "L1C" << "L1X" << "L1Z";
//    Galileo_L5a_list << "L5I" << "L5Q" << "L5X";
//    Galileo_L5b_list << "L7I" << "L7Q" << "L7X";
//    Galileo_L5_list  << "L8I" << "L8Q" << "L8X";
//    Galileo_L6_list  << "L6A" << "L6B" << "L6C" << "L6X" << "L6Z";
//    SBAS_P1_list << "C1C";
//    SBAS_P5_list << "C5I" << "C5Q" << "C5X";
//    SBAS_L1_list << "L1C";
//    SBAS_L5_list << "L5I" << "L5Q" << "L5X";
//    QZSS_P1_list << "C1C" << "C1S" << "C1L" << "C1X" << "C1Z";
//    QZSS_P2_list << "C2S" << "C2L" << "C2X";
//    QZSS_P5_list << "C5I" << "C5Q" << "C5X";
//    QZSS_P6_list << "C6S" << "C6L" << "C6X";
//    QZSS_L1_list << "L1C" << "L1S" << "L1L" << "L1X" << "L1Z";
//    QZSS_L2_list << "L2S" << "L2L" << "L2X";
//    QZSS_L5_list << "L5I" << "L5Q" << "L5X";
//    QZSS_L6_list << "L6S" << "L6L" << "L6X";

}
