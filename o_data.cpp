#include "o_data.h"
#include<QString>
#include<QVector>
#include<QList>
#include<QStringList>

o_heard::o_heard()
    :format_version(""),file_type(""),
    creating_program_name(""),creating_agency_name(""),creation_time(""),
    marker_name(""),marker_number(""),marke_type(""),
    observer_name(""),agency_name(""),
    receiver_number(""),receiver_type(""),receiver_version(""),
    antenna_number(""),antenna_type(""),
    position_X(0),position_Y(0),position_Z(0),
    antenna_H(0),antenna_E(0),antenna_N(0),
    interval(0),first_time(""),time_system(""),dbhz(""),GPS_number(0),observation_descriptor()
{

}

o_sate::o_sate()
    :satellite_infomation("")
{
    satellite_observation_value.resize(6);
    satellite_LLI.resize(6);
    satellite_signal_strength.resize(6);
}

o_epoch::o_epoch()
    :year(0),month(0),day(0),
    hour(0),minute(0),second(0),
    epoch_flag(0),number_of_satellite(0),clock_offset(0),
    GPSW(0),GPSS(0)
{

}

o_file::o_file()
{

}

phase_shift::phase_shift()
    :sate_type(""),cpoc(""),correction(0),satellite_number(0)
{

}


system_signal::system_signal()
    :GPS_P1(0),GPS_P2(0),GPS_P5(0),
     GPS_L1(0),GPS_L2(0),GPS_L5(0)


{
    GPS_P1_list << "C1C" << "C1S" << "C1L" << "C1X" << "C1P" << "C1W" << "C1Y" << "C1M";
    GPS_P2_list << "C2D" << "C2S" << "C2L" << "C2X" << "C2P" << "C2W" << "C2Y" << "C2M";
    GPS_P5_list << "C5I" << "C5Q" << "C5X";
    GPS_L1_list << "L1C" << "L1S" << "L1L" << "L1X" << "L1P" << "L1W" << "L1Y" << "L1M";
    GPS_L2_list << "L2D" << "L2S" << "L2L" << "L2X" << "L2P" << "L2W" << "L2Y" << "L2M";
    GPS_L5_list << "L5I" << "L5Q" << "L5X";
}
