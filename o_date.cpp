#include "o_date.h"
#include<QString>
#include<QVector>
#include<QList>
#include<QStringList>

system_record::system_record()
    :system_type("")
    ,observation_descriptor(0)
{

}

o_heard_date::o_heard_date()
    :format_version(""),file_type(""),satellite_system(""),
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

sate_date::sate_date()
    :satellite_infomation("")
{
    satellite_observation_value.reserve(6);
}

epoch_date::epoch_date()
    :year(0),month(0),day(0),
    hour(0),minute(0),second(0),
    epoch_flag(0),number_of_satellite(0),clock_offset(0)
{

}

file_date::file_date()
{

}

const char* obs_signal::P1_list = "C1C C1S C1L C1X C1P C1W C1Y C1M";
const char* obs_signal::P2_list = "C2D C2S C2L C2X C2P C2W C2Y C2M";
const char* obs_signal::P5_list = "C5I C5Q C5X";
const char* obs_signal::L1_list = "L1C L1S L1L L1X L1P L1W L1Y L1M L1N";
const char* obs_signal::L2_list = "L2D L2S L2L L2X L2P L2W L2Y L2M L2N";
const char* obs_signal::L5_list = "L5I L5Q L5X";

obs_signal::obs_signal()
{

}
