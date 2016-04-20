#include "clock_date.h"

clock_file::clock_file()
{

}



station_infomation::station_infomation()
    :station_name(""),station_id(""),station_x(0),station_y(0),station_z(0)
{

}


clock_date_heard::clock_date_heard()
    :rinex_format_version(""),rinex_file_type(""),satellite_system(""),
    creating_program_name(""),creating_agency_name(""),creation_time(""),
    clock_type_num(0),reference_num(0),reference_frame(""),satellite_num(0)

{

}



clock_info::clock_info()
    :sate_name("")

{

}


clock_epoch::clock_epoch()
    :year(0),month(0),day(0),
     hour(0),minute(0),second(0),GPSS(0),GPSW(0),number_of_data(0)
{
    GPS_epoch.reserve(30);
    BDS_epoch.reserve(20);
    GLONASS_epoch.reserve(30);
}
