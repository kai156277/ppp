#include "sp3_date.h"

sp3_heard_date::sp3_heard_date()
    :mode_flag(""),
    year(0),month(0),day(0),
    hour(0),minute(0),second(0),
    number_of_epoch(0),
    data_use(""),coordinate_system(""),orbit_type(""),agency(""),
    GPS_week_number(0),second_of_week(0),interval(0),
    modified_JD(0),fractional_day(0),number_of_satellites(0),
    file_type(""),time_type(""),
    p_v_base(0),c_c_base(0),
    BDS_satellites(0),Galileo_satellites(0),GPS_satellites(0),
    QZSS_satellites(0),GLONASS_satellites(0),SBAS_satellites(0)
{

}


sp3_sate_date::sp3_sate_date()
    :flag(""),sate_info(""),x(0),y(0),z(0),clock(0),x_SD(0),y_SD(0),z_SD(0),clock_SD(0)
{

}

bool sp3_sate_date::operator ==(const o_sate_date &left) const
{
    if(this->sate_info == left.satellite_infomation)
    {
        return true;
    }
    else
    {
        return false;
    }
}


sp3_epoch_date::sp3_epoch_date()
    :year(0),month(0),day(0),
    hour(0),minute(0),second(0),
    GPSW(0),GPSS(0)
{

}


sp3_file::sp3_file()
    :heard()
{

}
