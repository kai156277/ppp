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
    p_v_base(0),c_c_base(0)
{

}


sp3_sate_date::sp3_sate_date()
    :flag(""),sate_info(""),x(0),y(0),z(0),clock(0)
{

}


sp3_epoch_date::sp3_epoch_date()
    :year(0),month(0),day(0),
    hour(0),minute(0),second(0)
{

}


sp3_file::sp3_file()
{

}
