#include "antmod_data.h"


satellite_antmod::satellite_antmod()
    :antenna_type(""),sate_info(""),
    DAZI(0),ZEN1(0),ZEN2(0),DZEN(0), num_of_frequencies(0),
    start_time(""),end_time(""),frequency_type(""),
    start_year(0),start_month(0),start_day(0),start_hour(0),start_minute(0),start_second(0),
    end_year(0),end_month(0),end_day(0),end_hour(0),end_minute(0),end_second(0),
    frequency("G01"),APC_x(0),APC_y(0),APC_z(0)
{
    NOAZI.reserve(14);
}

bool satellite_antmod::operator ==(const satellite_antmod &right) const
{
    if(this->sate_info == right.sate_info)
    {
        return true;
    }
    else
    {
        return false;
    }
}


station_antmod::station_antmod()
    :antenna_type(""),receiver_info(""),
    DAZI(0),ZEN1(0),ZEN2(0),DZEN(0), num_of_frequencies(0),
    start_time(""),end_time(""),frequency_type(""),
    L1_frequency("G01"),L1_APC_N(0),L1_APC_E(0),L1_APC_U(0),
    L2_frequency("G02"),L2_APC_N(0),L2_APC_E(0),L2_APC_U(0)
{

}


antmod_file::antmod_file()
{

}
