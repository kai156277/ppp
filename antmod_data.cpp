#include "antmod_data.h"


satellite_antmod::satellite_antmod()
    :antenna_type(""),sate_info(""),
    DAZI(0),ZEN1(0),ZEN2(0),DZEN(0), num_of_frequencies(0),
    start_time(""),end_time(""),frequency_type(""),
    L1_frequency("G01"),L1_APC_x(0),L1_APC_y(0),L1_APC_z(0),
    L2_frequency("G02"),L2_APC_x(0),L2_APC_y(0),L2_APC_z(0)
{
    L1_NOAZI.reserve(14);
    L2_NOAZI.reserve(14);
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
