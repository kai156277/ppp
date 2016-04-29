#include "ocean_date.h"

ocean_file::ocean_file()
    :col(0)
{

}



ocean::ocean()
    :stationName(""),oceanModel(""),
      year(0),month(0),day(0),
      station_B(0),station_L(0),station_H(0)
{
    data[6][11] = {0};
}
