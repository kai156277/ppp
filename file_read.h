#ifndef FILE_READ_H
#define FILE_READ_H

#include<QString>

#include"o_data.h"
#include"sp3_data.h"
#include"clock_data.h"
#include"snx_data.h"
#include"antmod_data.h"

class file_read
{
public:
    file_read();
    void ppp_o_read(const QString &file_path,o_file &obs);
    void ppp_sp3_read(const QString &file_path,sp3_file &sp3);
    void ppp_clock_read(const QString &file_path,clock_file &clock);
    void ppp_snx_read(const QString &file_path,snx_data &snx,QString mark_name);
    void ppp_ant_read(const QString &file_path,antmod_file &ant);
private:
    void phase_matching(const QStringList &match_list, system_signal &sys_list);//选择观测信号
    int each_phase_matching(const QStringList & reference_list,const QStringList & obs_descriptor);//匹配最精确的观测信号
};

#endif // FILE_READ_H
