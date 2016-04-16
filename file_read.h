#ifndef FILE_READ_H
#define FILE_READ_H

#include<QString>

#include"o_date.h"
#include"sp3_date.h"
#include"clock_date.h"

class file_read
{
public:
    file_read();
    void ppp_o_read(const QString &file_path,o_file_date &obs);
    void ppp_sp3_read(const QString &file_path,sp3_file &sp3);
    void ppp_clock_read(const QString &file_path,clock_date &clock);
private:
    void phase_matching(const QVector<sys_record> &match_list, system_signal &sys_list);//选择观测信号
    int each_phase_matching(const QStringList & reference_list,const QStringList & obs_descriptor);//匹配最精确的观测信号
};

#endif // FILE_READ_H
