#include"max_reseetime.h"

//t_start t_end都是s
std::vector<double> max_reseetime(std::vector<std::vector<double>> reseetime_all,double t_start,double t_end, std::vector<double>& reseetime_max)
{
   
    std::vector<std::vector<double>> c;
   
    //std::vector<std::vector<double>> c;

    for (const auto& row : reseetime_all) {
        if (row.empty()) {
            c.push_back({(t_end - t_start)/3600.});
            continue;
        }

        std::vector<double> diffRow;
        diffRow.push_back((row[0] - t_start) / 3600.); // 第一个值减0
        //std::cout << diffRow[0] << " ";
        for (size_t j = 1; j <= row.size(); ++j) {
            if (j < row.size())
                diffRow.push_back((row[j] - row[j - 1]) / 3600.); // 相邻值的差
            else
                diffRow.push_back((t_end - row.back())/3600.);
             
           // std::cout << diffRow[j] << " ";
        }
      
      //  std::cout << "\n ";
        c.push_back(diffRow);
    }
    for (const auto& row : c) {
        if (!row.empty()) {
            double maxVal = *std::max_element(row.begin(), row.end());
            reseetime_max.push_back(maxVal);
        }
    }
    std::vector<double> result(reseetime_all.size(), 0.0);
    result[0] = 1;
    return result;
}

std::vector<double> max_reseetime(std::vector<std::vector<double>> reseetime_all, std::vector<double>& reseetime_max)
{
    double T_2day = 2. * 86400.;
    std::vector<std::vector<double>> c;

    for (const auto& row : reseetime_all) {
        if (row.empty()) {
            c.push_back({ T_2day/3600. });
            continue;
        }
        std::vector<double> diffRow;
        diffRow.push_back((row[0] - 0.0) / 3600.); // 第一个值减0
      
        for (size_t j = 1; j <= row.size(); ++j) {
            if (j < row.size())
                diffRow.push_back((row[j] - row[j - 1]) / 3600.); // 相邻值的差
            else
                diffRow.push_back((T_2day - row.back()) / 3600.);
        }
        c.push_back(diffRow);
    }
    for (const auto& row : c)
    {
        if (!row.empty())
        {
            double maxVal = *std::max_element(row.begin(), row.end());
            reseetime_max.push_back(maxVal);
        }
    }
    std::vector<double> result(reseetime_all.size(), 0.0);
    result[0] = 1;
    return result;
}


std::vector<double> max_reseetime(std::vector<std::vector<uint16_t>> reseetime_all, std::vector<double>& reseetime_max)
{
    double T_2day = 2. * 86400./60.0;
    std::vector<std::vector<double>> c;

    for (const auto& row : reseetime_all) {
        if (row.empty()) {
            c.push_back({ T_2day / 60. });
            continue;
        }
        std::vector<double> diffRow;
        diffRow.push_back((row[0] - 0.0) / 60.); // 第一个值减0

        for (size_t j = 1; j <= row.size(); ++j) {
            if (j < row.size())
                diffRow.push_back((row[j] - row[j - 1]) / 60.); // 相邻值的差
            else
                diffRow.push_back((T_2day - row.back()) / 60.);
        }
        c.push_back(diffRow);
    }
    for (const auto& row : c)
    {
        if (!row.empty())
        {
            double maxVal = *std::max_element(row.begin(), row.end());
            reseetime_max.push_back(maxVal);
        }
    }
    std::vector<double> result(reseetime_all.size(), 0.0);
    result[0] = 1;
    return result;
}

// 时刻表的单位：秒
void max_revisit_interval_beforeEnd(std::vector<double>& reseetime_max, const std::vector<std::vector<double>>& reseetime_all) {
    reseetime_max.clear();
    reseetime_max.resize(TargetNum);

    // 计算每个目标的重访时间
    for (int i = 0; i < TargetNum; i++) {
        double max_gap = 0.0;
        if (reseetime_all[i].empty()) continue;
        for (int j = 0; j < reseetime_all[i].size() - 1; j++) {
            double gap_temp = reseetime_all[i][j + 1] - reseetime_all[i][j];
            if (gap_temp > max_gap) max_gap = gap_temp;
        }

        reseetime_max[i] = max_gap;
    }
}