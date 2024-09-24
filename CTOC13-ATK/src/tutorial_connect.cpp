
#include "ATKConnectorDll.h" //包含头文件
#include <string>
#include <iostream>
#include <sstream>
using namespace std;

int tutorial_connect() {
    // 1.连接ATK软件
    int conID = atkOpen();

    // 2.添加卫星与敏感器
    if (atkConnect(conID, "DoesObjExist", "/ */Satellite/Satellite1") == "0\n")
    {
        atkConnect(conID, "New", "/ Satellite Satellite1"); // 新建卫星1
        atkConnect(conID, "New", "/ */Satellite/Satellite1/Sensor Sensor1"); // 新建敏感器1
        atkConnect(conID, "Define", "*/Satellite/Satellite1/Sensor/Sensor1 Rectangular 20 20"); // 定义敏感器1参数
    }
    if (atkConnect(conID, "DoesObjExist", "/ */Satellite/Satellite2") == "0\n")
    {
        atkConnect(conID, "New", "/ Satellite Satellite2"); // 新建卫星2
        atkConnect(conID, "New", "/ */Satellite/Satellite2/Sensor Sensor2"); // 新建敏感器2
        atkConnect(conID, "Define", "*/Satellite/Satellite2/Sensor/Sensor2 Rectangular 20 20"); // 定义敏感器2参数
    }
    if (atkConnect(conID, "DoesObjExist", "/ */Satellite/Satellite3") == "0\n")
    {
        atkConnect(conID, "New", "/ Satellite Satellite3");  // 新建卫星3
        atkConnect(conID, "New", "/ */Satellite/Satellite3/Sensor Sensor3"); // 新建敏感器3
        atkConnect(conID, "Define", "*/Satellite/Satellite3/Sensor/Sensor3 Rectangular 20 20"); // 定义敏感器3参数
    }

    // 3.修改卫星参数
    /// 给定卫星动力学模型、坐标系及历元
    string t1 = "*/Satellite/Satellite1 Classical J2Numerical \"26 Sep 2035 12:00:00.00\" \"28 Sep 2035 12:00:00.00\" 60 j2000 \"26 Sep 2035 12:00:00.00\" ";
    string t2 = "*/Satellite/Satellite2 Classical J2Numerical \"26 Sep 2035 12:00:00.00\" \"28 Sep 2035 12:00:00.00\" 60 j2000 \"26 Sep 2035 12:00:00.00\" ";
    string t3 = "*/Satellite/Satellite3 Classical J2Numerical \"26 Sep 2035 12:00:00.00\" \"28 Sep 2035 12:00:00.00\" 60 j2000 \"26 Sep 2035 12:00:00.00\" ";
    /// 给定卫星轨道根数
    string space, a, e,
        i1, Raan1, Perigee1, True1,
        i2, Raan2, Perigee2, True2,
        i3, Raan3, Perigee3, True3;
    space = " ";
    a = "7348136.00000014"; e = "0";
    i1 = "10"; Raan1 = "10"; Perigee1 = "10"; True1 = "10";
    i2 = "20"; Raan2 = "20"; Perigee2 = "20"; True2 = "20";
    i3 = "30"; Raan3 = "30"; Perigee3 = "30"; True3 = "30";
    string sat1 = t1 + a + space + e + space + i1 + space + Raan1 + space + Perigee1 + space + True1;
    string sat2 = t2 + a + space + e + space + i2 + space + Raan2 + space + Perigee2 + space + True2;
    string sat3 = t3 + a + space + e + space + i3 + space + Raan3 + space + Perigee3 + space + True3;
    /// 根据给定参数设置各卫星状态
    atkConnect(conID, "SetState", sat1);
    atkConnect(conID, "SetState", sat2);
    atkConnect(conID, "SetState", sat3);

    // 4.计算最大重访时间
    /// 多对多覆盖性分析设置
    atkConnect(conID, "CovMulti", "/ Assets */Satellite/Satellite1/Sensor/Sensor1 */Satellite/Satellite2/Sensor/Sensor2 */Satellite/Satellite3/Sensor/Sensor3");// 选择卫星敏感器
    atkConnect(conID, "CovMulti", "/ objects */Facility/Target1 */Facility/Target2 */Facility/Target3"); // 选择地面目标
    atkConnect(conID, "CovMulti", "/ Access Compute \"26 Sep 2035 12:00:00.00\" \"28 Sep 2035 12:00:00.00\""); // 计算指定时间段的覆盖性
    string result = atkConnect(conID, "CovMulti_RM", "/ multifomdefine Definition RevisitTime Compute Maximum"); // 计算最大重访时间品质值
    /// 从返回的result字符串中获取每个目标的最大重访时间
    bool Judge = false;
    vector<string> OutNumStrVec(3, "");
    int i = 0;
    int temp = 0;
    for (char c : result)
    {
        if (c == ':')
            Judge = true;
        else if (Judge)
        {
            OutNumStrVec.at(i) += c;
            temp += 1;
            if (temp == 12)
            {
                Judge = false;
                temp = 0;
                i++;
            }
        }
    }
    vector<double> maxRevisitTime;
    maxRevisitTime.resize(3);
    for (int i = 0; i < 3; i++)
    {
        int Num = 0;
        stringstream ss(OutNumStrVec.at(i));
        ss >> Num;
        maxRevisitTime[i] = Num;
    }
    // 输出每个目标的最大重访时间（-1代表该目标从未被重访）
    cout << "Target1MaxRevisitTime = " << maxRevisitTime[0] << endl;
    cout << "Target2MaxRevisitTime = " << maxRevisitTime[1] << endl;
    cout << "Target3MaxRevisitTime = " << maxRevisitTime[2] << endl;

    // 5.断开ATK软件
    atkClose(conID);

	return 0;
}