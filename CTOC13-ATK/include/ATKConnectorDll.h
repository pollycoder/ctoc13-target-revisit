#ifndef ATKConnectorDll_h
#define ATKConnectorDll_h
// 下列 ifdef 块是创建使从 DLL 导出更简单的
// 宏的标准方法。此 DLL 中的所有文件都是用命令行上定义的 ATKCONNECTORDLL_EXPORTS
// 符号编译的。在使用此 DLL 的
// 任何其他项目上不应定义此符号。这样，源文件中包含此文件的任何其他项目都会将
// PLUGINDLLTEST_API 函数视为是从 DLL 导入的，而此 DLL 则将用此宏定义的
// 符号视为是被导出的。

#if defined(_WIN64) || defined(WIN32)
#include <xstring>
#include <vector>
#include <process.h>
#include <WinSock2.h>

#else
#include <string>
#include <vector>
#include <semaphore.h>
//#include <QtCore/qglobal.h>
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <linux/rtc.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN64) || defined(WIN32)

#ifdef ATKCONNECTORDLL_EXPORTS
#define ATKCONNECTORDLL_API __declspec(dllexport)
#else
#define ATKCONNECTORDLL_API __declspec(dllimport)
#if defined(_WIN64) && defined(_DEBUG)
#pragma comment(lib, "ATKConnectorDllD64.lib")
#elif defined(_WIN64)
#pragma comment(lib, "ATKConnectorDll64.lib")
#elif defined(WIN32) && defined(_DEBUG)
#pragma comment(lib, "ATKConnectorDllD.lib")
#elif defined(WIN32)
#pragma comment(lib, "ATKConnectorDll.lib")
#endif
#endif

#else

#ifdef ATKCONNECTORDLL_LIBRARY
#define ATKCONNECTORDLL_API __attribute__((visibility("default")))
#else
#define ATKCONNECTORDLL_API __attribute__((visibility("default")))
#endif

#endif

	typedef struct CmdResult
	{
		std::vector<std::string> m_vectData;
		std::string Item(const int nIndex)
		{
			if (nIndex < m_vectData.size())
			{
				return m_vectData[nIndex];
			}
		}
	}CMDRESULT;

#define BUF_LENGTH 5000

	ATKCONNECTORDLL_API int atkOpen(const char* szIP = "127.0.0.1", const unsigned short& nPort = 6655);

	ATKCONNECTORDLL_API std::string atkConnect(const int conID, const std::string& strCommand, const std::string& strObjPath, const std::string& strCMDParam = "");

	ATKCONNECTORDLL_API void atkClose(const int conID);

	ATKCONNECTORDLL_API std::string atkConnectEx(const int conID, const std::string& strCommand, const std::string& striInput);

	ATKCONNECTORDLL_API int atkExecuteScript(std::string& striInput, std::string& strOutput);

	ATKCONNECTORDLL_API CmdResult atkExecuteCommand(const int conID, const std::string& strCommand, const std::string& strObjPath, const std::string& strCMDParam = "");

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

#endif /* ATKConnectorDll_h */
