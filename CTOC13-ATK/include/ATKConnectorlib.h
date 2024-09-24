#pragma once

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

#if defined(_WIN64) || defined(WIN32)
#if defined(_WIN64) && defined(_DEBUG)
#pragma comment(lib, "./Lib/ATKConnectorlibD64.lib")
#elif defined(_WIN64)
#pragma comment(lib, "./Lib/ATKConnectorlib64.lib")
#elif defined(WIN32) && defined(_DEBUG)
#pragma comment(lib, "./Lib/ATKConnectorlibD.lib")
#elif defined(WIN32)
#pragma comment(lib, "./Lib/ATKConnectorlib.lib")
#endif

#else


#endif

typedef struct CmdResult
{
	std::vector<std::string> m_vectData;
	std::string Item(const int nIndex)
	{
		if (nIndex <= m_vectData.size())
		{
			return m_vectData[nIndex];
		}
	}
}CMDRESULT;

#define BUF_LENGTH 5000

void InitConnector();

int atkOpen(const char* szIP="127.0.0.1", const unsigned short& nPort=6655);

std::string atkConnect(const int conID, const std::string& strCommand, const std::string& strObjPath, const std::string& strCMDParam="");

void atkClose(const int conID);

std::string atkConnectEx(const int conID, const std::string& strCommand, const std::string& striInput);

int atkExecuteScript(std::string& striInput, std::string& strOutput);

CmdResult atkExecuteCommand(const int conID, const std::string& strCommand, const std::string& strObjPath, const std::string& strCMDParam = "");
