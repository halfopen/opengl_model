
// MFC_OpenGL.h : MFC_OpenGL Ӧ�ó������ͷ�ļ�
//
#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"       // ������


// CMFC_OpenGLApp:
// �йش����ʵ�֣������ MFC_OpenGL.cpp
//

class CMFC_OpenGLApp : public CWinApp
{
public:
	CMFC_OpenGLApp();


// ��д
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();
	

// ʵ��
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CMFC_OpenGLApp theApp;
