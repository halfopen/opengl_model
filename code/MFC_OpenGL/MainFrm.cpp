
// MainFrm.cpp : CMainFrame ���ʵ��
//

#include "stdafx.h"
#include "MFC_OpenGL.h"

#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMainFrame

IMPLEMENT_DYNAMIC(CMainFrame, CMDIFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CMDIFrameWnd)
	ON_WM_CREATE()
	ON_COMMAND(ID_WINDOW_MYTILE_VERT, &CMainFrame::OnWindowMytileVert)
	ON_COMMAND(ID_WINDOW_MYCASCADE, &CMainFrame::OnWindowMycascade)
	ON_COMMAND(ID_WINDOW_MYARRANGE, &CMainFrame::OnWindowMyarrange)
END_MESSAGE_MAP()

static UINT indicators[] =
{
	ID_SEPARATOR,           // ״̬��ָʾ��
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

// CMainFrame ����/����

CMainFrame::CMainFrame()
{
	// TODO:  �ڴ���ӳ�Ա��ʼ������
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CMDIFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;

	if (!m_wndToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP | CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC) ||
		!m_wndToolBar.LoadToolBar(IDR_MAINFRAME))
	{
		TRACE0("δ�ܴ���������\n");
		return -1;      // δ�ܴ���
	}

	if (!m_wndStatusBar.Create(this))
	{
		TRACE0("δ�ܴ���״̬��\n");
		return -1;      // δ�ܴ���
	}
	m_wndStatusBar.SetIndicators(indicators, sizeof(indicators)/sizeof(UINT));

	// TODO:  �������Ҫ��ͣ������������ɾ��������
	m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
	EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_wndToolBar);


	return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CMDIFrameWnd::PreCreateWindow(cs) )
		return FALSE;
	// TODO:  �ڴ˴�ͨ���޸�
	//  CREATESTRUCT cs ���޸Ĵ��������ʽ
	//cs.cx = 100;
	//cs.cy = 100;
	//cs.x = 400;
	//cs.y = 400;
	return TRUE;
}

// CMainFrame ���

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CMDIFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CMDIFrameWnd::Dump(dc);
}
#endif //_DEBUG


// CMainFrame ��Ϣ�������



//void CMainFrame::OnWindowCascade()
//{
//	// TODO:  �ڴ���������������
//	SkinH_Detach();
//	this->MDITile(MDITILE_SKIPDISABLED);
//	SkinH_Attach();
//}

//ƽ��
void CMainFrame::OnWindowMytileVert()
{
	// TODO:  �ڴ���������������
	SkinH_Detach();
	this->MDITile(MDITILE_VERTICAL);
	OutputDebugString(_T("ver"));
	SkinH_Attach();
}

//���
void CMainFrame::OnWindowMycascade()
{
	// TODO:  �ڴ���������������
	// TODO:  �ڴ���������������
	SkinH_Detach();
	this->MDITile(MDITILE_HORIZONTAL);
	SkinH_Attach();
}


void CMainFrame::OnWindowMyarrange()
{
	// TODO:  �ڴ���������������
	SkinH_Detach();
	//this->MDITile(MDITILE_ZORDER);
	SkinH_Attach();
	OutputDebugString(_T("arrange\n"));
}
