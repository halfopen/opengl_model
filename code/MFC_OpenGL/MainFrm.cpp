
// MainFrm.cpp : CMainFrame 类的实现
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
	ID_SEPARATOR,           // 状态行指示器
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

// CMainFrame 构造/析构

CMainFrame::CMainFrame()
{
	// TODO:  在此添加成员初始化代码
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
		TRACE0("未能创建工具栏\n");
		return -1;      // 未能创建
	}

	if (!m_wndStatusBar.Create(this))
	{
		TRACE0("未能创建状态栏\n");
		return -1;      // 未能创建
	}
	m_wndStatusBar.SetIndicators(indicators, sizeof(indicators)/sizeof(UINT));

	// TODO:  如果不需要可停靠工具栏，则删除这三行
	m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
	EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_wndToolBar);


	return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CMDIFrameWnd::PreCreateWindow(cs) )
		return FALSE;
	// TODO:  在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式
	//cs.cx = 100;
	//cs.cy = 100;
	//cs.x = 400;
	//cs.y = 400;
	return TRUE;
}

// CMainFrame 诊断

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


// CMainFrame 消息处理程序



//void CMainFrame::OnWindowCascade()
//{
//	// TODO:  在此添加命令处理程序代码
//	SkinH_Detach();
//	this->MDITile(MDITILE_SKIPDISABLED);
//	SkinH_Attach();
//}

//平铺
void CMainFrame::OnWindowMytileVert()
{
	// TODO:  在此添加命令处理程序代码
	SkinH_Detach();
	this->MDITile(MDITILE_VERTICAL);
	OutputDebugString(_T("ver"));
	SkinH_Attach();
}

//层叠
void CMainFrame::OnWindowMycascade()
{
	// TODO:  在此添加命令处理程序代码
	// TODO:  在此添加命令处理程序代码
	SkinH_Detach();
	this->MDITile(MDITILE_HORIZONTAL);
	SkinH_Attach();
}


void CMainFrame::OnWindowMyarrange()
{
	// TODO:  在此添加命令处理程序代码
	SkinH_Detach();
	//this->MDITile(MDITILE_ZORDER);
	SkinH_Attach();
	OutputDebugString(_T("arrange\n"));
}
