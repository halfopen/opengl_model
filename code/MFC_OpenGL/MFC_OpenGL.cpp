
// MFC_OpenGL.cpp : 定义应用程序的类行为。
//

#include "stdafx.h"
#include "afxwinappex.h"
#include "afxdialogex.h"
#include "MFC_OpenGL.h"
#include "MainFrm.h"

#include "ChildFrm.h"
#include "MFC_OpenGLDoc.h"
#include "MFC_OpenGLView.h"
#include "afxwin.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMFC_OpenGLApp

BEGIN_MESSAGE_MAP(CMFC_OpenGLApp, CWinApp)
	ON_COMMAND(ID_APP_ABOUT, &CMFC_OpenGLApp::OnAppAbout)
	// 基于文件的标准文档命令
	ON_COMMAND(ID_FILE_NEW, &CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, &CWinApp::OnFileOpen)
	// 标准打印设置命令
	ON_COMMAND(ID_FILE_PRINT_SETUP, &CWinApp::OnFilePrintSetup)
END_MESSAGE_MAP()


// CMFC_OpenGLApp 构造

CMFC_OpenGLApp::CMFC_OpenGLApp()
{
	// 支持重新启动管理器
	m_dwRestartManagerSupportFlags = AFX_RESTART_MANAGER_SUPPORT_ALL_ASPECTS;
#ifdef _MANAGED
	// 如果应用程序是利用公共语言运行时支持(/clr)构建的，则: 
	//     1) 必须有此附加设置，“重新启动管理器”支持才能正常工作。
	//     2) 在您的项目中，您必须按照生成顺序向 System.Windows.Forms 添加引用。
	System::Windows::Forms::Application::SetUnhandledExceptionMode(System::Windows::Forms::UnhandledExceptionMode::ThrowException);
#endif

	// TODO:  将以下应用程序 ID 字符串替换为唯一的 ID 字符串；建议的字符串格式
	//为 CompanyName.ProductName.SubProduct.VersionInformation
	SetAppID(_T("MFC_OpenGL.AppID.NoVersion"));

	// TODO:  在此处添加构造代码，
	// 将所有重要的初始化放置在 InitInstance 中
}

// 唯一的一个 CMFC_OpenGLApp 对象

CMFC_OpenGLApp theApp;


// CMFC_OpenGLApp 初始化

BOOL CMFC_OpenGLApp::InitInstance()
{
	// 如果一个运行在 Windows XP 上的应用程序清单指定要
	// 使用 ComCtl32.dll 版本 6 或更高版本来启用可视化方式，
	//则需要 InitCommonControlsEx()。  否则，将无法创建窗口。
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// 将它设置为包括所有要在应用程序中使用的
	// 公共控件类。
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();


	// 初始化 OLE 库
	if (!AfxOleInit())
	{
		AfxMessageBox(IDP_OLE_INIT_FAILED);
		return FALSE;
	}

	AfxEnableControlContainer();

	EnableTaskbarInteraction(FALSE);

	// 使用 RichEdit 控件需要  AfxInitRichEdit2()	
	// AfxInitRichEdit2();

	// 标准初始化
	// 如果未使用这些功能并希望减小
	// 最终可执行文件的大小，则应移除下列
	// 不需要的特定初始化例程
	// 更改用于存储设置的注册表项
	// TODO:  应适当修改该字符串，
	// 例如修改为公司或组织名
	SetRegistryKey(_T("应用程序向导生成的本地应用程序"));
	LoadStdProfileSettings(4);  // 加载标准 INI 文件选项(包括 MRU)


	// 注册应用程序的文档模板。  文档模板
	// 将用作文档、框架窗口和视图之间的连接
	CMultiDocTemplate* pDocTemplate;
	pDocTemplate = new CMultiDocTemplate(IDR_MFC_OpenGLTYPE,
		RUNTIME_CLASS(CMFC_OpenGLDoc),
		RUNTIME_CLASS(CChildFrame), // 自定义 MDI 子框架
		RUNTIME_CLASS(CMFC_OpenGLView));
	if (!pDocTemplate)
		return FALSE;
	AddDocTemplate(pDocTemplate);

	// 创建主 MDI 框架窗口
	CMainFrame* pMainFrame = new CMainFrame;
	if (!pMainFrame || !pMainFrame->LoadFrame(IDR_MAINFRAME))
	{
		delete pMainFrame;
		return FALSE;
	}
	m_pMainWnd = pMainFrame;

	m_pMainWnd->DragAcceptFiles();
	
	// 分析标准 shell 命令、DDE、打开文件操作的命令行
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);



	// 调度在命令行中指定的命令。  如果
	// 用 /RegServer、/Register、/Unregserver 或 /Unregister 启动应用程序，则返回 FALSE。
	if (!ProcessShellCommand(cmdInfo))
		return FALSE;


	// 主窗口已初始化，因此显示它并对其进行更新
	pMainFrame->ShowWindow(m_nCmdShow);
	pMainFrame->UpdateWindow();

	

	//初始化显示帮助
	OnAppAbout();

	return TRUE;
}

int CMFC_OpenGLApp::ExitInstance()
{
	//TODO:  处理可能已添加的附加资源
	AfxOleTerm(FALSE);

	return CWinApp::ExitInstance();
}

// CMFC_OpenGLApp 消息处理程序


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };
	CFont font;
	RECT   m_source_pRectLink;
	RECT   m_blog_pRectLink;

	int cursor;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
	CStatic btn_source;
	CStatic btn_blog;
	afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	virtual BOOL OnInitDialog();
//	afx_msg void OnStnClickedStaticSource();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	afx_msg void OnNMThemeChangedStaticSource(NMHDR *pNMHDR, LRESULT *pResult);
//	afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
//	afx_msg void OnMouseHover(UINT nFlags, CPoint point);
//	afx_msg void OnDropFiles(HDROP hDropInfo);
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
	//GetDlgItem(IDC_STAC_TITLE)
	cursor = 0;
	SkinH_Attach(); //加载皮肤
	SkinH_SetAero(1); //透明
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_STATIC_SOURCE, btn_source);
	DDX_Control(pDX, IDC_STATIC_BLOG, btn_blog);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
	ON_WM_SHOWWINDOW()
	ON_WM_CTLCOLOR()
	ON_WM_MOUSEMOVE()
//	ON_STN_CLICKED(IDC_STATIC_SOURCE, &CAboutDlg::OnStnClickedStaticSource)
ON_WM_LBUTTONDOWN()
ON_NOTIFY(NM_THEMECHANGED, IDC_STATIC_SOURCE, &CAboutDlg::OnNMThemeChangedStaticSource)
ON_WM_SETCURSOR()
ON_WM_MOUSEHOVER()
ON_WM_DROPFILES()
END_MESSAGE_MAP()

// 用于运行对话框的应用程序命令
void CMFC_OpenGLApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
	
}

// CMFC_OpenGLApp 消息处理程序





void CAboutDlg::OnShowWindow(BOOL bShow, UINT nStatus)
{
	CDialogEx::OnShowWindow(bShow, nStatus);

	// TODO:  在此处添加消息处理程序代码
	font.CreatePointFont(200, _T("黑体"), NULL);
	GetDlgItem(IDC_STATIC_TITLE)->SetFont(&font);

}


HBRUSH CAboutDlg::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialogEx::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  在此更改 DC 的任何特性

	// TODO:  如果默认的不是所需画笔，则返回另一个画笔
	if (pWnd->GetDlgCtrlID() == IDC_STATIC_SOURCE){//源代码
		pDC->SetTextColor(RGB(0, 0, 255));//用RGB宏改变颜色       
	}
	if (pWnd->GetDlgCtrlID() == IDC_STATIC_BLOG){//博客
		pDC->SetTextColor(RGB(255, 0, 255));//用RGB宏改变颜色       
	}
	return hbr;
}


void CAboutDlg::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值
	if( (point.x>m_source_pRectLink.left &&point.x<m_source_pRectLink.right&&point.y>m_source_pRectLink.top&&point.y<m_source_pRectLink.bottom)
		
		)
	{
		GetDlgItem(IDC_STATIC_SOURCE)->SetWindowTextW(_T("打开"));
	}
	else if ((point.x>m_blog_pRectLink.left &&point.x<m_blog_pRectLink.right&&point.y>m_blog_pRectLink.top&&point.y < m_blog_pRectLink.bottom)){
		GetDlgItem(IDC_STATIC_BLOG)->SetWindowTextW(_T("打开"));
	}
	else
	//下面语句用来设置默认(箭头)鼠标形状,一般鼠标移开后窗口会自动恢复默认鼠标形状，可酌情添加 
	//此处酌情添加鼠标不在静态文本区的坐标算法,本例可不加 
	{
		
		GetDlgItem(IDC_STATIC_BLOG)->SetWindowTextW(_T("博客"));
		GetDlgItem(IDC_STATIC_SOURCE)->SetWindowTextW(_T("源码"));
	}
	CDialogEx::OnMouseMove(nFlags, point);
}


BOOL CAboutDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// TODO:  在此添加额外的初始化
	GetDlgItem(IDC_STATIC_BLOG)->GetWindowRect(&m_blog_pRectLink);
	//将静态文本的屏幕坐标存放在m-pRectLink中 
	ScreenToClient(&m_blog_pRectLink);
	//将屏幕坐标转换为客户坐标

	GetDlgItem(IDC_STATIC_SOURCE)->GetWindowRect(&m_source_pRectLink);
	//将静态文本的屏幕坐标存放在m-pRectLink中 
	ScreenToClient(&m_source_pRectLink);
	//将屏幕坐标转换为客户坐标
	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常:  OCX 属性页应返回 FALSE
}


//void CAboutDlg::OnStnClickedStaticSource()
//{
//	// TODO:  在此添加控件通知处理程序代码
//	ShellExecute(0, NULL, _T("www.baidu.com"), NULL, NULL, SW_NORMAL);
//	OutputDebugString(_T("点击源码"));
//}


void CAboutDlg::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值
	if( point.x>m_source_pRectLink.left &&point.x<m_source_pRectLink.right&&point.y>m_source_pRectLink.top&&point.y<m_source_pRectLink.bottom)
	{
		if (nFlags == MK_LBUTTON)//鼠标左键按下 
		{   //为改善鼠标效果，此处加入以上变换鼠标形状的代码 
			cursor = 1;
			ShellExecute(0, NULL, _T("https://github.com/halfopen/opengl_model"), NULL, NULL, SW_NORMAL);
			//也可以添加电子邮件的链接 
		}
	}
	if (point.x>m_blog_pRectLink.left &&point.x<m_blog_pRectLink.right&&point.y>m_blog_pRectLink.top&&point.y<m_blog_pRectLink.bottom)
	{
		if (nFlags == MK_LBUTTON)//鼠标左键按下 
		{   //为改善鼠标效果，此处加入以上变换鼠标形状的代码 
			//HCURSOR cursor = AfxGetApp()->LoadCursor(IDC_ARROW);
			//SetCursor(cursor);
			cursor = 1;
			ShellExecute(0, NULL, _T("http://blog.csdn.net/half_open"), NULL, NULL, SW_NORMAL);
			//也可以添加电子邮件的链接 
		}
	}
	CDialogEx::OnLButtonDown(nFlags, point);
}


BOOL CAboutDlg::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO:  在此添加专用代码和/或调用基类

	return CDialogEx::PreCreateWindow(cs);
}


void CAboutDlg::OnNMThemeChangedStaticSource(NMHDR *pNMHDR, LRESULT *pResult)
{
	// 该功能要求使用 Windows XP 或更高版本。
	// 符号 _WIN32_WINNT 必须 >= 0x0501。
	// TODO:  在此添加控件通知处理程序代码
	*pResult = 0;
}


//BOOL CAboutDlg::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message)
//{
//	// TODO:  在此添加消息处理程序代码和/或调用默认值
//	if (cursor == 1){
//		HCURSOR cursor = AfxGetApp()->LoadCursor(IDC_HAND);
//		SetCursor(cursor);
//	}
//	else{
//		HCURSOR cursor = AfxGetApp()->LoadCursor(IDC_ARROW);
//		SetCursor(cursor);
//	}
//	return TRUE;
//	//return CDialogEx::OnSetCursor(pWnd, nHitTest, message);
//}


//void CAboutDlg::OnMouseHover(UINT nFlags, CPoint point)
//{
//	// TODO:  在此添加消息处理程序代码和/或调用默认值
//	// TODO:  在此添加消息处理程序代码和/或调用默认值
//	if ((point.x>m_source_pRectLink.left &&point.x<m_source_pRectLink.right&&point.y>m_source_pRectLink.top&&point.y<m_source_pRectLink.bottom)
//		|| (point.x>m_blog_pRectLink.left &&point.x<m_blog_pRectLink.right&&point.y>m_blog_pRectLink.top&&point.y<m_blog_pRectLink.bottom)
//		)
//	{
//		HCURSOR   hCursor;
//		hCursor = AfxGetApp()->LoadCursor(IDC_HAND);
//		//将鼠标设为小手状 
//		cursor = 1;
//		SetCursor(hCursor);
//	}
//	else
//		//下面语句用来设置默认(箭头)鼠标形状,一般鼠标移开后窗口会自动恢复默认鼠标形状，可酌情添加 
//		//此处酌情添加鼠标不在静态文本区的坐标算法,本例可不加 
//	{
//		HCURSOR   hCursor;
//		hCursor = AfxGetApp()->LoadStandardCursor(IDC_ARROW);
//		//将光标设为默认值(箭头) 
//		cursor = 0;
//		SetCursor(hCursor);
//	}
//	CDialogEx::OnMouseHover(nFlags, point);
//}


//void CAboutDlg::OnDropFiles(HDROP hDropInfo)
//{
//	// TODO:  在此添加消息处理程序代码和/或调用默认值
//	// TODO: 在此添加消息处理程序代码和/或调用默认值  
//	int DropCount = DragQueryFile(hDropInfo, -1, NULL, 0);//取得被拖动文件的数目  
//	for (int i = 0; i< DropCount; i++)
//	{
//		WCHAR wcStr[MAX_PATH];
//		CString filepath;
//		DragQueryFile(hDropInfo, i, wcStr, MAX_PATH);//获得拖曳的第i个文件的文件名 
//		filepath = wcStr;
//		
//		break;//只打开第一个
//	}
//	DragFinish(hDropInfo);  //拖放结束后,释放内存  
//	CDialogEx::OnDropFiles(hDropInfo);
//}
