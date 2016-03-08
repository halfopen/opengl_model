
// MFC_OpenGL.cpp : ����Ӧ�ó��������Ϊ��
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
	// �����ļ��ı�׼�ĵ�����
	ON_COMMAND(ID_FILE_NEW, &CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, &CWinApp::OnFileOpen)
	// ��׼��ӡ��������
	ON_COMMAND(ID_FILE_PRINT_SETUP, &CWinApp::OnFilePrintSetup)
END_MESSAGE_MAP()


// CMFC_OpenGLApp ����

CMFC_OpenGLApp::CMFC_OpenGLApp()
{
	// ֧����������������
	m_dwRestartManagerSupportFlags = AFX_RESTART_MANAGER_SUPPORT_ALL_ASPECTS;
#ifdef _MANAGED
	// ���Ӧ�ó��������ù�����������ʱ֧��(/clr)�����ģ���: 
	//     1) �����д˸������ã�������������������֧�ֲ�������������
	//     2) ��������Ŀ�У������밴������˳���� System.Windows.Forms ������á�
	System::Windows::Forms::Application::SetUnhandledExceptionMode(System::Windows::Forms::UnhandledExceptionMode::ThrowException);
#endif

	// TODO:  ������Ӧ�ó��� ID �ַ����滻ΪΨһ�� ID �ַ�����������ַ�����ʽ
	//Ϊ CompanyName.ProductName.SubProduct.VersionInformation
	SetAppID(_T("MFC_OpenGL.AppID.NoVersion"));

	// TODO:  �ڴ˴���ӹ�����룬
	// ��������Ҫ�ĳ�ʼ�������� InitInstance ��
}

// Ψһ��һ�� CMFC_OpenGLApp ����

CMFC_OpenGLApp theApp;


// CMFC_OpenGLApp ��ʼ��

BOOL CMFC_OpenGLApp::InitInstance()
{
	// ���һ�������� Windows XP �ϵ�Ӧ�ó����嵥ָ��Ҫ
	// ʹ�� ComCtl32.dll �汾 6 ����߰汾�����ÿ��ӻ���ʽ��
	//����Ҫ InitCommonControlsEx()��  ���򣬽��޷��������ڡ�
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// ��������Ϊ��������Ҫ��Ӧ�ó�����ʹ�õ�
	// �����ؼ��ࡣ
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();


	// ��ʼ�� OLE ��
	if (!AfxOleInit())
	{
		AfxMessageBox(IDP_OLE_INIT_FAILED);
		return FALSE;
	}

	AfxEnableControlContainer();

	EnableTaskbarInteraction(FALSE);

	// ʹ�� RichEdit �ؼ���Ҫ  AfxInitRichEdit2()	
	// AfxInitRichEdit2();

	// ��׼��ʼ��
	// ���δʹ����Щ���ܲ�ϣ����С
	// ���տ�ִ���ļ��Ĵ�С����Ӧ�Ƴ�����
	// ����Ҫ���ض���ʼ������
	// �������ڴ洢���õ�ע�����
	// TODO:  Ӧ�ʵ��޸ĸ��ַ�����
	// �����޸�Ϊ��˾����֯��
	SetRegistryKey(_T("Ӧ�ó��������ɵı���Ӧ�ó���"));
	LoadStdProfileSettings(4);  // ���ر�׼ INI �ļ�ѡ��(���� MRU)


	// ע��Ӧ�ó�����ĵ�ģ�塣  �ĵ�ģ��
	// �������ĵ�����ܴ��ں���ͼ֮�������
	CMultiDocTemplate* pDocTemplate;
	pDocTemplate = new CMultiDocTemplate(IDR_MFC_OpenGLTYPE,
		RUNTIME_CLASS(CMFC_OpenGLDoc),
		RUNTIME_CLASS(CChildFrame), // �Զ��� MDI �ӿ��
		RUNTIME_CLASS(CMFC_OpenGLView));
	if (!pDocTemplate)
		return FALSE;
	AddDocTemplate(pDocTemplate);

	// ������ MDI ��ܴ���
	CMainFrame* pMainFrame = new CMainFrame;
	if (!pMainFrame || !pMainFrame->LoadFrame(IDR_MAINFRAME))
	{
		delete pMainFrame;
		return FALSE;
	}
	m_pMainWnd = pMainFrame;

	m_pMainWnd->DragAcceptFiles();
	
	// ������׼ shell ���DDE�����ļ�������������
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);



	// ��������������ָ�������  ���
	// �� /RegServer��/Register��/Unregserver �� /Unregister ����Ӧ�ó����򷵻� FALSE��
	if (!ProcessShellCommand(cmdInfo))
		return FALSE;


	// �������ѳ�ʼ���������ʾ����������и���
	pMainFrame->ShowWindow(m_nCmdShow);
	pMainFrame->UpdateWindow();

	

	//��ʼ����ʾ����
	OnAppAbout();

	return TRUE;
}

int CMFC_OpenGLApp::ExitInstance()
{
	//TODO:  �����������ӵĸ�����Դ
	AfxOleTerm(FALSE);

	return CWinApp::ExitInstance();
}

// CMFC_OpenGLApp ��Ϣ�������


// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// �Ի�������
	enum { IDD = IDD_ABOUTBOX };
	CFont font;
	RECT   m_source_pRectLink;
	RECT   m_blog_pRectLink;

	int cursor;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
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
	SkinH_Attach(); //����Ƥ��
	SkinH_SetAero(1); //͸��
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

// �������жԻ����Ӧ�ó�������
void CMFC_OpenGLApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
	
}

// CMFC_OpenGLApp ��Ϣ�������





void CAboutDlg::OnShowWindow(BOOL bShow, UINT nStatus)
{
	CDialogEx::OnShowWindow(bShow, nStatus);

	// TODO:  �ڴ˴������Ϣ����������
	font.CreatePointFont(200, _T("����"), NULL);
	GetDlgItem(IDC_STATIC_TITLE)->SetFont(&font);

}


HBRUSH CAboutDlg::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialogEx::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  �ڴ˸��� DC ���κ�����

	// TODO:  ���Ĭ�ϵĲ������軭�ʣ��򷵻���һ������
	if (pWnd->GetDlgCtrlID() == IDC_STATIC_SOURCE){//Դ����
		pDC->SetTextColor(RGB(0, 0, 255));//��RGB��ı���ɫ       
	}
	if (pWnd->GetDlgCtrlID() == IDC_STATIC_BLOG){//����
		pDC->SetTextColor(RGB(255, 0, 255));//��RGB��ı���ɫ       
	}
	return hbr;
}


void CAboutDlg::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ
	if( (point.x>m_source_pRectLink.left &&point.x<m_source_pRectLink.right&&point.y>m_source_pRectLink.top&&point.y<m_source_pRectLink.bottom)
		
		)
	{
		GetDlgItem(IDC_STATIC_SOURCE)->SetWindowTextW(_T("��"));
	}
	else if ((point.x>m_blog_pRectLink.left &&point.x<m_blog_pRectLink.right&&point.y>m_blog_pRectLink.top&&point.y < m_blog_pRectLink.bottom)){
		GetDlgItem(IDC_STATIC_BLOG)->SetWindowTextW(_T("��"));
	}
	else
	//���������������Ĭ��(��ͷ)�����״,һ������ƿ��󴰿ڻ��Զ��ָ�Ĭ�������״����������� 
	//�˴����������겻�ھ�̬�ı����������㷨,�����ɲ��� 
	{
		
		GetDlgItem(IDC_STATIC_BLOG)->SetWindowTextW(_T("����"));
		GetDlgItem(IDC_STATIC_SOURCE)->SetWindowTextW(_T("Դ��"));
	}
	CDialogEx::OnMouseMove(nFlags, point);
}


BOOL CAboutDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// TODO:  �ڴ���Ӷ���ĳ�ʼ��
	GetDlgItem(IDC_STATIC_BLOG)->GetWindowRect(&m_blog_pRectLink);
	//����̬�ı�����Ļ��������m-pRectLink�� 
	ScreenToClient(&m_blog_pRectLink);
	//����Ļ����ת��Ϊ�ͻ�����

	GetDlgItem(IDC_STATIC_SOURCE)->GetWindowRect(&m_source_pRectLink);
	//����̬�ı�����Ļ��������m-pRectLink�� 
	ScreenToClient(&m_source_pRectLink);
	//����Ļ����ת��Ϊ�ͻ�����
	return TRUE;  // return TRUE unless you set the focus to a control
	// �쳣:  OCX ����ҳӦ���� FALSE
}


//void CAboutDlg::OnStnClickedStaticSource()
//{
//	// TODO:  �ڴ���ӿؼ�֪ͨ����������
//	ShellExecute(0, NULL, _T("www.baidu.com"), NULL, NULL, SW_NORMAL);
//	OutputDebugString(_T("���Դ��"));
//}


void CAboutDlg::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ
	if( point.x>m_source_pRectLink.left &&point.x<m_source_pRectLink.right&&point.y>m_source_pRectLink.top&&point.y<m_source_pRectLink.bottom)
	{
		if (nFlags == MK_LBUTTON)//���������� 
		{   //Ϊ�������Ч�����˴��������ϱ任�����״�Ĵ��� 
			cursor = 1;
			ShellExecute(0, NULL, _T("https://github.com/halfopen/opengl_model"), NULL, NULL, SW_NORMAL);
			//Ҳ������ӵ����ʼ������� 
		}
	}
	if (point.x>m_blog_pRectLink.left &&point.x<m_blog_pRectLink.right&&point.y>m_blog_pRectLink.top&&point.y<m_blog_pRectLink.bottom)
	{
		if (nFlags == MK_LBUTTON)//���������� 
		{   //Ϊ�������Ч�����˴��������ϱ任�����״�Ĵ��� 
			//HCURSOR cursor = AfxGetApp()->LoadCursor(IDC_ARROW);
			//SetCursor(cursor);
			cursor = 1;
			ShellExecute(0, NULL, _T("http://blog.csdn.net/half_open"), NULL, NULL, SW_NORMAL);
			//Ҳ������ӵ����ʼ������� 
		}
	}
	CDialogEx::OnLButtonDown(nFlags, point);
}


BOOL CAboutDlg::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO:  �ڴ����ר�ô����/����û���

	return CDialogEx::PreCreateWindow(cs);
}


void CAboutDlg::OnNMThemeChangedStaticSource(NMHDR *pNMHDR, LRESULT *pResult)
{
	// �ù���Ҫ��ʹ�� Windows XP ����߰汾��
	// ���� _WIN32_WINNT ���� >= 0x0501��
	// TODO:  �ڴ���ӿؼ�֪ͨ����������
	*pResult = 0;
}


//BOOL CAboutDlg::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message)
//{
//	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ
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
//	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ
//	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ
//	if ((point.x>m_source_pRectLink.left &&point.x<m_source_pRectLink.right&&point.y>m_source_pRectLink.top&&point.y<m_source_pRectLink.bottom)
//		|| (point.x>m_blog_pRectLink.left &&point.x<m_blog_pRectLink.right&&point.y>m_blog_pRectLink.top&&point.y<m_blog_pRectLink.bottom)
//		)
//	{
//		HCURSOR   hCursor;
//		hCursor = AfxGetApp()->LoadCursor(IDC_HAND);
//		//�������ΪС��״ 
//		cursor = 1;
//		SetCursor(hCursor);
//	}
//	else
//		//���������������Ĭ��(��ͷ)�����״,һ������ƿ��󴰿ڻ��Զ��ָ�Ĭ�������״����������� 
//		//�˴����������겻�ھ�̬�ı����������㷨,�����ɲ��� 
//	{
//		HCURSOR   hCursor;
//		hCursor = AfxGetApp()->LoadStandardCursor(IDC_ARROW);
//		//�������ΪĬ��ֵ(��ͷ) 
//		cursor = 0;
//		SetCursor(hCursor);
//	}
//	CDialogEx::OnMouseHover(nFlags, point);
//}


//void CAboutDlg::OnDropFiles(HDROP hDropInfo)
//{
//	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ
//	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ  
//	int DropCount = DragQueryFile(hDropInfo, -1, NULL, 0);//ȡ�ñ��϶��ļ�����Ŀ  
//	for (int i = 0; i< DropCount; i++)
//	{
//		WCHAR wcStr[MAX_PATH];
//		CString filepath;
//		DragQueryFile(hDropInfo, i, wcStr, MAX_PATH);//�����ҷ�ĵ�i���ļ����ļ��� 
//		filepath = wcStr;
//		
//		break;//ֻ�򿪵�һ��
//	}
//	DragFinish(hDropInfo);  //�ϷŽ�����,�ͷ��ڴ�  
//	CDialogEx::OnDropFiles(hDropInfo);
//}
