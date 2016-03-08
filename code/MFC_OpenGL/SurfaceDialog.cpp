// SurfaceDialog.cpp : 实现文件
//

#include "stdafx.h"
#include "MFC_OpenGL.h"
#include "SurfaceDialog.h"
#include "afxdialogex.h"


// SurfaceDialog 对话框

IMPLEMENT_DYNAMIC(SurfaceDialog, CDialogEx)

SurfaceDialog::SurfaceDialog(CWnd* pParent /*=NULL*/)
	: CDialogEx(SurfaceDialog::IDD, pParent)
{
#ifndef _WIN32_WCE
	EnableActiveAccessibility();
#endif

}

SurfaceDialog::~SurfaceDialog()
{
}

void SurfaceDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_CHECK1, btn_check);
	DDX_Control(pDX, IDC_EDIT1, edit_times);
	DDX_Control(pDX, IDC_EDIT2, edit_step);
}


BEGIN_MESSAGE_MAP(SurfaceDialog, CDialogEx)
	ON_BN_CLICKED(IDOK, &SurfaceDialog::OnBnClickedOk)
	ON_WM_SHOWWINDOW()
END_MESSAGE_MAP()


// SurfaceDialog 消息处理程序


//void SurfaceDialog::OnEnChangeEdit1()
//{
//	// TODO:  如果该控件是 RICHEDIT 控件，它将不
//	// 发送此通知，除非重写 CDialogEx::OnInitDialog()
//	// 函数并调用 CRichEditCtrl().SetEventMask()，
//	// 同时将 ENM_CHANGE 标志“或”运算到掩码中。
//
//	// TODO:  在此添加控件通知处理程序代码
//}


void SurfaceDialog::OnBnClickedOk()
{
	// TODO:  在此添加控件通知处理程序代码
	CDialogEx::OnOK();
	CString info;
	// 获取选中状态
	m_check = this->btn_check.GetCheck();
	//info.Format(_T("check:%d"), m_check);
	//获取迭代次数
	this->edit_times.GetWindowTextW(info);
	_stscanf(info, _T("%d"), &m_times);
	//获取步长
	this->edit_step.GetWindowTextW(info);
	//OutputDebugString(info);
	_stscanf(info, _T("%lf"), &m_step_length);
	//输出获取的数据
	//info.Format(_T("check:%d m_times:%d step:%f\n"), m_check, m_times, m_step_length);
	//OutputDebugString(info);
	
}

void SurfaceDialog::OnShowWindow(BOOL bShow, UINT nStatus)
{
	CDialogEx::OnShowWindow(bShow, nStatus);

	// 给对话框默认值
	btn_check.SetCheck(1);
	edit_step.SetWindowTextW(_T("0.1"));
	edit_times.SetWindowTextW(_T("3"));
}


BOOL SurfaceDialog::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// TODO:  在此添加额外的初始化
	
	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常:  OCX 属性页应返回 FALSE
}
