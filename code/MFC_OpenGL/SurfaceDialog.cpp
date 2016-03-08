// SurfaceDialog.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "MFC_OpenGL.h"
#include "SurfaceDialog.h"
#include "afxdialogex.h"


// SurfaceDialog �Ի���

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


// SurfaceDialog ��Ϣ�������


//void SurfaceDialog::OnEnChangeEdit1()
//{
//	// TODO:  ����ÿؼ��� RICHEDIT �ؼ���������
//	// ���ʹ�֪ͨ��������д CDialogEx::OnInitDialog()
//	// ���������� CRichEditCtrl().SetEventMask()��
//	// ͬʱ�� ENM_CHANGE ��־�������㵽�����С�
//
//	// TODO:  �ڴ���ӿؼ�֪ͨ����������
//}


void SurfaceDialog::OnBnClickedOk()
{
	// TODO:  �ڴ���ӿؼ�֪ͨ����������
	CDialogEx::OnOK();
	CString info;
	// ��ȡѡ��״̬
	m_check = this->btn_check.GetCheck();
	//info.Format(_T("check:%d"), m_check);
	//��ȡ��������
	this->edit_times.GetWindowTextW(info);
	_stscanf(info, _T("%d"), &m_times);
	//��ȡ����
	this->edit_step.GetWindowTextW(info);
	//OutputDebugString(info);
	_stscanf(info, _T("%lf"), &m_step_length);
	//�����ȡ������
	//info.Format(_T("check:%d m_times:%d step:%f\n"), m_check, m_times, m_step_length);
	//OutputDebugString(info);
	
}

void SurfaceDialog::OnShowWindow(BOOL bShow, UINT nStatus)
{
	CDialogEx::OnShowWindow(bShow, nStatus);

	// ���Ի���Ĭ��ֵ
	btn_check.SetCheck(1);
	edit_step.SetWindowTextW(_T("0.1"));
	edit_times.SetWindowTextW(_T("3"));
}


BOOL SurfaceDialog::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// TODO:  �ڴ���Ӷ���ĳ�ʼ��
	
	return TRUE;  // return TRUE unless you set the focus to a control
	// �쳣:  OCX ����ҳӦ���� FALSE
}
