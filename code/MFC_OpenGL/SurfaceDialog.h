#pragma once
#include "afxwin.h"
#include "afxcmn.h"


// SurfaceDialog �Ի���

class SurfaceDialog : public CDialogEx
{
	DECLARE_DYNAMIC(SurfaceDialog)

public:
	SurfaceDialog(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~SurfaceDialog();

// �Ի�������
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

	DECLARE_MESSAGE_MAP()
public:
//	afx_msg void OnEnChangeEdit1();
	// ����
	double m_step_length;
	// ��������
	int m_times;
	//ѡ��״̬
	int m_check;
	CButton btn_check;
	CEdit edit_times;
	CEdit edit_step;
	afx_msg void OnBnClickedOk();
//	CProgressCtrl CProg_1;
	afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
	virtual BOOL OnInitDialog();
};
