#pragma once
#include "afxwin.h"
#include "afxcmn.h"


// SurfaceDialog 对话框

class SurfaceDialog : public CDialogEx
{
	DECLARE_DYNAMIC(SurfaceDialog)

public:
	SurfaceDialog(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~SurfaceDialog();

// 对话框数据
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
//	afx_msg void OnEnChangeEdit1();
	// 步长
	double m_step_length;
	// 迭代次数
	int m_times;
	//选中状态
	int m_check;
	CButton btn_check;
	CEdit edit_times;
	CEdit edit_step;
	afx_msg void OnBnClickedOk();
//	CProgressCtrl CProg_1;
	afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
	virtual BOOL OnInitDialog();
};
