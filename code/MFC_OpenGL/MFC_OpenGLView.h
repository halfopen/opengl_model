
// MFC_OpenGLView.h : CMFC_OpenGLView 类的接口
//

#pragma once
#include "GLEnabledView.h"
#include "SurfaceDialog.h"
#include "MeshData.h"
#include"PBCG.h"

DWORD WINAPI SmoothThread( //声明线程函数
	LPVOID lpParameter
	);

class CMFC_OpenGLView : public CGLEnabledView
{
protected: // 仅从序列化创建
	CMFC_OpenGLView();
	DECLARE_DYNCREATE(CMFC_OpenGLView)

// 特性
public:
	CMFC_OpenGLDoc* GetDocument() const;

// 操作
public:

// 重写
public:
//	virtual void OnDraw(CDC* pDC);  // 重写以绘制该视图
//	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
public://改为public,让子线程访问
	MeshData *mesh, *last_mesh;
	CPoint MouseDownPoint;
	SurfaceDialog *surfaceDialog;
	double X_Angle;
	double Y_Angle;

	// 绘制模式
	int style;	
	int cull;
	int zbuffer;
	int show_view;
	int PNormal;

	CString myText;//上部文字
	CString myInfo;//下面提示文字

	
// 实现
public:
	virtual ~CMFC_OpenGLView();
public:

	//////////////////////////////////////////////////////////////////////////////////////

	//改写CGLEnabledView类的场景绘制虚函数的声明

	void OnDrawGL(CDC* pDC);//绘制函数
	void drawMyText();	//显示提示文字,来自父类
	//////////////////////////////////////////////////////////////////////////////////////

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// 生成的消息映射函数
protected:
	DECLARE_MESSAGE_MAP()
public:
	//打开文件
	afx_msg void OnFileOpen();
	//保存文件
	afx_msg void OnFileSave();
	
	// 弹出对话框
	afx_msg void On32771();
	//初始化的时候，新建一个对话框
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	//绘制选项菜单
	afx_msg void OnLine();
	afx_msg void OnFill();	
	afx_msg void OnPoint();
	//法向
	afx_msg void OnVPoint();
	afx_msg void OnVFace();

	//按键
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	
	//绘制模型
	BOOL DrawModel(bool show, bool cull, int style, bool normalpoint, bool zbuffer);
	//平滑操作
	void Smooth(float dt);
	BOOL read_block = false;
	CWinThread* m_pThread;     // 线程对象指针
	HANDLE h1, h2; //定义线程句柄
	BOOL stopSmooth = FALSE;
	//显示面片信息
	afx_msg void OnShowMeshInfo();
	//切换为白色灯光
	afx_msg void OnWhiteLight();
	//切换为彩色灯光
	afx_msg void On32786();
	//初始化灯光
	void InitalLigt();
	virtual BOOL PreTranslateMessage(MSG* pMsg);
	
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	afx_msg void OnDropFiles(HDROP hDropInfo);
//	afx_msg void OnWindowCascade();
	afx_msg void OnStopsmooth();//停止平滑线程
};

#ifndef _DEBUG  // MFC_OpenGLView.cpp 中的调试版本
inline CMFC_OpenGLDoc* CMFC_OpenGLView::GetDocument() const
   { return reinterpret_cast<CMFC_OpenGLDoc*>(m_pDocument); }
#endif

