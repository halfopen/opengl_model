
// MFC_OpenGLView.h : CMFC_OpenGLView ��Ľӿ�
//

#pragma once
#include "GLEnabledView.h"
#include "SurfaceDialog.h"
#include "MeshData.h"
#include"PBCG.h"

DWORD WINAPI SmoothThread( //�����̺߳���
	LPVOID lpParameter
	);

class CMFC_OpenGLView : public CGLEnabledView
{
protected: // �������л�����
	CMFC_OpenGLView();
	DECLARE_DYNCREATE(CMFC_OpenGLView)

// ����
public:
	CMFC_OpenGLDoc* GetDocument() const;

// ����
public:

// ��д
public:
//	virtual void OnDraw(CDC* pDC);  // ��д�Ի��Ƹ���ͼ
//	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
public://��Ϊpublic,�����̷߳���
	MeshData *mesh, *last_mesh;
	CPoint MouseDownPoint;
	SurfaceDialog *surfaceDialog;
	double X_Angle;
	double Y_Angle;

	// ����ģʽ
	int style;	
	int cull;
	int zbuffer;
	int show_view;
	int PNormal;

	CString myText;//�ϲ�����
	CString myInfo;//������ʾ����

	
// ʵ��
public:
	virtual ~CMFC_OpenGLView();
public:

	//////////////////////////////////////////////////////////////////////////////////////

	//��дCGLEnabledView��ĳ��������麯��������

	void OnDrawGL(CDC* pDC);//���ƺ���
	void drawMyText();	//��ʾ��ʾ����,���Ը���
	//////////////////////////////////////////////////////////////////////////////////////

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()
public:
	//���ļ�
	afx_msg void OnFileOpen();
	//�����ļ�
	afx_msg void OnFileSave();
	
	// �����Ի���
	afx_msg void On32771();
	//��ʼ����ʱ���½�һ���Ի���
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	//����ѡ��˵�
	afx_msg void OnLine();
	afx_msg void OnFill();	
	afx_msg void OnPoint();
	//����
	afx_msg void OnVPoint();
	afx_msg void OnVFace();

	//����
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	
	//����ģ��
	BOOL DrawModel(bool show, bool cull, int style, bool normalpoint, bool zbuffer);
	//ƽ������
	void Smooth(float dt);
	BOOL read_block = false;
	CWinThread* m_pThread;     // �̶߳���ָ��
	HANDLE h1, h2; //�����߳̾��
	BOOL stopSmooth = FALSE;
	//��ʾ��Ƭ��Ϣ
	afx_msg void OnShowMeshInfo();
	//�л�Ϊ��ɫ�ƹ�
	afx_msg void OnWhiteLight();
	//�л�Ϊ��ɫ�ƹ�
	afx_msg void On32786();
	//��ʼ���ƹ�
	void InitalLigt();
	virtual BOOL PreTranslateMessage(MSG* pMsg);
	
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	afx_msg void OnDropFiles(HDROP hDropInfo);
//	afx_msg void OnWindowCascade();
	afx_msg void OnStopsmooth();//ֹͣƽ���߳�
};

#ifndef _DEBUG  // MFC_OpenGLView.cpp �еĵ��԰汾
inline CMFC_OpenGLDoc* CMFC_OpenGLView::GetDocument() const
   { return reinterpret_cast<CMFC_OpenGLDoc*>(m_pDocument); }
#endif

