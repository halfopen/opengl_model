#pragma once
#include "afxwin.h"


class CGLEnabledView :
	public CView
{
public:
	CGLEnabledView();
	~CGLEnabledView();
private:
	CDC* m_pDC;									// Windows�豸������
	HGLRC m_hRC;                                // OpenGL��Ⱦ������
	CRect m_ClientRect;							// �ͻ����Ĵ�С
	double m_dAspectRatio;						// ���ڵı���
	HPALETTE    m_hPalette;						//��ɫ��
	
public:
	DECLARE_MESSAGE_MAP()
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
//	virtual DROPEFFECT OnDragEnter(COleDataObject* pDataObject, DWORD dwKeyState, CPoint point);
	virtual void OnDrawGL(CDC * pDC);
	virtual void drawMyText();
	virtual void OnDraw(CDC* pDC);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	BOOL SetupPixelFormat();
	BOOL InitializeOpenGL(CDC* pDC);
	void SetLogicalPalette();
protected:		
	//�����Ʊ���
	float PI = 3.1415926;
	float c = PI / 180.0f;						//���ȺͽǶ�ת������
	int du = 90, oldmy = -1, oldmx = -1;		//du���ӵ���y��ĽǶ�,opengl��Ĭ��y�����Ϸ���
	float r = 1.5f, h = 0.0f;					//r���ӵ���y��İ뾶,h���ӵ�߶ȼ���y���ϵ�����
	GLfloat rotate_x = 0.0;
	//����
	GLfloat scale = 1.0;
	// ���̿���	��Ӧx,y,z������ת
	GLfloat r_x;
	GLfloat r_y;
	GLfloat r_z;
	int color_type;
public:
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);

	CCriticalSection m_Sec; //����ȫ�ֱ���m_Sec

};

