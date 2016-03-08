#pragma once
#include "afxwin.h"


class CGLEnabledView :
	public CView
{
public:
	CGLEnabledView();
	~CGLEnabledView();
private:
	CDC* m_pDC;									// Windows设备描述表
	HGLRC m_hRC;                                // OpenGL渲染描述表
	CRect m_ClientRect;							// 客户区的大小
	double m_dAspectRatio;						// 窗口的比例
	HPALETTE    m_hPalette;						//调色板
	
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
	//鼠标控制变量
	float PI = 3.1415926;
	float c = PI / 180.0f;						//弧度和角度转换参数
	int du = 90, oldmy = -1, oldmx = -1;		//du是视点绕y轴的角度,opengl里默认y轴是上方向
	float r = 1.5f, h = 0.0f;					//r是视点绕y轴的半径,h是视点高度即在y轴上的坐标
	GLfloat rotate_x = 0.0;
	//缩放
	GLfloat scale = 1.0;
	// 键盘控制	对应x,y,z方向旋转
	GLfloat r_x;
	GLfloat r_y;
	GLfloat r_z;
	int color_type;
public:
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);

	CCriticalSection m_Sec; //定义全局变量m_Sec

};

