
// MFC_OpenGLView.cpp : CMFC_OpenGLView 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "MFC_OpenGL.h"
#endif

#include "MFC_OpenGLDoc.h"
#include "MFC_OpenGLView.h"
#include "MeshData.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMFC_OpenGLView
IMPLEMENT_DYNCREATE(CMFC_OpenGLView, CGLEnabledView)

BEGIN_MESSAGE_MAP(CMFC_OpenGLView, CGLEnabledView)
	// 标准打印命令
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
	ON_COMMAND(ID_FILE_OPEN, &CMFC_OpenGLView::OnFileOpen)

	ON_COMMAND(ID_LINE, &CMFC_OpenGLView::OnLine)
	ON_COMMAND(ID_32771, &CMFC_OpenGLView::On32771)
	ON_WM_CREATE()

	ON_COMMAND(ID_FILL, &CMFC_OpenGLView::OnFill)
	ON_COMMAND(ID_POINT, &CMFC_OpenGLView::OnPoint)
	ON_COMMAND(ID_V_POINT, &CMFC_OpenGLView::OnVPoint)
	ON_COMMAND(ID_V_FACE, &CMFC_OpenGLView::OnVFace)
	ON_WM_KEYDOWN()
	ON_COMMAND(ID_32783, &CMFC_OpenGLView::OnShowMeshInfo)
	ON_COMMAND(ID_32785, &CMFC_OpenGLView::OnWhiteLight)
	ON_COMMAND(ID_32786, &CMFC_OpenGLView::On32786)
	ON_COMMAND(ID_FILE_SAVE, &CMFC_OpenGLView::OnFileSave)
	ON_WM_DROPFILES()

	ON_COMMAND(ID_STOPSMOOTH, &CMFC_OpenGLView::OnStopsmooth)
END_MESSAGE_MAP()

// CMFC_OpenGLView 构造/析构

CMFC_OpenGLView::CMFC_OpenGLView()
{
	mesh = NULL;//读取的面片
	last_mesh = NULL;	//备份面片
	style = 2;//线框模式			
	cull = 1;//cull face
	zbuffer = 1;//启用z缓存
	show_view = 1;//默认显示图形
	PNormal = 0;//法向类型
}

CMFC_OpenGLView::~CMFC_OpenGLView()
{
}


// CMFC_OpenGLView 打印
BOOL CMFC_OpenGLView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 默认准备
	return DoPreparePrinting(pInfo);
}

void CMFC_OpenGLView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  添加额外的打印前进行的初始化过程
}

void CMFC_OpenGLView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  添加打印后进行的清理过程
}

void CMFC_OpenGLView::OnDrawGL(CDC *pDC)
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	InitalLigt();		///初始化光照信息
	glPushMatrix();
	//glLoadIdentity();
	glTranslatef(0.0f, 0.0f, scale);		//滚轮缩放
	//printf("At:%.2f %.2f %.2f\n",r*cos(c*du),h,r*sin(c*du)); //这就是视点的坐标
	gluLookAt(r*cos(c*du), h, r*sin(c*du), 0, 0, 0, 0, 1, 0); //从视点看远点,y轴方向(0,1,0)是上方向，鼠标拖动
	glRotatef(this->rotate_x, 1.0, 0.0, 0.0);

	if (mesh)	//如果打开了模型文件
	{	cull = 0;	//允许观察内部
		DrawModel(show_view, cull,style, PNormal, zbuffer);
		//OutputDebugString(_T("Draw Model\n"));
	}
	else{		//默认显示茶壶
		glutWireTeapot(0.5f);
		//pDC->DrawText(_T("画一个茶壶"),);
		myText = "默认画一个茶壶\n请打开一个ply2文件";
		//OutputDebugString(_T("茶壶\n"));
	}
	glPopMatrix();
}
// CMFC_OpenGLView 诊断

#ifdef _DEBUG
void CMFC_OpenGLView::AssertValid() const
{
	CView::AssertValid();
}

void CMFC_OpenGLView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMFC_OpenGLDoc* CMFC_OpenGLView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMFC_OpenGLDoc)));
	return (CMFC_OpenGLDoc*)m_pDocument;
}
#endif //_DEBUG

// CMFC_OpenGLView 消息处理程序
void CMFC_OpenGLView::OnFileOpen()
{
	CString filter = _T("PLY2 File (*.PLY2)|*.ply2||");
	CString str;
	CFileDialog fileDlg(TRUE, NULL, NULL, NULL, filter, NULL);
	fileDlg.m_ofn.Flags |= OFN_FILEMUSTEXIST;
	fileDlg.m_ofn.lpstrTitle = _T("读取ply2文件");
	if (fileDlg.DoModal() == IDOK)
	{
		if (mesh)			//删除之前模型
			delete mesh;
		mesh = new MeshData();	// 实例化一个新模型
		BOOL flag = mesh->ReadModelFile(fileDlg.GetPathName());	//读取文件
		//smoother->setMeshData(mesh);
		mesh->computeFaceNormal();	//计算面法向
		mesh->computeNormal();	//计算点法向
		myText.Format(_T("打开文件成功\n点数:%d \n面片数:%d"), mesh->getVertexCount(), mesh->getFaceCount());
	};
	
	this->InvalidateRect(NULL, FALSE);//刷新界面
}

void CMFC_OpenGLView::drawMyText()
{
	CDC *pDC = GetWindowDC();
	pDC->DrawText(myText, CRect(6, 6, 160, 420), DT_WORDBREAK); ReleaseDC(pDC);
	pDC->DrawText(myInfo, CRect(160, 420, 320,840), DT_WORDBREAK); ReleaseDC(pDC);
}

// 弹出对话框
void CMFC_OpenGLView::On32771()
{
	CString info;
	surfaceDialog = new SurfaceDialog;

	if (surfaceDialog->DoModal() == IDOK)
	{
		if (mesh)
		{			
			h1 = ::CreateThread(NULL, 0, SmoothThread, this, 0, NULL); //创建线程1

			if (NULL == m_pThread)
			{
				TRACE("创建新的线程出错！\n");
				return;
			}

		}
		else{
			myText = "请先打开一个ply2文件";
		}

	}	
	this->OnDraw(this->GetDC());	//重绘界面	
}

void CMFC_OpenGLView::InitalLigt()
{
	GLfloat light_position1[4] = { -52, -16, -50, 0 };
	GLfloat light_position2[4] = { -26, -48, -50, 0 };
	GLfloat light_position3[4] = { 16, -52, -50, 0 };

	GLfloat direction1[3] = { 52, 16, 50 };
	GLfloat direction2[3] = { 26, 48, 50 };
	GLfloat direction3[3] = { -16, 52, 50 };

	GLfloat light_position4[4] = { 52, 16, 50, 0 };
	GLfloat light_position5[4] = { 26, 48, 50, 0 };
	GLfloat light_position6[4] = { -16, 52, 50, 0 };

	GLfloat direction4[3] = { -52, -16, -50 };
	GLfloat direction5[3] = { -26, -48, -50 };
	GLfloat direction6[3] = { 16, -52, -50 };

	GLfloat color1[4], color2[4], color3[4], color4[4], color5[4], color6[4];

	glClearColor(1, 1, 1, 0);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	if (color_type == 0){	//彩色灯光
		color1[0] = 1; color1[1] = 0; color1[2] = 0; color1[3] = 1;
		color2[0] = 0; color2[1] = 1; color2[2] = 0; color2[3] = 1;
		color3[0] = 0; color3[1] = 0; color3[2] = 1; color3[3] = 1;

		color4[0] = 1; color4[1] = 0; color4[2] = 0; color4[3] = 1;
		color5[0] = 0; color5[1] = 1; color5[2] = 0; color5[3] = 1;
		color6[0] = 0; color6[1] = 0; color6[2] = 1; color6[3] = 1;

		GLfloat ambient[4] = { 0.3f, 0.3f, 0.3f, 1.0f };

		GLfloat material_color[4] = { 1, 1, 1, 0.5f };
		GLfloat material_specular[4] = { 0.5f, 0.5f, 0.5f, 0.5f };
		GLfloat material_ambient[4] = { 0.0, 0.0, 0.0, 0.0 };

		glLightfv(GL_LIGHT3, GL_POSITION, light_position4);
		glLightfv(GL_LIGHT3, GL_SPOT_DIRECTION, direction4);
		glLightfv(GL_LIGHT3, GL_DIFFUSE, color4);
		glLightfv(GL_LIGHT3, GL_SPECULAR, color4);

		glLightfv(GL_LIGHT4, GL_POSITION, light_position5);
		glLightfv(GL_LIGHT4, GL_SPOT_DIRECTION, direction5);
		glLightfv(GL_LIGHT4, GL_DIFFUSE, color5);
		glLightfv(GL_LIGHT4, GL_SPECULAR, color5);

		glLightfv(GL_LIGHT5, GL_POSITION, light_position6);
		glLightfv(GL_LIGHT5, GL_SPOT_DIRECTION, direction6);
		glLightfv(GL_LIGHT5, GL_DIFFUSE, color6);
		glLightfv(GL_LIGHT5, GL_SPECULAR, color6);

		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);		

		/*
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);
		glEnable(GL_LIGHT2);
		*/
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHTING);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT3);
		glEnable(GL_LIGHT4);
		glEnable(GL_LIGHT5);

		glDisable(GL_COLOR_MATERIAL);
		return;
	}
	else{
		//白色灯光
		glDisable(GL_LIGHT3);
		glDisable(GL_LIGHT4);
		glDisable(GL_LIGHT5);
		glDisable(GL_LIGHTING);
		GLfloat m_LightPostion[4] = { 0.0f, 10.0f, 10.0f, 1.0f };

		GLfloat ambientLight[] = { 0.25f, 0.25f, 0.25f, 1.0f };
		GLfloat diffuseLight[] = { 0.5, 0.5f, 0.5f, 1.0f };
		GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };

		glEnable(GL_LIGHTING);
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
		glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
		glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
		glLightfv(GL_LIGHT0, GL_POSITION, m_LightPostion);
		glEnable(GL_LIGHT0);

		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	}
}

int CMFC_OpenGLView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CGLEnabledView::OnCreate(lpCreateStruct) == -1)
		return -1;

	surfaceDialog = new SurfaceDialog;
	return 0;
}

void CMFC_OpenGLView::OnFill()
{
	if (mesh)
	{
		style = 2;
		myText = "Fill 模式";
	}
	else{
		myText = "请先打开一个ply2文件";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::OnPoint()
{
	// TODO:  在此添加命令处理程序代码
	if (mesh)
	{
		style = 0;
		myText = "Point 模式";
	}
	else{
		myText = "请先打开一个ply2文件";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::OnLine()
{
	if (mesh)
	{
		style = 1;
		myText = "Line 模式";
	}
	else{
		myText = "请先打开一个ply2文件";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::OnVPoint()
{
	if (mesh)
	{
		PNormal = 1;
		myText = "点法向";
	}
	else{
		myText = "请先打开一个ply2文件";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::OnVFace()
{
	if (mesh)
	{
		PNormal = 0;
		myText = "面法向";
	}
	else{
		myText = "请先打开一个ply2文件";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{

	CGLEnabledView::OnKeyDown(nChar, nRepCnt, nFlags);
}

BOOL CMFC_OpenGLView::PreTranslateMessage(MSG* pMsg)
{
	if (pMsg->message == WM_KEYDOWN)  // If a keydown message
	{

		if (pMsg->wParam == _T('S'))  // 按s可以直接平滑化一次，方便快速测试
		{
			if (this->mesh){
				if (last_mesh != NULL)
					delete last_mesh;
				last_mesh = new MeshData;
				last_mesh->copyMesh(mesh);
				double init_vol = mesh->getVolume();
				mesh->generateVertexLink();
				Smooth(0.1);
				mesh->computeFaceNormal();
				mesh->computeNormal();
				this->InvalidateRect(NULL, FALSE);
			};
		}
		else if (pMsg->wParam == _T('K'))  // 按s可以直接平滑化一次，方便快速测试
		{
			
			this->rotate_x += 2.0;
			if (this->rotate_x > 360)this->rotate_x = -360;
			myInfo.Format(_T("rotate_x:%f"), rotate_x);
			OutputDebugString(myInfo);
			this->InvalidateRect(NULL, FALSE);
		}
		//MessageBox((LPCTSTR)L"Key Down In PreTanslate", (LPCTSTR)L"Message", MB_OK);
	}

	return CView::PreTranslateMessage(pMsg);
}

BOOL CMFC_OpenGLView::DrawModel(bool show, bool cull, int style, bool normalpoint, bool zbuffer)
{
	CString info;

	if (!show) return TRUE;
	if (style == 0)	//绘制点图
	{
		glBegin(GL_POINTS);
		float p[3];
		for (long point = 0; point < mesh->getVertexCount(); point++){
			mesh->getVertex(point, p);
			glVertex3f(p[0], p[1], p[2]);
		}
		glEnd();
		return TRUE;
	}

	if (cull)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);

	if (zbuffer)
		glEnable(GL_DEPTH_TEST);
	else
		glDisable(GL_DEPTH_TEST);
	if (style == 2)
		//glPolygonMode(GL_FRONT_AND_BACK ,GL_FILL );
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//	glColor3f(1.0f,0.0f,1.0f);
	
	for (long face = 0; face < mesh->getFaceCount(); face++)
	{
		float p[3][3];
		int index[3];
		mesh->getFace(face, index);
		mesh->getVertex(index[0], p[0]);
		mesh->getVertex(index[1], p[1]);
		mesh->getVertex(index[2], p[2]);

		if (!normalpoint)	//画三角面片
		{
			glNormal3f(mesh->normal_f[face][0], mesh->normal_f[face][1], mesh->normal_f[face][2]);
			glBegin(GL_TRIANGLES);
			//glVertexPointer //绑定索引数组
			//info.Format(_T("point:%d (%f %f %f)"), face, p[0][0], p[1][1], p[2][2]);
			//OutputDebugString(info);
			glVertex3f(p[0][0], p[0][1], p[0][2]);
			glVertex3f(p[1][0], p[1][1], p[1][2]);
			glVertex3f(p[2][0], p[2][1], p[2][2]);
			glEnd();
		}
		else
		{
			glBegin(GL_TRIANGLES);
			glNormal3f(mesh->normal[index[0]][0], mesh->normal[index[0]][1], mesh->normal[index[0]][2]);
			glVertex3f(p[0][0], p[0][1], p[0][2]);
			glNormal3f(mesh->normal[index[1]][0], mesh->normal[index[1]][1], mesh->normal[index[1]][2]);
			glVertex3f(p[1][0], p[1][1], p[1][2]);
			glNormal3f(mesh->normal[index[2]][0], mesh->normal[index[2]][1], mesh->normal[index[2]][2]);
			glVertex3f(p[2][0], p[2][1], p[2][2]);
			glEnd();
		}
		//CSingleLock singleLock(&m_csSyncObj, FALSE);
	}
	return TRUE;
}

void CMFC_OpenGLView::Smooth(float dt)
{
	PBCG *pbcg = new PBCG;
	int n = last_mesh->vertex_N;
	int *degree = last_mesh->degree_v;
	int size = n;
	for (int i = 0; i<n; i++)size += degree[i];
	pbcg->sa = new double[size + 2];
	pbcg->ija = new unsigned long[size + 2];
	int i, j;
	unsigned long k;
	
	int **link = last_mesh->vertex_link_v;
	float(*vertex)[3] = last_mesh->vertex;
	BOOL *isBound = last_mesh->isBound;

	double *sa = pbcg->sa;
	unsigned long *ija = pbcg->ija;

	ija[1] = n + 2;
	k = n + 1;
	for (int i = 1; i <= n; i++)//遍历所有点
	{
		int ii = i - 1;	//ii从0开始计数
		int deg = degree[ii];	//点的度:相邻点的个数
		//myText.Format(_T("deg%d=%d\n"), i, deg);
		//OutputDebugString(myText);
		if (deg == 0 || isBound[ii]){	//如果没有相邻点，或者是边界
			sa[i] = 1.0;
			ija[i + 1] = k + 1;
			continue;
		}
		int *l = link[ii];
		double *w = new double[deg];
		for (int j = 0; j<deg; j++)w[j] = 0;
		double total_cot = 0;
		double A = 0;	//总面积
		//myText.Format(_T("\npoint i : %d\n"), ii);
		//OutputDebugString(myText);
		for (j = 0; j<deg; j++){	//遍历相邻点
			int t = l[j];	// t=j, s=j+1  排序后的相邻点
			int s = l[(j + 1) % deg];	
			//myText.Format(_T("[%d]->%d "), i, j, t);	
			//OutputDebugString(myText);
			//PQ为后面两个向量之间的向量差 该点与相邻两点构成一个三角形
			float PQ1[3], PQ2[3], QQ[3];
			MeshData::VEC(PQ1, vertex[ii], vertex[s]); //相邻边1
			MeshData::VEC(PQ2, vertex[ii], vertex[t]);	//相邻边2
			MeshData::VEC(QQ, vertex[s], vertex[t]);	//对边

			//三角形的面积
			double Ai = MeshData::AREA(vertex[ii], vertex[s], vertex[t]);

			if (Ai > 0){
				A += Ai;
				//向量相乘
				double dot1 = -MeshData::DOT(PQ1, QQ);	//
				double dot2 = MeshData::DOT(PQ2, QQ);	//
				//除以面积
				double cot1 = dot1 / Ai;
				double cot2 = dot2 / Ai;

				w[j] += cot1;
				w[(j + 1) % deg] += cot2;
				total_cot += cot1 + cot2;
			}
		}
		A *= 2.0;
		if (A > 0){
			sa[i] = 1.0 + dt*total_cot / A;
			for (int k = 0; k<deg; k++)
				w[k] /= A;
		}
		else
			sa[i] = 1.0;

		int *tmp = new int[deg];	//暂时保存一份点的索引
		for (j = 0; j<deg; j++)
			tmp[j] = l[j] + 1;

		//排序索引,用的是冒泡排序 = =
		for (int s = 0; s<deg; s++){
			for (int t = s; t<deg; t++){
				if (tmp[s] > tmp[t]){
					int i_tmp = tmp[s];
					tmp[s] = tmp[t];
					tmp[t] = i_tmp;
					double w_tmp = w[s];
					w[s] = w[t];
					w[t] = w_tmp;
				}
			}
		}

		for (j = 0; j<deg; j++){
			if (w[j] == 0.0)
				continue;
			int index = tmp[j];
			k++;
			sa[k] = -dt*w[j];	//
			ija[k] = index;		//
		}
		ija[i + 1] = k + 1;		//
		delete[] tmp;
		delete[] w;
	}
	//解线性方程组 Ax = b
	int iter;
	double err;
	double *b_vertex = new double[n + 1];
	double *a_vertex = new double[n + 1];
	/*
	n:大小
	b_vertex:原矩阵
	a_vertex:求解出来的矩阵
	*/
	for (i = 0; i<n; i++)b_vertex[i + 1] = a_vertex[i + 1] = vertex[i][0];
	pbcg->linbcg(n, b_vertex, a_vertex, 1, 1e-5, 100, &iter, &err);
	for (i = 0; i<n; i++)vertex[i][0] = (float)a_vertex[i + 1];//更新x坐标

	for (i = 0; i<n; i++)b_vertex[i + 1] =a_vertex[i + 1] = vertex[i][1];
	pbcg->linbcg(n, b_vertex, a_vertex, 1, 1e-5, 100, &iter, &err);
	for (i = 0; i<n; i++)vertex[i][1] = (float)a_vertex[i + 1];//更新y坐标

	for (i = 0; i<n; i++)b_vertex[i + 1] =a_vertex[i + 1] = vertex[i][2];
	pbcg->linbcg(n, b_vertex, a_vertex, 1, 1e-5, 100, &iter, &err);
	for (i = 0; i<n; i++)vertex[i][2] = (float)a_vertex[i + 1];//更新z坐标
	delete[] b_vertex;
	delete[] a_vertex;
	delete pbcg;
}

void CMFC_OpenGLView::OnShowMeshInfo()
{
	// TODO:  在此添加命令处理程序代码
	if (mesh){
		myText.Format(_T("模型信息\n点数:%d \n面片数:%d\n体积:%.2lf"),
			mesh->getVertexCount(),
			mesh->getFaceCount(),
			mesh->getVolume()
		);
	}
	this->InvalidateRect(NULL, FALSE);
}

void CMFC_OpenGLView::OnWhiteLight()
{
	// TODO:  在此添加命令处理程序代码
	
	myInfo = "白色灯光";
	color_type = 1;
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::On32786()
{
	myInfo = "彩色灯光";
	color_type = 0;
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//重绘界面
}

void CMFC_OpenGLView::OnFileSave()
{
	if (mesh == NULL){
		myText = "请先打开一个模型";
		drawMyText();
		OutputDebugString(myText);
		this->InvalidateRect(NULL, FALSE);
		return;
	}
	myText = "准备保存模型";
	drawMyText();
	OutputDebugString(myText);
	this->InvalidateRect(NULL, FALSE);
	CString buffer;
	TCHAR szFilter[] = _T("Ply2文件(*.ply2)|*.ply2|所有文件(*.*)|*.*||");
	// 构造保存文件对话框    
	CFileDialog fileDlg(FALSE, _T("ply2"), _T("newModel"), OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, this);
	CString strFilePath;
	CStdioFile file;

	// 显示保存文件对话框    
	if (IDOK == fileDlg.DoModal())
	{
		// 如果点击了文件对话框上的“保存”按钮，则将选择的文件路径显示到编辑框里    
		strFilePath = fileDlg.GetPathName();
		file.Open(strFilePath, CFile::modeCreate | CFile::modeWrite | CFile::typeText);
		int v_count = mesh->getVertexCount();
		int f_count = mesh->getFaceCount();
		buffer.Format(_T("%d\n"), v_count);
		file.WriteString(buffer);  //写入点数
		buffer.Format(_T("%d\n"), f_count);
		file.WriteString(buffer);  //写入面片数
		float p[3];
		//写入点的坐标
		for (int i = 0; i < v_count; ++i){
			mesh->getVertex(i, p);
			buffer.Format(_T("%f %f %f\n"), p[0], p[1], p[2]);
			file.WriteString(buffer);
		}
		//写入面片索引
		int f[3];
		for (int i = 0; i < f_count; ++i){
			mesh->getFace(i, f);
			buffer.Format(_T("3 %d %d %d\n"), f[0], f[1], f[2]);
			file.WriteString(buffer);
		}
		file.Close();
		myText = "保存模型成功";
		drawMyText();
		OutputDebugString(myText);
		this->InvalidateRect(NULL, FALSE);
		//SetDlgItemText(IDC_SAVE_EDIT, strFilePath);    
	}

}

BOOL CMFC_OpenGLView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO:  在此添加专用代码和/或调用基类

	if (!CGLEnabledView::PreCreateWindow(cs)){
		return FALSE;
	}
	//cs.cx = 100;
	//cs.cy = 100;
	//cs.x = 400;
	//cs.y = 400;
	return TRUE;
}

void CMFC_OpenGLView::OnDropFiles(HDROP hDropInfo)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值

	if (mesh)			//删除之前模型
		delete mesh;
	mesh = new MeshData();	// 实例化一个新模型

	WCHAR wcStr[MAX_PATH];
	DragQueryFile(hDropInfo, 0, wcStr, MAX_PATH);//获得拖曳的第一个文件的文件名  
	CString filepath = wcStr;
	OutputDebugString(filepath);
	BOOL flag = mesh->ReadModelFile(filepath);	//读取文件
	//smoother->setMeshData(mesh);
	mesh->computeFaceNormal();
	mesh->computeNormal();
	myText.Format(_T("打开文件成功\n点数:%d \n面片数:%d"), mesh->getVertexCount(), mesh->getFaceCount());


	this->InvalidateRect(NULL, FALSE);
	//CGLEnabledView::OnDropFiles(hDropInfo);
}

DWORD WINAPI SmoothThread(LPVOID lpParameter) {

	
	CString info;
	CMFC_OpenGLView* pThis = (CMFC_OpenGLView*)lpParameter;
	pThis->stopSmooth = FALSE;
	pThis->m_Sec.Lock();
	double ori_v;
	float step = pThis->surfaceDialog->m_step_length;
	int times = pThis->surfaceDialog->m_times;
	float hold_v = pThis->surfaceDialog->m_check;
	//pThis->m_Sec.Lock();
	if (hold_v == 1)ori_v = pThis->mesh->getVolume();	//获取原始体积
	pThis->m_Sec.Unlock();
	
	while (times-- && pThis->stopSmooth != TRUE)
	{
		//mesh->Smooth(1, step);
		//mesh->Test();
		pThis->m_Sec.Lock();
		if (pThis->last_mesh != NULL)//清除上次结果
			delete pThis->last_mesh;

		pThis->last_mesh = new MeshData;
		pThis->last_mesh->copyMesh(pThis->mesh);  //保存当前mesh结果,进行运算

		pThis->m_Sec.Unlock(); //锁定临界区
		//pThis->myText.Format(_T("迭代%d / %d"), pThis->surfaceDialog->m_times - times, pThis->surfaceDialog->m_times);
		//pThis->drawMyText();

		//OutputDebugString(_T("1:over\n"));
		//mesh->generateFaceLink();
		//OutputDebugString(_T("2:over\n"));	
		//AfxBeginThread(Smooth, step);
		
		pThis->last_mesh->generateVertexLink();
		OutputDebugString(_T("平滑"));
		pThis->Smooth(step);//最耗费时间的操作不加锁
		OutputDebugString(_T("平滑结束"));
		pThis->m_Sec.Lock(); //锁定临界区
		OutputDebugString(_T("加锁计算向量"));
		pThis->mesh->copyMesh(pThis->last_mesh);  //保存当前mesh结果,进行运算

		pThis->mesh->computeFaceNormal();
		
		pThis->mesh->computeNormal();
		pThis->m_Sec.Unlock();

		//myText.Format(_T("迭代第%d次完成"), surfaceDialog->m_times - times);
		pThis->myText.Format(_T("迭代%d / %d完成"), pThis->surfaceDialog->m_times - times, pThis->surfaceDialog->m_times);
		OutputDebugString(pThis->myText);
		//drawMyText();
		//pThis->OnDraw(pThis->GetDC());	//重绘界面
		pThis->PostMessage(WM_SIZE, NULL);
		//pThis->m_Sec.Unlock();
		//Sleep(2000);
		//OnDrawGL();
		//this->InvalidateRect(NULL, FALSE);
		//OutputDebugString(_T("刷新界面\n"));
	}
	pThis->m_Sec.Lock();
	pThis->myText = "迭代完成";
	pThis->drawMyText();
	if (hold_v == 1){
		//pThis->m_Sec.Lock();
		double changed_v = pThis->mesh->getVolume();
		//this->scale *= pow((ori_v / changed_v), 1 / 3.0);
		pThis->mesh->rescale(pow((ori_v / changed_v), 1 / 3.0));
		//pThis->m_Sec.Unlock();
		pThis->myText = "恢复体积";
		//drawMyText();
		OutputDebugString(pThis->myText);
	}

	pThis->m_Sec.Unlock();

	return 0;
}

void CMFC_OpenGLView::OnStopsmooth()
{
	// TODO:  在此添加命令处理程序代码
	if (h1 != NULL){
		stopSmooth = TRUE;
		myText = "中断平滑操作";
		
	}
	else{
		myText = "当前没有进行平滑操作";
	}

	drawMyText();
}
