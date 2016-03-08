
// MFC_OpenGLView.cpp : CMFC_OpenGLView ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
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
	// ��׼��ӡ����
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

// CMFC_OpenGLView ����/����

CMFC_OpenGLView::CMFC_OpenGLView()
{
	mesh = NULL;//��ȡ����Ƭ
	last_mesh = NULL;	//������Ƭ
	style = 2;//�߿�ģʽ			
	cull = 1;//cull face
	zbuffer = 1;//����z����
	show_view = 1;//Ĭ����ʾͼ��
	PNormal = 0;//��������
}

CMFC_OpenGLView::~CMFC_OpenGLView()
{
}


// CMFC_OpenGLView ��ӡ
BOOL CMFC_OpenGLView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// Ĭ��׼��
	return DoPreparePrinting(pInfo);
}

void CMFC_OpenGLView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  ��Ӷ���Ĵ�ӡǰ���еĳ�ʼ������
}

void CMFC_OpenGLView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  ��Ӵ�ӡ����е��������
}

void CMFC_OpenGLView::OnDrawGL(CDC *pDC)
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	InitalLigt();		///��ʼ��������Ϣ
	glPushMatrix();
	//glLoadIdentity();
	glTranslatef(0.0f, 0.0f, scale);		//��������
	//printf("At:%.2f %.2f %.2f\n",r*cos(c*du),h,r*sin(c*du)); //������ӵ������
	gluLookAt(r*cos(c*du), h, r*sin(c*du), 0, 0, 0, 0, 1, 0); //���ӵ㿴Զ��,y�᷽��(0,1,0)���Ϸ�������϶�
	glRotatef(this->rotate_x, 1.0, 0.0, 0.0);

	if (mesh)	//�������ģ���ļ�
	{	cull = 0;	//����۲��ڲ�
		DrawModel(show_view, cull,style, PNormal, zbuffer);
		//OutputDebugString(_T("Draw Model\n"));
	}
	else{		//Ĭ����ʾ���
		glutWireTeapot(0.5f);
		//pDC->DrawText(_T("��һ�����"),);
		myText = "Ĭ�ϻ�һ�����\n���һ��ply2�ļ�";
		//OutputDebugString(_T("���\n"));
	}
	glPopMatrix();
}
// CMFC_OpenGLView ���

#ifdef _DEBUG
void CMFC_OpenGLView::AssertValid() const
{
	CView::AssertValid();
}

void CMFC_OpenGLView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMFC_OpenGLDoc* CMFC_OpenGLView::GetDocument() const // �ǵ��԰汾��������
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMFC_OpenGLDoc)));
	return (CMFC_OpenGLDoc*)m_pDocument;
}
#endif //_DEBUG

// CMFC_OpenGLView ��Ϣ�������
void CMFC_OpenGLView::OnFileOpen()
{
	CString filter = _T("PLY2 File (*.PLY2)|*.ply2||");
	CString str;
	CFileDialog fileDlg(TRUE, NULL, NULL, NULL, filter, NULL);
	fileDlg.m_ofn.Flags |= OFN_FILEMUSTEXIST;
	fileDlg.m_ofn.lpstrTitle = _T("��ȡply2�ļ�");
	if (fileDlg.DoModal() == IDOK)
	{
		if (mesh)			//ɾ��֮ǰģ��
			delete mesh;
		mesh = new MeshData();	// ʵ����һ����ģ��
		BOOL flag = mesh->ReadModelFile(fileDlg.GetPathName());	//��ȡ�ļ�
		//smoother->setMeshData(mesh);
		mesh->computeFaceNormal();	//�����淨��
		mesh->computeNormal();	//����㷨��
		myText.Format(_T("���ļ��ɹ�\n����:%d \n��Ƭ��:%d"), mesh->getVertexCount(), mesh->getFaceCount());
	};
	
	this->InvalidateRect(NULL, FALSE);//ˢ�½���
}

void CMFC_OpenGLView::drawMyText()
{
	CDC *pDC = GetWindowDC();
	pDC->DrawText(myText, CRect(6, 6, 160, 420), DT_WORDBREAK); ReleaseDC(pDC);
	pDC->DrawText(myInfo, CRect(160, 420, 320,840), DT_WORDBREAK); ReleaseDC(pDC);
}

// �����Ի���
void CMFC_OpenGLView::On32771()
{
	CString info;
	surfaceDialog = new SurfaceDialog;

	if (surfaceDialog->DoModal() == IDOK)
	{
		if (mesh)
		{			
			h1 = ::CreateThread(NULL, 0, SmoothThread, this, 0, NULL); //�����߳�1

			if (NULL == m_pThread)
			{
				TRACE("�����µ��̳߳���\n");
				return;
			}

		}
		else{
			myText = "���ȴ�һ��ply2�ļ�";
		}

	}	
	this->OnDraw(this->GetDC());	//�ػ����	
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

	if (color_type == 0){	//��ɫ�ƹ�
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
		//��ɫ�ƹ�
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
		myText = "Fill ģʽ";
	}
	else{
		myText = "���ȴ�һ��ply2�ļ�";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::OnPoint()
{
	// TODO:  �ڴ���������������
	if (mesh)
	{
		style = 0;
		myText = "Point ģʽ";
	}
	else{
		myText = "���ȴ�һ��ply2�ļ�";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::OnLine()
{
	if (mesh)
	{
		style = 1;
		myText = "Line ģʽ";
	}
	else{
		myText = "���ȴ�һ��ply2�ļ�";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::OnVPoint()
{
	if (mesh)
	{
		PNormal = 1;
		myText = "�㷨��";
	}
	else{
		myText = "���ȴ�һ��ply2�ļ�";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::OnVFace()
{
	if (mesh)
	{
		PNormal = 0;
		myText = "�淨��";
	}
	else{
		myText = "���ȴ�һ��ply2�ļ�";
	}
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{

	CGLEnabledView::OnKeyDown(nChar, nRepCnt, nFlags);
}

BOOL CMFC_OpenGLView::PreTranslateMessage(MSG* pMsg)
{
	if (pMsg->message == WM_KEYDOWN)  // If a keydown message
	{

		if (pMsg->wParam == _T('S'))  // ��s����ֱ��ƽ����һ�Σ�������ٲ���
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
		else if (pMsg->wParam == _T('K'))  // ��s����ֱ��ƽ����һ�Σ�������ٲ���
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
	if (style == 0)	//���Ƶ�ͼ
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

		if (!normalpoint)	//��������Ƭ
		{
			glNormal3f(mesh->normal_f[face][0], mesh->normal_f[face][1], mesh->normal_f[face][2]);
			glBegin(GL_TRIANGLES);
			//glVertexPointer //����������
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
	for (int i = 1; i <= n; i++)//�������е�
	{
		int ii = i - 1;	//ii��0��ʼ����
		int deg = degree[ii];	//��Ķ�:���ڵ�ĸ���
		//myText.Format(_T("deg%d=%d\n"), i, deg);
		//OutputDebugString(myText);
		if (deg == 0 || isBound[ii]){	//���û�����ڵ㣬�����Ǳ߽�
			sa[i] = 1.0;
			ija[i + 1] = k + 1;
			continue;
		}
		int *l = link[ii];
		double *w = new double[deg];
		for (int j = 0; j<deg; j++)w[j] = 0;
		double total_cot = 0;
		double A = 0;	//�����
		//myText.Format(_T("\npoint i : %d\n"), ii);
		//OutputDebugString(myText);
		for (j = 0; j<deg; j++){	//�������ڵ�
			int t = l[j];	// t=j, s=j+1  ���������ڵ�
			int s = l[(j + 1) % deg];	
			//myText.Format(_T("[%d]->%d "), i, j, t);	
			//OutputDebugString(myText);
			//PQΪ������������֮��������� �õ����������㹹��һ��������
			float PQ1[3], PQ2[3], QQ[3];
			MeshData::VEC(PQ1, vertex[ii], vertex[s]); //���ڱ�1
			MeshData::VEC(PQ2, vertex[ii], vertex[t]);	//���ڱ�2
			MeshData::VEC(QQ, vertex[s], vertex[t]);	//�Ա�

			//�����ε����
			double Ai = MeshData::AREA(vertex[ii], vertex[s], vertex[t]);

			if (Ai > 0){
				A += Ai;
				//�������
				double dot1 = -MeshData::DOT(PQ1, QQ);	//
				double dot2 = MeshData::DOT(PQ2, QQ);	//
				//�������
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

		int *tmp = new int[deg];	//��ʱ����һ�ݵ������
		for (j = 0; j<deg; j++)
			tmp[j] = l[j] + 1;

		//��������,�õ���ð������ = =
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
	//�����Է����� Ax = b
	int iter;
	double err;
	double *b_vertex = new double[n + 1];
	double *a_vertex = new double[n + 1];
	/*
	n:��С
	b_vertex:ԭ����
	a_vertex:�������ľ���
	*/
	for (i = 0; i<n; i++)b_vertex[i + 1] = a_vertex[i + 1] = vertex[i][0];
	pbcg->linbcg(n, b_vertex, a_vertex, 1, 1e-5, 100, &iter, &err);
	for (i = 0; i<n; i++)vertex[i][0] = (float)a_vertex[i + 1];//����x����

	for (i = 0; i<n; i++)b_vertex[i + 1] =a_vertex[i + 1] = vertex[i][1];
	pbcg->linbcg(n, b_vertex, a_vertex, 1, 1e-5, 100, &iter, &err);
	for (i = 0; i<n; i++)vertex[i][1] = (float)a_vertex[i + 1];//����y����

	for (i = 0; i<n; i++)b_vertex[i + 1] =a_vertex[i + 1] = vertex[i][2];
	pbcg->linbcg(n, b_vertex, a_vertex, 1, 1e-5, 100, &iter, &err);
	for (i = 0; i<n; i++)vertex[i][2] = (float)a_vertex[i + 1];//����z����
	delete[] b_vertex;
	delete[] a_vertex;
	delete pbcg;
}

void CMFC_OpenGLView::OnShowMeshInfo()
{
	// TODO:  �ڴ���������������
	if (mesh){
		myText.Format(_T("ģ����Ϣ\n����:%d \n��Ƭ��:%d\n���:%.2lf"),
			mesh->getVertexCount(),
			mesh->getFaceCount(),
			mesh->getVolume()
		);
	}
	this->InvalidateRect(NULL, FALSE);
}

void CMFC_OpenGLView::OnWhiteLight()
{
	// TODO:  �ڴ���������������
	
	myInfo = "��ɫ�ƹ�";
	color_type = 1;
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::On32786()
{
	myInfo = "��ɫ�ƹ�";
	color_type = 0;
	//this->InvalidateRect(NULL, FALSE);
	this->OnDraw(this->GetDC());	//�ػ����
}

void CMFC_OpenGLView::OnFileSave()
{
	if (mesh == NULL){
		myText = "���ȴ�һ��ģ��";
		drawMyText();
		OutputDebugString(myText);
		this->InvalidateRect(NULL, FALSE);
		return;
	}
	myText = "׼������ģ��";
	drawMyText();
	OutputDebugString(myText);
	this->InvalidateRect(NULL, FALSE);
	CString buffer;
	TCHAR szFilter[] = _T("Ply2�ļ�(*.ply2)|*.ply2|�����ļ�(*.*)|*.*||");
	// ���챣���ļ��Ի���    
	CFileDialog fileDlg(FALSE, _T("ply2"), _T("newModel"), OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, this);
	CString strFilePath;
	CStdioFile file;

	// ��ʾ�����ļ��Ի���    
	if (IDOK == fileDlg.DoModal())
	{
		// ���������ļ��Ի����ϵġ����桱��ť����ѡ����ļ�·����ʾ���༭����    
		strFilePath = fileDlg.GetPathName();
		file.Open(strFilePath, CFile::modeCreate | CFile::modeWrite | CFile::typeText);
		int v_count = mesh->getVertexCount();
		int f_count = mesh->getFaceCount();
		buffer.Format(_T("%d\n"), v_count);
		file.WriteString(buffer);  //д�����
		buffer.Format(_T("%d\n"), f_count);
		file.WriteString(buffer);  //д����Ƭ��
		float p[3];
		//д��������
		for (int i = 0; i < v_count; ++i){
			mesh->getVertex(i, p);
			buffer.Format(_T("%f %f %f\n"), p[0], p[1], p[2]);
			file.WriteString(buffer);
		}
		//д����Ƭ����
		int f[3];
		for (int i = 0; i < f_count; ++i){
			mesh->getFace(i, f);
			buffer.Format(_T("3 %d %d %d\n"), f[0], f[1], f[2]);
			file.WriteString(buffer);
		}
		file.Close();
		myText = "����ģ�ͳɹ�";
		drawMyText();
		OutputDebugString(myText);
		this->InvalidateRect(NULL, FALSE);
		//SetDlgItemText(IDC_SAVE_EDIT, strFilePath);    
	}

}

BOOL CMFC_OpenGLView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO:  �ڴ����ר�ô����/����û���

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
	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ

	if (mesh)			//ɾ��֮ǰģ��
		delete mesh;
	mesh = new MeshData();	// ʵ����һ����ģ��

	WCHAR wcStr[MAX_PATH];
	DragQueryFile(hDropInfo, 0, wcStr, MAX_PATH);//�����ҷ�ĵ�һ���ļ����ļ���  
	CString filepath = wcStr;
	OutputDebugString(filepath);
	BOOL flag = mesh->ReadModelFile(filepath);	//��ȡ�ļ�
	//smoother->setMeshData(mesh);
	mesh->computeFaceNormal();
	mesh->computeNormal();
	myText.Format(_T("���ļ��ɹ�\n����:%d \n��Ƭ��:%d"), mesh->getVertexCount(), mesh->getFaceCount());


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
	if (hold_v == 1)ori_v = pThis->mesh->getVolume();	//��ȡԭʼ���
	pThis->m_Sec.Unlock();
	
	while (times-- && pThis->stopSmooth != TRUE)
	{
		//mesh->Smooth(1, step);
		//mesh->Test();
		pThis->m_Sec.Lock();
		if (pThis->last_mesh != NULL)//����ϴν��
			delete pThis->last_mesh;

		pThis->last_mesh = new MeshData;
		pThis->last_mesh->copyMesh(pThis->mesh);  //���浱ǰmesh���,��������

		pThis->m_Sec.Unlock(); //�����ٽ���
		//pThis->myText.Format(_T("����%d / %d"), pThis->surfaceDialog->m_times - times, pThis->surfaceDialog->m_times);
		//pThis->drawMyText();

		//OutputDebugString(_T("1:over\n"));
		//mesh->generateFaceLink();
		//OutputDebugString(_T("2:over\n"));	
		//AfxBeginThread(Smooth, step);
		
		pThis->last_mesh->generateVertexLink();
		OutputDebugString(_T("ƽ��"));
		pThis->Smooth(step);//��ķ�ʱ��Ĳ���������
		OutputDebugString(_T("ƽ������"));
		pThis->m_Sec.Lock(); //�����ٽ���
		OutputDebugString(_T("������������"));
		pThis->mesh->copyMesh(pThis->last_mesh);  //���浱ǰmesh���,��������

		pThis->mesh->computeFaceNormal();
		
		pThis->mesh->computeNormal();
		pThis->m_Sec.Unlock();

		//myText.Format(_T("������%d�����"), surfaceDialog->m_times - times);
		pThis->myText.Format(_T("����%d / %d���"), pThis->surfaceDialog->m_times - times, pThis->surfaceDialog->m_times);
		OutputDebugString(pThis->myText);
		//drawMyText();
		//pThis->OnDraw(pThis->GetDC());	//�ػ����
		pThis->PostMessage(WM_SIZE, NULL);
		//pThis->m_Sec.Unlock();
		//Sleep(2000);
		//OnDrawGL();
		//this->InvalidateRect(NULL, FALSE);
		//OutputDebugString(_T("ˢ�½���\n"));
	}
	pThis->m_Sec.Lock();
	pThis->myText = "�������";
	pThis->drawMyText();
	if (hold_v == 1){
		//pThis->m_Sec.Lock();
		double changed_v = pThis->mesh->getVolume();
		//this->scale *= pow((ori_v / changed_v), 1 / 3.0);
		pThis->mesh->rescale(pow((ori_v / changed_v), 1 / 3.0));
		//pThis->m_Sec.Unlock();
		pThis->myText = "�ָ����";
		//drawMyText();
		OutputDebugString(pThis->myText);
	}

	pThis->m_Sec.Unlock();

	return 0;
}

void CMFC_OpenGLView::OnStopsmooth()
{
	// TODO:  �ڴ���������������
	if (h1 != NULL){
		stopSmooth = TRUE;
		myText = "�ж�ƽ������";
		
	}
	else{
		myText = "��ǰû�н���ƽ������";
	}

	drawMyText();
}
