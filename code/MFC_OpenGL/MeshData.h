#if !defined(AFX_MESHDATA_H__958F025D_682B_4971_A477_DC00A2379F28__INCLUDED_)
#define AFX_MESHDATA_H__958F025D_682B_4971_A477_DC00A2379F28__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "math.h"
#include "Node.h"
#include "SVD.h"

#define INNER 0
#define BOUNDARY 1
#define NON_MANIFOLD 2

class RList{
public:
	float p1[3], p2[3];
	double k1, k2, d1, d2;
	RList* next;
	bool vis;
};

class RPoint{
public:
	float p[3];
	double k, k2;
	bool strong;
};

class REdge{
public:
	float p1[3], p2[3];
	double k;
	REdge* next;
};

class RHead{
public:
	REdge* head;
	double strength;
	RHead* next;
};

class MeshData  
{
public:
/*
	class Edges{
		int index1;
		int index2;
		Edges* next;
	}
*/
	//The number of vertex and face (File data)
	int face_N, vertex_N;

	//Vertex coordinates
	float (*vertex)[3];
	//Face index
	int (*face)[3];

	//Vertex link data (vertex to 1-ring neighbor verticse)
	//This data is not sorted (random order)
	//Node** vertex_link;
	int** vertex_link_v;
	int** vertex_link_f;
	int* degree_v;	//µ„
	int* degree_f;	//√Ê∆¨

	//Face link data (3 faces, linked along common edge)
	int (*face_link_E)[3];

	//Face link data (# = face_ad+N, linkd at common vertex)
	int *face_link_V_N;
	int **face_link_V;

	//Is the vertex specified a index boundary?
	//value = 0 if not boundary 
	//velue = index of adjacent boundary face
	int *isBound;

	//Is the vertex specified a index ridge or ravine?
	BOOL *isRidge, *isRavine;

	//Is the face specified a index ridge or ravine?
	BOOL *ridge_tri, *ravine_tri;
	
	//Vertex normal
	float (*normal)[3];

	//Face normal
	float (*normal_f)[3];

	//Mean curvature, median, varience
	float *K, midK, varience;

	//principal curvature
	double *k_max, *k_min;
	//principal direction
	double (*t_max)[3], (*t_min)[3];

	//derivative of principal curvatures associated with curvature lines
	double *r_max, *r_min, *d_max, *d_min;
	RList *ridge_L, *ravine_L;
	RHead *ridge_S, *ravine_S;

	//ridge and ravine direction
	double (*ridge_dir)[3], (*ravine_dir)[3];

	//directinal curvature on edge (How bends the edges
	double (*k_edge)[3];

	//ridge and ravines edges
	Node **ridge_edge, **ravine_edge;

	float (*center)[3];

	//edge* edge_list;

	double (*ridge_T)[3], (*ravine_T)[3];
	BOOL *isRidge_T, *isRavine_T;

	//NPR data
	struct TRI_POINT{
		int id;
		double u, v;
	}typedef tri_point;
	tri_point *sample_point;
	int sample_N;
	tri_point **t_max_line;
	int *t_max_line_N;
	tri_point **t_min_line;
	int *t_min_line_N;

	tri_point **hatch_line;
	int *hatch_line_N;

	tri_point *ridge_point;
	int ridge_point_N;
	tri_point **ridge_line;
	int *ridge_line_N;
	
	tri_point *ravine_point;
	int ravine_point_N;
	tri_point **ravine_line;
	int *ravine_line_N;

	BOOL *is_hatch;

	//Bounding Box Data for Normalize Model Size
	float original_max[3], original_min[3];
	double normalize_scale;

	int *triangle_id;
	Node *selected_triangle;

public:
	float averageOfEdgeLength();
	int selectTriangle(float s[3], float d[3]);
	void computePrincipalHeckbert();
	void computeCenter();
	void meanCurvatureNormalWithoutA(int index, double *Hn);
	void MeshData::faceNormal(float n[3], int f);
	void undoNormalizeScale();
	void tangent(double t[3], int from, int to);
	void computeHatchTriangle(float T);
	void quickSort(int *index, float *w, int start, int end);
	void paralellTrans(double vj[3], double vi[3], int i, int j);
	void traceTmin2(tri_point seed, float length, tri_point *&line, int &point_N, float step, float T);
	void traceTmin(tri_point seed, float length, tri_point *&line, int &point_N, float step);
	void traceTmax2(tri_point seed, float length, tri_point *&line, int &point_N, float step, float T);
	void traceTmax(tri_point seed, float length, tri_point *&line, int &point_N, float step);
	void computePrincipal2();
	void worldcoor(tri_point bary, float p[3]);
	void traceRavineDir2(tri_point seed, float length, tri_point *&line, int &point_N, float step, double T);
	void traceRidgeDir2(tri_point seed, float length, tri_point *&line, int &point_N, float step, double T);
	bool barycentricVector(int id, double t[3], double &u, double &v);
	void traceRavineDir(tri_point seed, float length, tri_point *&line, int &point_N, float step, double T);
	void traceRidgeDir(struct TRI_POINT seed, float length, struct TRI_POINT*& line, int& point_N, float step, double T);
	int barycentricCoor(int id, float p[3], double &u, double &v);

	void copyMesh2(MeshData* mesh);
	void computeEdgeCurvature();
	double faceArea(int f);
	void faceCenter(double c[3], int f);
	void faceCenter(float c[3], int f);
	int sortVertexLink(int *link, int N, int *&link_v, int *&link_f, int &v, int &f);
	void deleteLinkData();
	int sortVertexLink(Node* link, int *&link_v, int *&link_f, int &v, int &f);
	void rescale(double rate);
	void allocateTagEdges();
	void copyMesh(MeshData* mesh);
	float minEdge();
	void generateFaceLinkV();
	void computeRidgeDirection();
	void deleteFace(int index);
	BOOL getConsistent1Ring(int index, int *&nei_v, int *&nei_f, int &nei_N);
	void getConsistentNRings(int index, int n, int **&nei, int *&nei_N);
	double getVolume();
	int previousFace(int face_id, int vertex_id);
	int nextVertex(int face_id, int vertex_id);
	int nextFace(int face_id, int vertex_id);
	void getConsistent1Ring(int index, int *&nei, int &nei_N);
	void computePrincipal();
	void meanCurvatureNormal(int index, double* Hn);
	void laplacian(int index, double *L);
	void computeNormal();
	void computeFaceNormal();
	void get1Ring(int index, int *&nei, int &nei_N);
	//Generate face link (face_link_E)
	void generateFaceLink();
	//Generate vertex link (vertex_link)
	void generateVertexLink();

	//Count valid # face 
	int countValidFace();
	//Count valid # vertex (eliminated vertex is not included)
	int countValidVertex();

	//Scale mesh size into 20x20x20 cube 
	void normalizeScale();

	//get & set data (These data can be acessed directly, so no need)
	void getFace(int index, int vertex_index[3]);
	void getVertex(int index, float coordinates[3]);
	void setFace(int index, int vertex_index[3]);
	void setVertex(int index, float coordinates[3]);
	int getFaceCount();
	int getVertexCount();
	void setFaceCount(int face_N);
	void setVertexCount(int vertex_N);
	BOOL MeshData::ReadModelFile(CString filePath);
	
	BOOL MeshData::Smooth(int times, float step);
	void MeshData::Test();
	MeshData();
	virtual ~MeshData();

	static inline double MeshData::LENGTH(float v[3]){
		return sqrt((double)v[0]*(double)v[0] + (double)v[1]*(double)v[1] + (double)v[2]*(double)v[2]);
	}

	static inline double MeshData::LENGTH(double v[3]){
		return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	}

	static inline void MeshData::VEC(float v[3], float p1[3], float p2[3]){
		v[0] = p2[0]-p1[0];
		v[1] = p2[1]-p1[1];
		v[2] = p2[2]-p1[2];
	}

	static inline void MeshData::VEC(double v[3], double p1[3], double p2[3]){
		v[0] = p2[0]-p1[0];
		v[1] = p2[1]-p1[1];
		v[2] = p2[2]-p1[2];
	}

	static inline void MeshData::VEC(double v[3], double p1[3], float p2[3]){
		v[0] = (double)p2[0]-p1[0];
		v[1] = (double)p2[1]-p1[1];
		v[2] = (double)p2[2]-p1[2];
	}

	static inline void MeshData::VEC(double v[3], float p1[3], double p2[3]){
		v[0] = p2[0]-(double)p1[0];
		v[1] = p2[1]-(double)p1[1];
		v[2] = p2[2]-(double)p1[2];
	}
	
	static inline void MeshData::TIMES(float kv[3], float k, float v[3]){
		kv[0] = k*v[0];
		kv[1] = k*v[1];
		kv[2] = k*v[2];
	}

	static inline double MeshData::DIST(float p1[3], float p2[3]){
		float v[3];
		VEC(v,p1,p2);
		return LENGTH(v);
	}

	static inline double MeshData::DIST(double p1[3], double p2[3]){
		double v[3];
		VEC(v,p1,p2);
		return LENGTH(v);
	}

	static inline double MeshData::DIST2(float p1[3], float p2[3]){
		float v[3];
		VEC(v,p1,p2);
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	static inline double MeshData::DIST2(double p1[3], double p2[3]){
		double v[3];
		VEC(v,p1,p2);
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	static inline double MeshData::DIST2(double p1[3], float p2[3]){
		double v[3];
		VEC(v,p1,p2);
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	static inline double MeshData::DOT(float v1[3], float v2[3]){
		return (double)v1[0]*(double)v2[0] + (double)v1[1]*(double)v2[1] + (double)v1[2]*(double)v2[2];
	}

	static inline double MeshData::DOT(double v1[3], float v2[3]){
		return v1[0]*(double)v2[0] + v1[1]*(double)v2[1] + v1[2]*(double)v2[2];
	}

	static inline double MeshData::DOT(float v1[3], double v2[3]){
		return (double)v1[0]*v2[0] + (double)v1[1]*v2[1] + (double)v1[2]*v2[2];
	}

	static inline double MeshData::DOT(double v1[3], double v2[3]){
		return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	}

	static inline void MeshData::CROSS(float n[3], float v1[3], float v2[3]){
		n[0] = v1[1]*v2[2] - v1[2]*v2[1];
		n[1] = v1[2]*v2[0] - v1[0]*v2[2];
		n[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	static inline void MeshData::CROSS(double n[3], double v1[3], double v2[3]){
		n[0] = v1[1]*v2[2] - v1[2]*v2[1];
		n[1] = v1[2]*v2[0] - v1[0]*v2[2];
		n[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	static inline void MeshData::CROSS(double n[3], float v1[3], float v2[3]){
		n[0] = (double)v1[1]*(double)v2[2] - (double)v1[2]*(double)v2[1];
		n[1] = (double)v1[2]*(double)v2[0] - (double)v1[0]*(double)v2[2];
		n[2] = (double)v1[0]*(double)v2[1] - (double)v1[1]*(double)v2[0];
	}

	static inline void MeshData::CROSS(double n[3], double v1[3], float v2[3]){
		n[0] = v1[1]*(double)v2[2] - v1[2]*(double)v2[1];
		n[1] = v1[2]*(double)v2[0] - v1[0]*(double)v2[2];
		n[2] = v1[0]*(double)v2[1] - v1[1]*(double)v2[0];
	}

	static inline void MeshData::CROSS(double n[3], float v1[3], double v2[3]){
		n[0] = (double)v1[1]*v2[2] - (double)v1[2]*v2[1];
		n[1] = (double)v1[2]*v2[0] - (double)v1[0]*v2[2];
		n[2] = (double)v1[0]*v2[1] - (double)v1[1]*v2[0];
	}

	static inline void MeshData::CROSS(float n[3], double v1[3], double v2[3]){
		n[0] = (float)(v1[1]*v2[2] - v1[2]*v2[1]);
		n[1] = (float)(v1[2]*v2[0] - v1[0]*v2[2]);
		n[2] = (float)(v1[0]*v2[1] - v1[1]*v2[0]);
	}

	static inline void MeshData::CROSS(float n[3], double v1[3], float v2[3]){
		n[0] = (float)(v1[1]*v2[2] - v1[2]*v2[1]);
		n[1] = (float)(v1[2]*v2[0] - v1[0]*v2[2]);
		n[2] = (float)(v1[0]*v2[1] - v1[1]*v2[0]);
	}

	static inline void MeshData::CROSS(float n[3], float v1[3], double v2[3]){
		n[0] = (float)(v1[1]*v2[2] - v1[2]*v2[1]);
		n[1] = (float)(v1[2]*v2[0] - v1[0]*v2[2]);
		n[2] = (float)(v1[0]*v2[1] - v1[1]*v2[0]);
	}

	static inline double MeshData::AREA(float p1[3], float p2[3], float p3[3])
	{
		double n[3];
		float v1[3], v2[3];
		VEC(v1, p2, p1);
		VEC(v2, p3, p1);
		CROSS(n, v1, v2);

		return 0.5*LENGTH(n);
	}

	static inline double MeshData::AREA(double p1[3], double p2[3], double p3[3])
	{
		double n[3];
		double v1[3], v2[3];
		VEC(v1, p2, p1);
		VEC(v2, p3, p1);
		CROSS(n, v1, v2);

		return 0.5*LENGTH(n);
	}

	static inline double MeshData::AREA(float p1[3], float p2[3], double p3[3])
	{
		double n[3];
		float v1[3];
		double v2[3];
		VEC(v1, p2, p1);
		VEC(v2, p3, p1);
		CROSS(n, v1, v2);

		return 0.5*LENGTH(n);
	}

		static inline void MAT_VEC(double y[3], double A[3][3], double x[3]){
		y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
		y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
		y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
	}

	static inline void MAT_TIMES(double C[3][3], double A[3][3], double B[3][3]){
		C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
		C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
		C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];

		C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
		C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
		C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];

		C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
		C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
		C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
	}

	static inline void GENERATE_MAT(double R[3][3], double angle, double v[3]){
		double l = sqrt(v[0]*v[0] + v[1]*v[1]);
		double a;
		if(l == 0)
			a = 0;
		else
			a = acos(v[1]/l);
		if(v[0] < 0)
			a = -a;
		double b = acos(v[2]/MeshData::LENGTH(v));

		double Rz[3][3], Rx[3][3], Rt[3][3]; 
		Rz[0][0] = Rz[1][1] = (float)cos(a);
		Rz[0][1] = sin(a);
		Rz[1][0] = -Rz[0][1];
		Rz[2][0] = Rz[2][1] = Rz[0][2] = Rz[1][2] = 0;
		Rz[2][2] = 1;

		Rt[0][0] = Rt[1][1] = cos(angle);
		Rt[0][1] = -sin(angle);
		Rt[1][0] = -Rt[0][1];
		Rt[2][0] = Rt[2][1] = Rt[0][2] = Rt[1][2] = 0;
		Rt[2][2] = 1;

		Rx[0][0] = 1;
		Rx[0][1] = Rx[0][2] = Rx[1][0] = Rx[2][0] = 0;
		Rx[1][1] = Rx[2][2] = cos(b);
		Rx[1][2] = sin(b);
		Rx[2][1] = -Rx[1][2];

		double tmp[3][3];
		MAT_TIMES(tmp, Rz, Rx);
		MAT_TIMES(R, tmp, Rt);
		Rz[1][0] *= -1;
		Rz[0][1] *= -1;
		Rx[1][2] *= -1;
		Rx[2][1] *= -1;
		MAT_TIMES(tmp, R, Rx);
		MAT_TIMES(R, tmp, Rz);

		MAT_TIMES(tmp, Rx, Rz);
		double x[3];
		MAT_VEC(x, tmp, v);
	}

	static inline void GENERATE_MAT(double R[3][3], double angle, float v1[3]){
		double v[3];
		v[0] = v1[0];
		v[1] = v1[1];
		v[2] = v1[2];

		double l = sqrt(v[0]*v[0] + v[1]*v[1]);
		double a;
		if(l == 0)
			a = 0;
		else
			a = acos(v[1]/l);
		if(v[0] < 0)
			a = -a;
		double b = acos(v[2]/MeshData::LENGTH(v));

		double Rz[3][3], Rx[3][3], Rt[3][3]; 
		Rz[0][0] = Rz[1][1] = (float)cos(a);
		Rz[0][1] = sin(a);
		Rz[1][0] = -Rz[0][1];
		Rz[2][0] = Rz[2][1] = Rz[0][2] = Rz[1][2] = 0;
		Rz[2][2] = 1;

		Rt[0][0] = Rt[1][1] = cos(angle);
		Rt[0][1] = -sin(angle);
		Rt[1][0] = -Rt[0][1];
		Rt[2][0] = Rt[2][1] = Rt[0][2] = Rt[1][2] = 0;
		Rt[2][2] = 1;

		Rx[0][0] = 1;
		Rx[0][1] = Rx[0][2] = Rx[1][0] = Rx[2][0] = 0;
		Rx[1][1] = Rx[2][2] = cos(b);
		Rx[1][2] = sin(b);
		Rx[2][1] = -Rx[1][2];

		double tmp[3][3];
		MAT_TIMES(tmp, Rz, Rx);
		MAT_TIMES(R, tmp, Rt);
		Rz[1][0] *= -1;
		Rz[0][1] *= -1;
		Rx[1][2] *= -1;
		Rx[2][1] *= -1;
		MAT_TIMES(tmp, R, Rx);
		MAT_TIMES(R, tmp, Rz);

		MAT_TIMES(tmp, Rx, Rz);
		double x[3];
		MAT_VEC(x, tmp, v);
	}
};

#endif // !defined(AFX_MESHDATA_H__958F025D_682B_4971_A477_DC00A2379F28__INCLUDED_)
