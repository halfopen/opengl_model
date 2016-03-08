#include "stdafx.h"
#include "MeshData.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define PI 3.14159265
int i, j, k, n;

MeshData::MeshData()
{
	face_N = 0;
	vertex_N = 0;

	vertex = NULL;
	face = NULL;

	vertex_link_v = NULL;
	vertex_link_f = NULL;

	degree_v = NULL;
	degree_f = NULL;

	face_link_E = NULL;

	face_link_V_N = NULL;
	face_link_V = NULL;

	isBound = NULL;

	isRidge = NULL;
	isRavine = NULL;

	ridge_tri = NULL;
	ravine_tri = NULL;

	ridge_edge = NULL;
	ravine_edge = NULL;
	
	normal = NULL;
	normal_f = NULL;

	K = NULL;
	midK = 0;
	varience = 0;

	k_max = NULL;
	k_min = NULL;
	
	t_max = NULL;
	t_min = NULL;

	k_edge = NULL;

	ridge_dir = NULL;
	ravine_dir = NULL;

	//edge_list = NULL;

	isRidge_T = NULL;
	isRavine_T = NULL;

	ridge_T = NULL;
	ravine_T = NULL;

	sample_point = NULL;
	sample_N = 0;

	center = NULL;

	triangle_id = NULL;
	selected_triangle = NULL;

	ridge_S = NULL;
	ravine_S = NULL;
}

MeshData::~MeshData()
{
	this->deleteLinkData();
}

void MeshData::setVertexCount(int vertex_N)
{
	if(vertex != NULL)
		delete[] vertex;
	this->vertex_N = vertex_N;
	vertex = new float[vertex_N][3];
}

void MeshData::setFaceCount(int face_N)
{
	if(face != NULL)
		delete[] face;
	this->face_N = face_N;
	face = (int (*)[3])new int[face_N][3];
}

int MeshData::getVertexCount()
{
	return vertex_N;
}

int MeshData::getFaceCount()
{
	return face_N;
}

void MeshData::setVertex(int index, float coordinates[3])
{
	vertex[index][0] = coordinates[0];
	vertex[index][1] = coordinates[1];
	vertex[index][2] = coordinates[2];
}

void MeshData::setFace(int index, int vertex_index[])
{
	face[index][0] = vertex_index[0];
	face[index][1] = vertex_index[1];
	face[index][2] = vertex_index[2];
}

void MeshData::getVertex(int index, float coordinates[3])
{
	coordinates[0] = vertex[index][0];
	coordinates[1] = vertex[index][1];
	coordinates[2] = vertex[index][2];
}

void MeshData::getFace(int index, int vertex_index[])
{
	vertex_index[0] = face[index][0];
	vertex_index[1] = face[index][1];
	vertex_index[2] = face[index][2];
}

void MeshData::normalizeScale()
{
	float max[3], min[3];
	max[0] = max[1] = max[2] = -1000000;
	min[0] = min[1] = min[2] = 1000000;
	for(int i=0; i<vertex_N; i++){
		for(int j=0; j<3; j++){
			if(max[j] < vertex[i][j])
				max[j] = vertex[i][j];
			else if(min[j] > vertex[i][j])
				min[j] = vertex[i][j];
		}
	}
	float sizeX = max[0] - min[0];
	float sizeY = max[1] - min[1];
	float sizeZ = max[2] - min[2];

	double scale = sizeX;
	if(scale < sizeY)
		scale = sizeY;
	if(scale < sizeZ)
		scale = sizeZ;
	normalize_scale = scale*0.05;
	scale = 20.0/scale; //20.0f/scale;

	for(i=0; i<vertex_N; i++){
		vertex[i][0] = (vertex[i][0]-0.5f*(max[0]+min[0]))*scale;
		vertex[i][1] = (vertex[i][1]-0.5f*(max[1]+min[1]))*scale;
		vertex[i][2] = (vertex[i][2]-0.5f*(max[2]+min[2]))*scale;// + 5.0f; //10.0f;
	}

	original_max[0] = max[0];
	original_max[1] = max[1];
	original_max[2] = max[2];

	original_min[0] = min[0];
	original_min[1] = min[1];
	original_min[2] = min[2];

}

int MeshData::countValidVertex()
{
	int counter = 0;
	for(int i=0; i<vertex_N; i++)
		if(degree_v[i] != 0)
			counter++;
	return counter;
}

int MeshData::countValidFace()
{
	int counter = 0;
	for(int i=0; i<face_N; i++)
		if(face[i][0] != -1)
			counter++;
	return counter;
}

void MeshData::generateVertexLink()
{
	int* link_N = new int[vertex_N];
	for(int i=0; i<vertex_N; i++)
		link_N[i] = 0;
	for(i=0; i<face_N; i++){
		for(int j=0; j<3; j++){
			if(face[i][j] < 0 || face[i][j] >= vertex_N)
				face[i][0] = face[i][1] = face[i][2] = -1;
			if(face[i][j] == face[i][(j+1)%3])
				face[i][0] = face[i][1] = face[i][2] = -1;
		}
	}

	for(i=0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;
		link_N[face[i][0]]++;
		link_N[face[i][1]]++;
		link_N[face[i][2]]++;
	}
	int** link = new int*[vertex_N];
	for(i=0; i<vertex_N; i++){
		link[i] = new int[link_N[i]*3];
		link_N[i] = 0;
	}
	for(i=0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;
		for(int j=0; j<3; j++){
			int index = face[i][j];
			if(index < 0)
				continue;
			link[index] [3*link_N[index]] = face[i][(j+1)%3];
			link[index] [3*link_N[index]+1] = face[i][(j+2)%3];
			link[index] [3*link_N[index]+2] = i;
			link_N[index]++;
		}
	}

	if(vertex_link_v != NULL){
		for(i=0; i<vertex_N; i++)
			if(degree_v[i] != 0)
				delete[] vertex_link_v[i];
		delete[] vertex_link_v;
		delete[] degree_v;
	}
	vertex_link_v = new int*[vertex_N];
	degree_v = new int[vertex_N];

	if(vertex_link_f != NULL){
		for(i=0; i<vertex_N; i++)
			if(degree_f[i] != 0)
				delete[] vertex_link_f[i];
			delete[] vertex_link_f;
	}
	vertex_link_f = new int*[vertex_N];
	degree_f = new int[vertex_N];

	if(isBound != NULL)
		delete[] isBound;
	isBound = new BOOL[vertex_N];

		
	//OutputDebugString(_T("begin delete\n"));
	for(int i=0; i<vertex_N; i++){
		isBound[i] = this->sortVertexLink(link[i], link_N[i], vertex_link_v[i], vertex_link_f[i], 
										  degree_v[i], degree_f[i]);
	}
	//OutputDebugString(_T("not delete\n"));

	delete[] link;
	delete[] link_N;
	//delete link;
	//delete[] vertex_link;
}

void MeshData::generateFaceLink()
{
	if(face_link_E != NULL)
		delete[] face_link_E;

	face_link_E = (int (*)[3])new int[face_N][3];
	for(int i=0; i<face_N; i++)
		face_link_E[i][0] = face_link_E[i][1] = face_link_E[i][2] = -1;

	for(i=0; i<vertex_N; i++){
		if(isBound[i] == INNER){
			int* link_f = vertex_link_f[i];
			int* link_v = vertex_link_v[i];
			int f = degree_f[i];
			for(int j=0; j<f; j++){
				if(face[link_f[j]][0] == i)
					face_link_E[link_f[j]][2] = link_f[(j+1)%f];
				else if(face[link_f[j]][1] == i)
					face_link_E[link_f[j]][0] = link_f[(j+1)%f];
				else
					face_link_E[link_f[j]][1] = link_f[(j+1)%f];

				if(isBound[link_v[(j+1)%f]] == NON_MANIFOLD){
					if(face[link_f[(j+1)%f]][0] == i)
						face_link_E[link_f[(j+1)%f]][0] = link_f[j];
					else if(face[link_f[(j+1)%f]][1] == i)
						face_link_E[link_f[(j+1)%f]][1] = link_f[j];
					else
						face_link_E[link_f[(j+1)%f]][2] = link_f[j];
				}
			}
		}
		else if(isBound[i] == BOUNDARY){
			int* link_f = vertex_link_f[i];
			int* link_v = vertex_link_v[i];
			int f = degree_f[i];
			for(int j=0; j<f-1; j++){
				if(face[link_f[j]][0] == i)
					face_link_E[link_f[j]][2] = link_f[j+1];
				else if(face[link_f[j]][1] == i)
					face_link_E[link_f[j]][0] = link_f[j+1];
				else
					face_link_E[link_f[j]][1] = link_f[j+1];

				if(isBound[link_v[j+1]] == NON_MANIFOLD){
					if(face[link_f[j+1]][0] == i)
						face_link_E[link_f[j+1]][0] = link_f[j];
					else if(face[link_f[j+1]][1] == i)
						face_link_E[link_f[j+1]][1] = link_f[j];
					else
						face_link_E[link_f[j+1]][2] = link_f[j];
				}
			}
		}
	}
}


void MeshData::get1Ring(int index, int *&nei, int &nei_N)
{
	nei_N = degree_v[index];
	nei = vertex_link_v[index];
}
//计算面法向
void MeshData::computeFaceNormal()
{
	if(normal_f != NULL)
		delete[] normal_f;
	normal_f = new float[face_N][3];
	for(int i=0; i<face_N; i++){
		float v1[3], v2[3];
		double n[3];
		double length;
		VEC(v1, vertex[face[i][0]], vertex[face[i][1]]);//获得两个向量
		VEC(v2, vertex[face[i][0]], vertex[face[i][2]]);
		CROSS(n, v1, v2);// n = (p2-p0)*(p1-p0)
		length = LENGTH(n); 
		if((float)length != 0){
			normal_f[i][0] = (float)(n[0]/length);
			normal_f[i][1] = (float)(n[1]/length);
			normal_f[i][2] = (float)(n[2]/length);
		}
		else
			normal_f[i][0] = normal_f[i][1] = normal_f[i][2] = 0.0f;
	}
}

void MeshData::computeNormal()
{
	if(normal != NULL)
		delete[] normal;
	normal = new float[vertex_N][3];

	double (*tmp)[3] = new double[vertex_N][3];
	for(int i=0; i<vertex_N; i++)
		tmp[i][0] = tmp[i][1] = tmp[i][2] = 0;

	for(i=0; i<face_N; i++){
		int *f = face[i];
		double area = AREA(vertex[f[0]], vertex[f[1]], vertex[f[2]]);
		float v1[3], v2[3];
		double c[3];
		VEC(v1, vertex[f[1]], vertex[f[0]]);
		VEC(v2, vertex[f[2]], vertex[f[0]]);
		CROSS(c, v1, v2);

		tmp[f[0]][0] += c[0];
		tmp[f[0]][1] += c[1];
		tmp[f[0]][2] += c[2];

		tmp[f[1]][0] += c[0];
		tmp[f[1]][1] += c[1];
		tmp[f[1]][2] += c[2];

		tmp[f[2]][0] += c[0];
		tmp[f[2]][1] += c[1];
		tmp[f[2]][2] += c[2];
	}

	for(i=0; i<vertex_N; i++){
		double length = LENGTH(tmp[i]);
		if((float)length != 0){
			normal[i][0] = (float)(tmp[i][0]/length);
			normal[i][1] = (float)(tmp[i][1]/length);
			normal[i][2] = (float)(tmp[i][2]/length);
		}
		else
			normal[i][0] = normal[i][1] = normal[i][2] = 0.0f;
	}

	delete[] tmp;
}

void MeshData::laplacian(int index, double* L)
{
	L[0] = L[1] = L[2] = 0;
	int n = degree_v[index];
	if(n == 0)
		return;
	int* link_v = vertex_link_v[index];
	for(int i=0; i<n; i++){
		L[0] += vertex[link_v[i]][0];
		L[1] += vertex[link_v[i]][1];
		L[2] += vertex[link_v[i]][2];
	}
	L[0] = L[0]/(double)n - (double)vertex[index][0];
	L[1] = L[1]/(double)n - (double)vertex[index][1];
	L[2] = L[2]/(double)n - (double)vertex[index][2];
}

void MeshData::meanCurvatureNormal(int index, double *Hn)
{
	Hn[0] = Hn[1] = Hn[2] = 0;
	if(isBound[index] != INNER)
		return;

	double A = 0;
	int n = degree_v[index];
	int* link_v = vertex_link_v[index];
	for(int i=0; i<n; i++){
		int t = link_v[i];
		int s = link_v[(i+1)%n];

		//vectors of edges of triangles
		float PQ1[3], PQ2[3], QQ[3];
		VEC(PQ1, vertex[index], vertex[s]);
		VEC(PQ2, vertex[index], vertex[t]);
		VEC(QQ, vertex[s], vertex[t]);

		//normal vector of triangle
		double n[3];
		CROSS(n, PQ1, PQ2); 

		//area of triangle
		double Ai = AREA(vertex[index], vertex[s], vertex[t]);
		
		if((float)Ai != 0){ //> 0.001){
			A += 0.5*Ai;

			double dot1 = DOT(PQ1, QQ);
			double dot2 = -DOT(PQ2,QQ);

			double cot1 = dot2/Ai;
			double cot2 = dot1/Ai;

			Hn[0] += cot1*PQ1[0] + cot2*PQ2[0];
			Hn[1] += cot1*PQ1[1] + cot2*PQ2[1];
			Hn[2] += cot1*PQ1[2] + cot2*PQ2[2];
		}
	}
	if((float)A != 0){ //> 0.006){
		Hn[0] /= A;
		Hn[1] /= A;
		Hn[2] /= A;
	}
	else{
		Hn[0] = 0;
		Hn[1] = 0;
		Hn[2] = 0;
	}
}

void MeshData::computePrincipal()
{
	if(k_max != NULL){
		delete[] k_max;
	}
	k_max = new double[vertex_N];

	if(k_min != NULL){
		delete[] k_min;
	}
	k_min = new double[vertex_N];

	if(t_max != NULL){
		delete[] t_max;
	}
	t_max = new double[vertex_N][3];

	if(t_min != NULL){
		delete[] t_min;
	}
	t_min = new double[vertex_N][3];

	if(normal == NULL)
		computeNormal();

	for(int i=0; i<vertex_N; i++){
		int nei_N, *nei;
		getConsistent1Ring(i, nei, nei_N);
		if(nei_N == 0){
			k_max[i] = k_min[i] = 0;
			t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
			t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
			continue;
		}

		if(nei_N < 3){
			k_max[i] = k_min[i] = 0;
			t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
			t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
			continue;
		}

		//generate weights 
		double total_w = 0;
		double* w = new double[nei_N];
		if(isBound[i]){
			w[0] = AREA(vertex[i], vertex[nei[0]], vertex[nei[1]]);
			total_w += w[0];
			for(int j=1; j<nei_N-1; j++){
				w[j] = AREA(vertex[i], vertex[nei[j]], vertex[nei[j-1]])
					+ AREA(vertex[i], vertex[nei[j]], vertex[nei[(j+1)%nei_N]]);
				total_w += w[j];
			}
			w[nei_N-1] = AREA(vertex[i], vertex[nei[nei_N-2]], vertex[nei[nei_N-1]]);
			total_w += w[nei_N-1]; 
		}
		else{
			for(int j=0; j<nei_N; j++){
				w[j] = AREA(vertex[i], vertex[nei[j]], vertex[nei[(j+nei_N-1)%nei_N]])
					+ AREA(vertex[i], vertex[nei[j]], vertex[nei[(j+1)%nei_N]]);
				total_w += w[j];
			}
		}

		if(total_w == 0){
			k_max[i] = k_min[i] = 0;
			t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
			t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
			continue;
		}

		double M[3][3];
		M[0][0] = M[0][1] = M[0][2] = M[1][0] = M[1][1] = M[1][2] = M[2][0] = M[2][1] = M[2][2] = 0;

		for(int m=0; m<nei_N; m++){
			int j = nei[m];
			float v[3];
			VEC(v, vertex[j], vertex[i]);

			//directinal curvature
			
			double len2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
			if((float)len2 == 0)
				continue;
			double ki = 2.0*DOT(v,normal[i])/len2;
			
			

			
			//another approximation of directinal curvature
			
			//float in = DOT(normal[i], normal[j]);
			/*
			if(in > 1)
				in = 1;
			else if (in < -1)
				in = -1;
			float ki = -acos(in)/LENGTH(v);
			*/
			
			/*
			float ki;
			if(in > 1)
				ki = 0;
			else
				ki = -(2.0*sqrt(0.5*(1.0-in)))/LENGTH(v);*/
			
			//tangent 
			double t[3];
			t[0] = (1.0-normal[i][0]*normal[i][0])*v[0] 
				- normal[i][0]*normal[i][1]*v[1] - normal[i][0]*normal[i][2]*v[2];
			t[1] = -normal[i][0]*normal[i][1]*v[0] + (1.0-normal[i][1]*normal[i][1])*v[1] 
				- normal[i][1]*normal[i][2]*v[2];
			t[2] = -normal[i][0]*normal[i][2]*v[0] - normal[i][1]*normal[i][2]*v[1] 
				+ (1.0-normal[i][2]*normal[i][2])*v[2];

			
			//another approximation of directinal curvature
			//if concave then -1 times
			//if(DOT(normal[j],t) < 0)
				//ki = -ki;

			double len = LENGTH(t);
			if((float)len != 0){
				t[0] /= len;
				t[1] /= len;
				t[2] /= len;

				for(int row=0; row<3; row++)
					for(int column=0; column<3; column++)
						M[row][column] += w[m]*ki*t[row]*t[column];
			}			
		}
		delete[] w;

		//normalize weight
		if(total_w == 0){
			k_max[i] = k_min[i] = 0;
			t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
			t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
			continue;
		}
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
				M[j][k] /= total_w;

		double e1[3] = {1.0, 0.0, 0.0};
		double W[3];
		if(normal[i][0] < 0){
			W[0] = (1.0 - normal[i][0])/sqrt(2.0-2.0*normal[i][0]);
			W[1] = -normal[i][1]/sqrt(2.0-2.0*normal[i][0]);
			W[2] = -normal[i][2]/sqrt(2.0-2.0*normal[i][0]);
		}
		else{
			W[0] = (1.0 + normal[i][0])/sqrt(2.0+2.0*normal[i][0]);
			W[1] = normal[i][1]/sqrt(2.0+2.0*normal[i][0]);
			W[2] = normal[i][2]/sqrt(2.0+2.0*normal[i][0]);
		}
		double Q[3][3];
		for(int row=0; row<3; row++)
			for(int column=0; column<3; column++)
				if(row != column)
					Q[row][column] = -2.0*W[row]*W[column];
				else
					Q[row][column] = 1.0 - 2.0*W[row]*W[column];
		double T1[3] = {Q[1][0], Q[1][1], Q[1][2]};
		double T2[3] = {Q[2][0], Q[2][1], Q[2][2]};

		double M2[2][2];
		M2[0][0] = M[0][0]*Q[0][1]*Q[0][1] + 2.0f*M[0][1]*Q[0][1]*Q[1][1]
			     + M[1][1]*Q[1][1]*Q[1][1] + 2.0f*M[0][2]*Q[0][1]*Q[1][2]
			     + 2.0f*M[1][2]*Q[1][1]*Q[1][2] + M[2][2]*Q[1][2]*Q[1][2];

		M2[0][1] = Q[0][2]*(M[0][0]*Q[0][1] + M[0][1]*Q[1][1] + M[0][2]*Q[1][2])
			     + Q[1][2]*(M[0][1]*Q[0][1] + M[1][1]*Q[1][1] + M[1][2]*Q[1][2])
			     + Q[2][2]*(M[0][2]*Q[0][1] + M[1][2]*Q[1][1] + M[2][2]*Q[1][2]);
		M2[1][0] = M2[0][1];

		M2[1][1] = M[0][0]*Q[0][2]*Q[0][2] + 2.0f*M[0][1]*Q[0][2]*Q[1][2] 
			     + M[1][1]*Q[1][2]*Q[1][2] + 2.0f*M[0][2]*Q[0][2]*Q[2][2]
				 + M[2][2]*Q[2][2]*Q[2][2] + 2.0f*M[1][2]*Q[1][2]*Q[2][2];

		double th = 0;
		if(M2[0][1] != 0)
			th = 0.5*(M2[1][1]-M2[0][0])/M2[0][1];
		double t = 1.0/(fabs(th)+sqrt(th*th+1));
		if(th < 0)
			t = -t;
		double cos = 1.0/sqrt(t*t+1);
		double sin = t*cos;
		
		double K1_tmp = M2[0][0] - t*M2[0][1];
		double K2_tmp = M2[1][1] + t*M2[0][1];

		double K1 = 3.0*K1_tmp - K2_tmp;
		double K2 = 3.0*K2_tmp - K1_tmp;
		
		if(K1 > K2){
			k_max[i] = K1;
			t_max[i][0] = cos*T1[0] - sin*T2[0];
			t_max[i][1] = cos*T1[1] - sin*T2[1];
			t_max[i][2] = cos*T1[2] - sin*T2[2];

			k_min[i] = K2;
			t_min[i][0] = sin*T1[0] + cos*T2[0];
			t_min[i][1] = sin*T1[1] + cos*T2[1];
			t_min[i][2] = sin*T1[2] + cos*T2[2];
		}
		else{
			k_min[i] = K1;
			t_min[i][0] = cos*T1[0] - sin*T2[0];
			t_min[i][1] = cos*T1[1] - sin*T2[1];
			t_min[i][2] = cos*T1[2] - sin*T2[2];

			k_max[i] = K2;
			t_max[i][0] = sin*T1[0] + cos*T2[0];
			t_max[i][1] = sin*T1[1] + cos*T2[1];
			t_max[i][2] = sin*T1[2] + cos*T2[2];
		}
	}
	for(i=0; i<vertex_N; i++){
		double n[3];
		CROSS(n, normal[i], t_max[i]);
		if(DOT(n, t_min[i])< 0){
			t_min[i][0] = -t_min[i][0];
			t_min[i][1] = -t_min[i][1];
			t_min[i][2] = -t_min[i][2];
		}	
	}
}

void MeshData::computeRidgeDirection()
{
	/*
	ridge_dir = t_min;
	ravine_dir = t_max;
	return;
	*/

	int i;
	double* dk_max = new double[vertex_N];
	double* dk_min = new double[vertex_N];

	if(ridge_dir == NULL)
		ridge_dir = new double[vertex_N][3];
	if(ravine_dir == NULL)
		ravine_dir = new double[vertex_N][3];

	for(i=0; i<vertex_N; i++){
		if(isBound[i]){
			continue;
		}
	
		int *nei_v, nei_N;
		this->getConsistent1Ring(i, nei_v, nei_N);
		if(nei_N == 0)
			continue;

		double k_max_max1, k_max_max2;
		double k_min_min1, k_min_min2;
		double l_max1, l_max2, l_min1, l_min2;

		double n1[3];
		CROSS(n1, normal[i], t_max[i]);
		double n2[3];
		CROSS(n2, normal[i], t_min[i]);

		for(int m=0; m<nei_N; m++){
			int j = nei_v[m];
			int k = nei_v[(m+1)%nei_N];

			float t1[3];
			VEC(t1, vertex[j], vertex[i]);
			float t2[3];
			VEC(t2, vertex[k], vertex[i]);

			double in1 = DOT(n1, t1);
			double in2 = DOT(n1, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = AREA(o, t1, v);
				double area2 = AREA(o, t2, v);
				double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
			
				if(isOne){
					k_max_max1 = k_max_in;
					l_max1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
				else{
					k_max_max2 = k_max_in;
					l_max2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
			}

			in1 = DOT(n2, t1);
			in2 = DOT(n2, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = AREA(o, t1, v);
				double area2 = AREA(o, t2, v);
				double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

				if(isOne){
					k_min_min1 = k_min_in;
					l_min1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
				else{
					k_min_min2 = k_min_in;
					l_min2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}	
			}
		}

		double d1 = sqrt(l_max1);
		double d2 = sqrt(l_max2);
		dk_max[i] = (k_max_max2-k_max[i])/d2 + (k_max[i]-k_max_max1)/d1 - (k_max_max2-k_max_max1)/(d1+d2);
	
		d1 = sqrt(l_min1);
		d2 = sqrt(l_min2);
		dk_min[i] = (k_min_min2-k_min[i])/d2 + (k_min[i]-k_min_min1)/d1 - (k_min_min2-k_min_min1)/(d1+d2);
	}

	for(i=0; i<vertex_N; i++){
		if(isBound[i])
			continue;

		int *nei, nei_N;
		this->getConsistent1Ring(i, nei, nei_N);
		if(nei_N == 0)
			continue;

		double k_max_max1, k_max_max2, k_max_min1, k_max_min2;
		double k_min_min1, k_min_min2, k_min_max1, k_min_max2;
		double l_max1, l_max2, l_min1, l_min2;

		double n1[3];
		CROSS(n1, normal[i], t_max[i]);
		double n2[3];
		CROSS(n2, normal[i], t_min[i]);

		for(int m=0; m<nei_N; m++){
			int j = nei[m];
			int k = nei[(m+1)%nei_N];

			float t1[3];
			VEC(t1, vertex[j], vertex[i]);
			float t2[3];
			VEC(t2, vertex[k], vertex[i]);

			
			double in1 = DOT(n1, t1);
			double in2 = DOT(n1, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = AREA(o, t1, v);
				double area2 = AREA(o, t2, v);

				double dk_max1 = dk_max[j];
				double dk_min1 = dk_min[j];
				if(DOT(t_max[i], t_max[j])<0){
					dk_max1 = -dk_max1;
					dk_min1 = -dk_min1;
				}
				
				double dk_max2 = dk_max[k];
				double dk_min2 = dk_min[k];
				if(DOT(t_max[i], t_max[k])<0){
					dk_max2 = -dk_max2;
					dk_min2 = -dk_min2;
				}

				double k_max_in = (area2*dk_max1 + area1*dk_max2)/(area1+area2);
				double k_min_in = (area2*dk_min1 + area1*dk_min2)/(area1+area2);
			
				if(isOne){
					k_max_max1 = k_max_in;
					k_min_max1 = k_min_in;
					l_max1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
				else{
					k_max_max2 = k_max_in;
					k_min_max2 = k_min_in;
					l_max2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
			}

			in1 = DOT(n2, t1);
			in2 = DOT(n2, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = AREA(o, t1, v);
				double area2 = AREA(o, t2, v);

				double dk_max1 = dk_max[j];
				double dk_min1 = dk_min[j];
				if(DOT(t_min[i], t_min[j])<0){
					dk_max1 = -dk_max1;
					dk_min1 = -dk_min1;
				}
				
				double dk_max2 = dk_max[k];
				double dk_min2 = dk_min[k];
				if(DOT(t_min[i], t_min[k])<0){
					dk_max2 = -dk_max2;
					dk_min2 = -dk_min2;
				}

				double k_max_in = (area2*dk_max1 + area1*dk_max2)/(area1+area2);
				double k_min_in = (area2*dk_min1 + area1*dk_min2)/(area1+area2);

				if(isOne){
					k_min_min1 = k_min_in;
					k_max_min1 = k_max_in;
					l_min1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
				else{
					k_min_min2 = k_min_in;
					k_max_min2 = k_max_in;
					l_min2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}	
			}
		}

		if(true){
			double d1 = sqrt(l_max1);
			double d2 = sqrt(l_max2);
		
			double a  = (k_max_max2)/d2 + (-k_max_max1)/d1 
				- (k_max_max2-k_max_max1)/(d1+d2);

			d1 = sqrt(l_min1);
			d2 = sqrt(l_min2);
		
			double b = (k_max_min2)/d2 + (-k_max_min1)/d1 
				- (k_max_min2-k_max_min1)/(d1+d2);

			ridge_dir[i][0] = b*t_max[i][0] - a*t_min[i][0];
			ridge_dir[i][1] = b*t_max[i][1] - a*t_min[i][1];
			ridge_dir[i][2] = b*t_max[i][2] - a*t_min[i][2];
	
			double len = LENGTH(ridge_dir[i]);
			ridge_dir[i][0] /= len;
			ridge_dir[i][1] /= len;
			ridge_dir[i][2] /= len;
		}

		if(true){
			double d1 = sqrt(l_min1);
			double d2 = sqrt(l_min2);
			double a = d2*(k_min_min1-dk_min[i])/(d1*(d1+d2)) - d1*(k_min_min2-dk_min[i])/(d2*(d1+d2));
			
			d1 = sqrt(l_max1);
			d2 = sqrt(l_max2);
			double b = d2*(k_min_max1-dk_min[i])/(d1*(d1+d2)) - d1*(k_min_max2-dk_min[i])/(d2*(d1+d2));

			ravine_dir[i][0] = b*t_min[i][0] - a*t_max[i][0];
			ravine_dir[i][1] = b*t_min[i][1] - a*t_max[i][1];
			ravine_dir[i][2] = b*t_min[i][2] - a*t_max[i][2];
		
			double len = LENGTH(ravine_dir[i]);
			ravine_dir[i][0] /= len;
			ravine_dir[i][1] /= len;
			ravine_dir[i][2] /= len;
		}
	}
	delete[] dk_max;
	delete[] dk_min;
}

void MeshData::getConsistent1Ring(int index, int *&nei, int &nei_N)
{
	nei_N = degree_v[index];
	nei = vertex_link_v[index];
}

int MeshData::nextFace(int face_id, int vertex_id)
{
	if(face_id < 0 || face_id >= face_N)
		return -1; 
	else if(vertex_id < 0 || vertex_id >= vertex_N)
		return -1;
	else if(face[face_id][0] == vertex_id)
		return face_link_E[face_id][2];
	else if(face[face_id][1] == vertex_id)
		return face_link_E[face_id][0];
	else if(face[face_id][2] == vertex_id)
		return face_link_E[face_id][1];
	else
		return -1;
}

int MeshData::nextVertex(int face_id, int vertex_id)
{
	if(face_id < 0 || face_id >= face_N)
		return -1; 
    else if(face[face_id][0] == vertex_id)
		return face[face_id][1];
	else if(face[face_id][1] == vertex_id)
		return face[face_id][2];
	else if(face[face_id][2] == vertex_id)
		return  face[face_id][0];
	else
		return -1;
}

int MeshData::previousFace(int face_id, int vertex_id)
{
	if(face_id < 0 || face_id >= face_N)
		return -1; 
	else if(vertex_id < 0 || vertex_id >= vertex_N)
		return -1;
	else if(face[face_id][0] == vertex_id)
		return face_link_E[face_id][0];
	else if(face[face_id][1] == vertex_id)
		return face_link_E[face_id][1];
	else if(face[face_id][2] == vertex_id)
		return face_link_E[face_id][2];
	else
		return -1;
}
//求空间体积
double MeshData::getVolume()
{
	double vol = 0;
	double g[3],n[3],v1[3],v2[3];

	for(int i=0; i<face_N; i++){
		g[0] = (vertex[ face[i][0] ][0] + vertex[ face[i][1] ][0] + vertex[ face[i][2] ][0])/3.0;
		g[1] = (vertex[ face[i][0] ][1] + vertex[ face[i][1] ][1] + vertex[ face[i][2] ][1])/3.0;
		g[2] = (vertex[ face[i][0] ][2] + vertex[ face[i][1] ][2] + vertex[ face[i][2] ][2])/3.0;

		v1[0] = vertex[ face[i][1] ][0] - vertex[ face[i][0] ][0];
		v1[1] = vertex[ face[i][1] ][1] - vertex[ face[i][0] ][1];
		v1[2] = vertex[ face[i][1] ][2] - vertex[ face[i][0] ][2];

		v2[0] = vertex[ face[i][2] ][0] - vertex[ face[i][0] ][0];
		v2[1] = vertex[ face[i][2] ][1] - vertex[ face[i][0] ][1];
		v2[2] = vertex[ face[i][2] ][2] - vertex[ face[i][0] ][2];

		n[0] = v1[1]*v2[2] - v1[2]*v2[1];
		n[1] = v1[2]*v2[0] - v1[0]*v2[2];
		n[2] = v1[0]*v2[1] - v1[1]*v2[0];

		vol += (g[0]*n[0] + g[1]*n[1] + g[2]*n[2])/6.0;
	}
	return vol;
}

void MeshData::getConsistentNRings(int index, int n, int **&nei, int *&nei_N)
{
	/*
	Node* que = new Node;
	que->v = index;

	//In this list, f is distance(the number of edges) from the target vertex 
	Node* visit = new Node;
	visit->append(index, 0);
	visit->append(index, 0);
	
	Node* currentQ = que;
	//breath first numbering
	//this code can not compute the exact rings near boundary vertex.
	while(true){
		for(Node* currentL = vertex_link[currentQ->v]; currentL->next != NULL; currentL = currentL->next->next){
			//check if currentL->v is visited
			for(Node* currentV = visit; currentV->next != NULL; currentV = currentV->next->next){
				if(currentV->v == currentL->v)
					break;
			}
			if(currentV->next == NULL){
				//search currentQ->v from visit-list for computing distance
				for(currentV = visit; currentV->next != NULL; currentV = currentV->next->next)
					if(currentV->v == currentQ->v)
						break;
				if(currentV->f < n){
					int dist = currentV->f+1;
					for(currentV = visit; currentV->next != NULL; currentV = currentV->next->next)
						if(currentV->v == currentL->next->v)
							break;
					if(currentV->next == NULL || currentV->f == dist){
						visit->append(currentL->v, dist);
						//additinal element is appended for generating triangle data 
						visit->append(currentL->next->v, dist);

						//prepend node into que
						Node* head = new Node;
						head->v = currentL->v;
						head->next = que;
						que = head;
					}
				}
			}
		}
		if(que == currentQ){
			delete que;
			break;
		}
		for(Node* currentQ2 = que; currentQ2->next != currentQ; currentQ2 = currentQ2->next);
		currentQ2->next = NULL;
		delete currentQ;

		// currentQ is setted to tail of que
		for(currentQ = que; currentQ->next != NULL; currentQ = currentQ->next);
	}

	nei_N = new int[n];
	for(int i=0; i<n; i++)
		nei_N[i] = 0;

	for(Node* currentV = visit->next->next; currentV->next != NULL; currentV = currentV->next)
		nei_N[currentV->f-1]++;

	nei = (int **)new int[n];
	int* counter = new int[n];
	for(i=0; i<n; i++){
		nei[i] = new int[nei_N[i]];
		counter[i] = 0;
	}

	for(currentV = visit->next->next; currentV->next != NULL; currentV = currentV->next){
		nei[currentV->f-1][counter[currentV->f-1]++] = currentV->v;
	}
	delete visit;
	delete[] counter;

	//sort ring index consistently
	for(i=0; i<n; i++){
		int bound_flag = -1;
		for(int j=0; j<nei_N[i]; j+=2){
			if(isBound[nei[i][j]]){
				bound_flag = j;
				break;
			}
		}
		if(bound_flag >= 0){
			int tmp1 = nei[i][0];
			int tmp2 = nei[i][1];
			nei[i][0] = nei[i][bound_flag];
			nei[i][1] = nei[i][bound_flag+1];
			nei[i][bound_flag] = tmp1;
			nei[i][bound_flag+1] = tmp2;
		}
		for(j=0; j<nei_N[i]-2; j+=2){
			int pair = nei[i][j+1];
			for(int k=j+2; k<nei_N[i]; k+= 2){
				if(pair == nei[i][k]){
					int tmp1 = nei[i][k];
					int tmp2 = nei[i][k+1];
					nei[i][k] = nei[i][j+2];
					nei[i][k+1] = nei[i][j+3];
					nei[i][j+2] = tmp1;
					nei[i][j+3] = tmp2;
				}
			}
		}
		int *nei_tmp = new int[nei_N[i]/2];
		for(j=0; j<nei_N[i]; j+=2)
			nei_tmp[j/2] = nei[i][j];
		delete nei[i];
		nei[i] = nei_tmp;
		nei_N[i] /= 2;
	}
	*/
}	

BOOL MeshData::getConsistent1Ring(int index, int *&nei_v, int *&nei_f, int &nei_N)
{
	nei_N = degree_v[index];
	nei_v = vertex_link_v[index];
	nei_f = vertex_link_f[index];

	if(isBound[index] == NON_MANIFOLD)
		return false;
	else
		return true;
}

void MeshData::deleteFace(int index)
{
	if(index < 0 || index >= face_N)
		return;

	for(int i=0; i<3; i++){
		int pair = face_link_E[index][i];
		if(pair >= 0){
			if(face_link_E[pair][0] == index)
				face_link_E[pair][0] = -1;
			else if(face_link_E[pair][1] == index)
				face_link_E[pair][1] = -1;
			else if(face_link_E[pair][2] == index)
				face_link_E[pair][2] = -1;
			face_link_E[index][i] = -1;
		}
	}
	face[index][0] = -1;
	face[index][1] = -1;
	face[index][2] = -1;
}

void MeshData::generateFaceLinkV()
{
	if(face_link_V_N != NULL)
		delete[] face_link_V_N;

	if(face_link_V != NULL){
		for(int i=0; i<vertex_N; i++)
			if(face_link_V[i] != NULL)
				delete[] face_link_V[i];
		delete[] face_link_V;
	}

	face_link_V_N = new int[face_N];
	face_link_V = new int*[face_N];
	for(int i=0; i<vertex_N; i++)
		face_link_V[i] = NULL;

	for(i=0; i<face_N; i++){
		face_link_V_N[i] = 0;
		Node* list = new Node;
		for(int j=0; j<3; j++){
			int v = face[i][j];
			int *nei_v, *nei_f, nei_N;
			getConsistent1Ring(v, nei_v, nei_f, nei_N);
			if(nei_N == 0)
				continue;
			if(isBound[v])
				nei_N--;
			for(int k=0; k<nei_N; k++){
				if(nei_f[k] != i 
					&& nei_f[k] != face_link_E[i][0]
						&& nei_f[k] != face_link_E[i][1]
							&& nei_f[k] != face_link_E[i][2]){
					list->append(-1, nei_f[k]);
					face_link_V_N[i]++;
				}
				
			}
		}
		if(face_link_V_N[i] != 0){
			face_link_V[i] = new int[face_link_V_N[i]];
			Node* current = list;
			for(j=0; j<face_link_V_N[i]; j++){
				face_link_V[i][j] = current->f;
				current = current->next;
			}
		}
		delete list;
	}
}

float MeshData::minEdge()
{
	double min = 1000000;
	for(int i=0; i<vertex_N; i++){
		int n = degree_v[i];
		int* link_v = vertex_link_v[i];
		for(int j=0; j<n; j++){
			if(i > link_v[j]){
				float v[3];
				VEC(v, vertex[i], vertex[link_v[j]]);
				double len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
				if((float)len != 0 && min > len)
					min = len;
			}
		}
	}
	return (float)min;
}

void MeshData::copyMesh(MeshData *mesh)
{
	setVertexCount(mesh->vertex_N);
	for(int i=0; i<vertex_N; i++){
		vertex[i][0] = mesh->vertex[i][0];
		vertex[i][1] = mesh->vertex[i][1];
		vertex[i][2] = mesh->vertex[i][2];
	}

	setFaceCount(mesh->face_N);
	for(i=0; i<face_N; i++){
		face[i][0] = mesh->face[i][0];
		face[i][1] = mesh->face[i][1];
		face[i][2] = mesh->face[i][2];
	}
	this->normalize_scale = mesh->normalize_scale;
	this->original_max[0] = mesh->original_max[0];
	this->original_max[1] = mesh->original_max[1];
	this->original_max[2] = mesh->original_max[2];
	this->original_min[0] = mesh->original_min[0];
	this->original_min[1] = mesh->original_min[1];
	this->original_min[2] = mesh->original_min[2];
}

void MeshData::allocateTagEdges()
{
	if(ridge_edge != NULL){
		for(int i=0; i<vertex_N; i++)
			delete ridge_edge[i];
	}
	delete[] ridge_edge;
	if(ravine_edge != NULL){
		for(int i=0; i<vertex_N; i++)
			delete ravine_edge[i];
	}
	delete[] ravine_edge;

	ridge_edge = new Node*[vertex_N];
	ravine_edge = new Node*[vertex_N];

	for(int i=0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}
}
//把坐标按比例缩放
void MeshData::rescale(double rate)
{
	for(int i=0; i<vertex_N; i++){
		vertex[i][0] = (float)(rate*vertex[i][0]);
		vertex[i][1] = (float)(rate*vertex[i][1]);;
		vertex[i][2] = (float)(rate*vertex[i][2]);
	}
}

int MeshData::sortVertexLink(Node *link, int *&link_v, int *&link_f, int &v, int &f)
{
	//count incident face
	f = 0;
	Node* current = link;
	for(; current->next!=NULL; current=current->next->next)
		f++;
	//if no face is incident, then the point is isolated.
	if(f == 0){
		v = 0;
		link_v = NULL;
		link_f = NULL;
		return NON_MANIFOLD;
	}

	//un-sorted adjacent vertex data
	link_f = new int[f];
	int *link_v_tmp = new int[2*f];
	int i = 0;
	for(current=link; current->next!=NULL; current=current->next->next){
		link_v_tmp[2*i] = current->v;
		link_v_tmp[2*i+1] = current->next->v;
		link_f[i] = current->f;
		i++;
	}

	//check topology
	int b_even = 0;
	int b_odd = 0;
	int b_index;
	for(i=0; i<f; i++){
		int target = link_v_tmp[2*i];
		int consist = 0;
		for(int j=0; j<f; j++){
			if(target == link_v_tmp[2*j+1]){
				consist++;
				
			}
		}
		if(consist != 1){
			b_index = i;
			b_even++;
		}

		target = link_v_tmp[2*i+1];
		consist = 0;
		for(j=0; j<f; j++){
			if(target == link_v_tmp[2*j])
				consist++;
		}
		if(consist != 1)
			b_odd++;
	}
	int type;
	if(b_even == 0 && b_odd == 0){
		//the vertex is inner
		type = INNER;
		v = f;

		//sort tmp data by using boble sort
		for(i=0; i<f-1; i++){
			int target = link_v_tmp[2*i+1];
			for(int j=i+1; j<f; j++){
				if(target == link_v_tmp[2*j]){
					int tmp = link_v_tmp[2*i+2];
					link_v_tmp[2*i+2] = link_v_tmp[2*j];
					link_v_tmp[2*j] = tmp;

					tmp = link_v_tmp[2*i+3];
					link_v_tmp[2*i+3] = link_v_tmp[2*j+1];
					link_v_tmp[2*j+1] = tmp;

					tmp = link_f[i+1];
					link_f[i+1] = link_f[j];
					link_f[j] = tmp;

					break;
				}
			}
		}
		link_v = new int[v];
		for(i=0; i<v; i++)
			link_v[i] = link_v_tmp[2*i];
	}
	else if(b_even == 1 && b_odd == 1){
		//the vertex is boundary
		type = BOUNDARY;
		v = f+1;
		
		//swap head data
		int tmp = link_v_tmp[0];
		link_v_tmp[0] = link_v_tmp[2*b_index];
		link_v_tmp[2*b_index] = tmp;

		tmp = link_v_tmp[1];
		link_v_tmp[1] = link_v_tmp[2*b_index+1];
		link_v_tmp[2*b_index+1] = tmp;

		tmp = link_f[0];
		link_f[0] = link_f[b_index];
		link_f[b_index] = tmp;

		//sort tmp data by using boble sort
		for(i=0; i<f-1; i++){
			int target = link_v_tmp[2*i+1];
			for(int j=i+1; j<f; j++){
				if(target == link_v_tmp[2*j]){
					tmp = link_v_tmp[2*i+2];
					link_v_tmp[2*i+2] = link_v_tmp[2*j];
					link_v_tmp[2*j] = tmp;

					tmp = link_v_tmp[2*i+3];
					link_v_tmp[2*i+3] = link_v_tmp[2*j+1];
					link_v_tmp[2*j+1] = tmp;

					tmp = link_f[i+1];
					link_f[i+1] = link_f[j];
					link_f[j] = tmp;

					break;
				}
			}
		}

		link_v = new int[v];
		for(i=0; i<f; i++)
			link_v[i] = link_v_tmp[2*i];
		link_v[i] = link_v_tmp[2*f-1];
	}
	else{
		//the vertex is a non-manifold point
		type =  NON_MANIFOLD;
		v = 2*f;
		link_v = new int[v];
		for(i=0; i<v; i++)
			link_v[i] = link_v_tmp[i];
	}

	delete[] link_v_tmp;
	return type;
}

void MeshData::deleteLinkData()
{
	if(vertex != NULL)
		delete[] vertex;

	if(face != NULL)
		delete[] face;

	if(isBound != NULL)
		delete[] isBound;

	if(vertex_link_v != NULL){
		for(int i=0; i<vertex_N; i++){
			if(vertex_link_v[i] != NULL)
				delete[] vertex_link_v[i];
		}
		delete[] vertex_link_v;
		delete[] degree_v;
	}

	if(vertex_link_f != NULL){
		for(int i=0; i<vertex_N; i++){
			if(vertex_link_f[i] != NULL)
				delete[] vertex_link_f[i];
		}
		delete[] vertex_link_f;
		delete[] degree_f;
	}

	if(face_link_E != NULL)
		delete[] face_link_E;

	if(face_link_V_N != NULL)
		delete[] face_link_V_N;

	if(face_link_V != NULL){
		for(int i=0; i<face_N; i++)
			if(face_link_V[i] != NULL)
				delete[] face_link_V[i];
		delete[] face_link_V;
	}

	if(normal != NULL)
		delete[] normal;

	if(normal_f != NULL)
		delete[] normal_f;

	if(k_max != NULL)
		delete[] k_max;

	if(k_min != NULL)
		delete[] k_min;
	
	if(t_max != NULL)
		delete[] t_max;

	if(t_min != NULL)
		delete[] t_min;

	if(isRidge != NULL)
		delete[] isRidge;

	if(isRavine != NULL)
		delete[] isRavine;

	if(ridge_dir != NULL)
		delete[] ridge_dir;

	if(ravine_dir != NULL)
		delete[] ravine_dir;

	if(k_edge != NULL)
		delete[] k_edge;

	if(ridge_edge != NULL){
		for(int i=0; i<vertex_N; i++)
			delete ridge_edge[i];
	}
	delete[] ridge_edge;
	if(ravine_edge != NULL){
		for(int i=0; i<vertex_N; i++)
			delete ravine_edge[i];
	}
	delete[] ravine_edge;

	if(sample_point != NULL){
		delete[] sample_point;

		for(int i=0; i<sample_N; i++)
			if(t_max_line_N != 0)
				delete[] t_max_line[i];
		delete[] t_max_line;
		delete[] t_max_line_N;

		for(i=0; i<sample_N; i++)
			if(t_min_line_N != 0)
				delete[] t_min_line[i];
		delete[] t_min_line;
		delete[] t_min_line_N;
	}

	if(center != NULL)
		delete [] center;

	if(triangle_id != NULL)
		delete[] triangle_id;

	if(selected_triangle != NULL)
		delete selected_triangle;

	face_N = 0;
	vertex_N = 0;

	vertex = NULL;
	face = NULL;

	vertex_link_v = NULL;
	vertex_link_f = NULL;

	degree_v = NULL;
	degree_f = NULL;

	face_link_E = NULL;

	face_link_V_N = NULL;
	face_link_V = NULL;

	isBound = NULL;

	isRidge = NULL;
	isRavine = NULL;

	ridge_tri = NULL;
	ravine_tri = NULL;

	ridge_edge = NULL;
	ravine_edge = NULL;
	
	normal = NULL;
	normal_f = NULL;

	K = NULL;
	midK = 0;
	varience = 0;

	k_max = NULL;
	k_min = NULL;
	
	t_max = NULL;
	t_min = NULL;

	k_edge = NULL;

	ridge_dir = NULL;
	ravine_dir = NULL;

	sample_point = NULL;

	center = NULL;

	triangle_id = NULL;
	selected_triangle = NULL;
}

int MeshData::sortVertexLink(int *link, int N, int *&link_v, int *&link_f, int &v, int &f)
{
	//count incident face
	f = N;
	//if no face is incident, then the point is isolated.
	if(f == 0){
		v = 0;
		link_v = NULL;
		link_f = NULL;
		return NON_MANIFOLD;
	}

	//TODO 此处内存泄漏 
	//if (link_f != NULL)delete[] link_f;  //防止内存泄漏
	//un-sorted adjacent vertex data
	link_f = new int[f];
	int *link_v_tmp = new int[2*f];
	for(int i=0; i<N; i++){
		link_v_tmp[2*i] = link[3*i];
		link_v_tmp[2*i+1] = link[3*i+1];
		link_f[i] = link[3*i+2];
	}

	int b_index = -1;
	for(i=0; i<f; i++){
		int target = link_v_tmp[2*i];
		for(int j=0; j<f; j++)
			if(target == link_v_tmp[2*j+1])
				break;
		if(j == f)
			b_index = i;
	}
	if(b_index != -1){
		int tmp = link_v_tmp[0];
		link_v_tmp[0] = link_v_tmp[2*b_index];
		link_v_tmp[2*b_index] = tmp;

		tmp = link_v_tmp[1];
		link_v_tmp[1] = link_v_tmp[2*b_index+1];
		link_v_tmp[2*b_index+1] = tmp;

		tmp = link_f[0];
		link_f[0] = link_f[b_index];
		link_f[b_index] = tmp;
	}
	//sort tmp data by using boble sort
	for(i=0; i<f-1; i++){
		int target = link_v_tmp[2*i+1];
		for(int j=i+1; j<f; j++){
			if(target == link_v_tmp[2*j]){
				int tmp = link_v_tmp[2*i+2];
				link_v_tmp[2*i+2] = link_v_tmp[2*j];
				link_v_tmp[2*j] = tmp;

				tmp = link_v_tmp[2*i+3];
				link_v_tmp[2*i+3] = link_v_tmp[2*j+1];
				link_v_tmp[2*j+1] = tmp;

				tmp = link_f[i+1];
				link_f[i+1] = link_f[j];
				link_f[j] = tmp;

				break;
			}
		}
		if(j == f)
			b_index = -2;
	}
	if(b_index < 0 && link_v_tmp[0] != link_v_tmp[2*f-1])
		b_index = -2;
	if(b_index == -1){
		for(i=0; i<f; i++){
			int target = link_v_tmp[2*i];
			int consist = 0;
			for(int j=0; j<f; j++){
				if(target == link_v_tmp[2*j+1]){
					consist++;		
				}
			}
			if(consist != 1){
				b_index = -2;
				break;
			}
		}
	}
	else if(b_index >= 0){
		for(i=1; i<f; i++){
			int target = link_v_tmp[2*i];
			int consist = 0;
			for(int j=0; j<f; j++){
				if(target == link_v_tmp[2*j+1]){
					consist++;		
				}
			}
			if(consist != 1){
				b_index = -2;
				break;
			}
		}
	}

/*
	//check topology
	int b_even = 0;
	int b_odd = 0;
	int b_index;
	for(i=0; i<f; i++){
		int target = link_v_tmp[2*i];
		int consist = 0;
		for(int j=0; j<f; j++){
			if(target == link_v_tmp[2*j+1]){
				consist++;
				
			}
		}
		if(consist != 1){
			b_index = i;
			b_even++;
		}

		target = link_v_tmp[2*i+1];
		consist = 0;
		for(j=0; j<f; j++){
			if(target == link_v_tmp[2*j])
				consist++;
		}
		if(consist != 1)
			b_odd++;
	}
	*/
	int type;
	//if(b_even == 0 && b_odd == 0){
	if(b_index == -1){
		//the vertex is inner
		type = INNER;
		v = f;
/*
		//sort tmp data by using boble sort
		for(i=0; i<f-1; i++){
			int target = link_v_tmp[2*i+1];
			for(int j=i+1; j<f; j++){
				if(target == link_v_tmp[2*j]){
					int tmp = link_v_tmp[2*i+2];
					link_v_tmp[2*i+2] = link_v_tmp[2*j];
					link_v_tmp[2*j] = tmp;

					tmp = link_v_tmp[2*i+3];
					link_v_tmp[2*i+3] = link_v_tmp[2*j+1];
					link_v_tmp[2*j+1] = tmp;

					tmp = link_f[i+1];
					link_f[i+1] = link_f[j];
					link_f[j] = tmp;

					break;
				}
			}
		}*/
		link_v = new int[v];
		for(i=0; i<v; i++)
			link_v[i] = link_v_tmp[2*i];
	}
	else if(b_index == -2){
		//the vertex is a non-manifold point
		type =  NON_MANIFOLD;
		v = 2*f;
		link_v = new int[v];
		for(i=0; i<v; i++)
			link_v[i] = link_v_tmp[i];
	}
	//else if(b_even == 1 && b_odd == 1){
	else{
		//the vertex is boundary
		type = BOUNDARY;
		v = f+1;
		/*
		//swap head data
		int tmp = link_v_tmp[0];
		link_v_tmp[0] = link_v_tmp[2*b_index];
		link_v_tmp[2*b_index] = tmp;

		tmp = link_v_tmp[1];
		link_v_tmp[1] = link_v_tmp[2*b_index+1];
		link_v_tmp[2*b_index+1] = tmp;

		tmp = link_f[0];
		link_f[0] = link_f[b_index];
		link_f[b_index] = tmp;

		//sort tmp data by using boble sort
		for(i=0; i<f-1; i++){
			int target = link_v_tmp[2*i+1];
			for(int j=i+1; j<f; j++){
				if(target == link_v_tmp[2*j]){
					tmp = link_v_tmp[2*i+2];
					link_v_tmp[2*i+2] = link_v_tmp[2*j];
					link_v_tmp[2*j] = tmp;

					tmp = link_v_tmp[2*i+3];
					link_v_tmp[2*i+3] = link_v_tmp[2*j+1];
					link_v_tmp[2*j+1] = tmp;

					tmp = link_f[i+1];
					link_f[i+1] = link_f[j];
					link_f[j] = tmp;

					break;
				}
			}
		}*/

		link_v = new int[v];
		for(i=0; i<f; i++)
			link_v[i] = link_v_tmp[2*i];
		link_v[i] = link_v_tmp[2*f-1];
	}
	/*
	else{
		//the vertex is a non-manifold point
		type =  NON_MANIFOLD;
		v = 2*f;
		link_v = new int[v];
		for(i=0; i<v; i++)
			link_v[i] = link_v_tmp[i];
	}*/

	delete[] link_v_tmp;
	return type;
}


void MeshData::faceCenter(double c[], int f)
{
	int i1 = face[f][0];
	int i2 = face[f][1];
	int i3 = face[f][2];
	c[0] = (vertex[i1][0] + vertex[i2][0] + vertex[i3][0])/3.0;
	c[1] = (vertex[i1][1] + vertex[i2][1] + vertex[i3][1])/3.0;
	c[2] = (vertex[i1][2] + vertex[i2][2] + vertex[i3][2])/3.0;
}

void MeshData::faceCenter(float c[], int f)
{
	int i1 = face[f][0];
	int i2 = face[f][1];
	int i3 = face[f][2];
	c[0] = (vertex[i1][0] + vertex[i2][0] + vertex[i3][0])/3.0f;
	c[1] = (vertex[i1][1] + vertex[i2][1] + vertex[i3][1])/3.0f;
	c[2] = (vertex[i1][2] + vertex[i2][2] + vertex[i3][2])/3.0f;
}


double MeshData::faceArea(int f)
{
	return AREA(vertex[face[f][0]], vertex[face[f][1]], vertex[face[f][2]]);
}

void MeshData::faceNormal(float n[3], int f){
	float v1[3], v2[3];
	VEC(v1, vertex[face[f][0]], vertex[face[f][1]]);
	VEC(v2, vertex[face[f][0]], vertex[face[f][2]]);
	CROSS(n, v1, v2);
	double len = LENGTH(n);
	if((float)len != 0){
		n[0] /= len;
		n[1] /= len;
		n[2] /= len;
	}
}

void MeshData::computeEdgeCurvature()
{
	if(k_edge != NULL)
		delete[] k_edge;
	k_edge = new double[face_N][3];
	if(normal_f == NULL)
		computeFaceNormal();

	for(int i=0; i<face_N; i++){
		for(int j=0; j<3; j++){
			double c[3];
			faceCenter(c, i);
			double a = faceArea(i);
			int pair = face_link_E[i][j];
			if(pair < 0){
				k_edge[i][j] = 0;
				continue;
			}
			else if(i > pair){
				double dot = DOT(normal_f[i], normal_f[pair]);
				if(dot > 1.0)
					dot = 1.0;
				else if(dot < -1.0)
					dot = -1.0;
				double a1 = faceArea(pair);
				double l = DIST(vertex[face[i][j]], vertex[face[i][(j+1)%3]]);
				if((float)a != 0 && (float)a1 != 0)
					k_edge[i][j] = acos(dot)*l/(a+a1);
				else
					k_edge[i][j] = 0;

				double c1[3];
				faceCenter(c1, pair);
				double v[3];
				MeshData::VEC(v, c, c1);
				double len = MeshData::LENGTH(v);
				double dot1, dot2;
				if((float)len != 0){
					dot1 = MeshData::DOT(v, normal_f[i])/len;
					dot2 = -MeshData::DOT(v, normal_f[pair])/len;
				}
				else{
					dot1 = 1;
					dot2 = 1;
				}

				if(dot1 > 1.0)
					dot1 = 1.0;
				else if(dot1 < -1.0)
					dot1 = -1.0;
				if(dot2 > 1.0)
					dot2 = 1.0;
				else if(dot2 < -1.0)
					dot2 = -1.0;
				double angle = acos(dot1) + acos(dot2);
				if(angle < PI)
					k_edge[i][j] = -k_edge[i][j];

				if(face_link_E[pair][0] == i)
					k_edge[pair][0] = k_edge[i][j];
				else if(face_link_E[pair][1] == i)
					k_edge[pair][1] = k_edge[i][j];
				else
					k_edge[pair][2] = k_edge[i][j];
			}
		}
	}
}

void MeshData::copyMesh2(MeshData *mesh)
{
	int f_N = mesh->countValidFace();
	int v_N = mesh->countValidVertex();
	int* vertex_table = new int[mesh->vertex_N];
	
	setVertexCount(v_N);
	int counter = 0;
	for(int i=0; i<mesh->vertex_N; i++){
		if(mesh->degree_v[i] != 0){
			vertex[counter][0] = mesh->vertex[i][0];
			vertex[counter][1] = mesh->vertex[i][1];
			vertex[counter][2] = mesh->vertex[i][2];
			vertex_table[i] = counter;
			counter++;
		}
		else
			vertex_table[i] = -1;
	}
	setFaceCount(f_N);
	counter = 0;
	for(i=0; i<mesh->face_N; i++){
		if(mesh->face[i][0] >= 0){
			face[counter][0] = vertex_table[mesh->face[i][0]];
			face[counter][1] = vertex_table[mesh->face[i][1]];
			face[counter][2] = vertex_table[mesh->face[i][2]];
			counter++;
		}
	}
	delete[] vertex_table;

	this->normalize_scale = mesh->normalize_scale;
	this->original_max[0] = mesh->original_max[0];
	this->original_max[1] = mesh->original_max[1];
	this->original_max[2] = mesh->original_max[2];
	this->original_min[0] = mesh->original_min[0];
	this->original_min[1] = mesh->original_min[1];
	this->original_min[2] = mesh->original_min[2];
}

void MeshData::quickSort(int *index, float *w, int start, int end)
{
	if(start < end){
		float v = w[end];
		int i = start-1;
		int j = end;
		while(j > i){
			for(i = i+1; w[i] < v; i++);
			for(j = j-1; w[j] > v; j--);
			float t = w[i]; 
			w[i] = w[j]; 
			w[j] = t;

			int tmp = index[i];
			index[i] = index[j];
			index[j] = tmp;
		}
		float t = w[j];
		w[j] = w[i];
		w[i] = w[end];
		w[end] = t;

		int tmp = index[j];
		index[j] = index[i];
		index[i] = index[end];
		index[end] = tmp;

		quickSort(index, w, start, i-1);
		quickSort(index, w, i+1, end);
	}
	else
		return;
}

int MeshData::barycentricCoor(int id, float p[], double &u, double &v)
{
	double a = faceArea(id);
	double a0 = AREA(p, vertex[face[id][1]], vertex[face[id][2]]);
	double a1 = AREA(p, vertex[face[id][2]], vertex[face[id][0]]);
	double a2 = AREA(p, vertex[face[id][0]], vertex[face[id][1]]);

	if(a < a1 + a2){
		u = a1/a; //(a1+a2);
		v = a2/a; //(a1+a2);
		return -1;
	}
	if(a < a2 + a0){
		u = -a1/a;
		v = a2/a; //(a2+a0);
		return -2;
	}
	if(a < a0 + a1){
		u = a1/a; //(a0+a1);
		v = -a2/a;
		return -3;
	}

	u = a1/a;
	v = a2/a;
	return 0;
}

void MeshData::traceRidgeDir(struct TRI_POINT seed, float length, struct TRI_POINT *&line, int &point_N, float step, double T)
{
	int max = 200;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = ridge_dir[f[0]][0];
	dir0[1] = ridge_dir[f[0]][1];
	dir0[2] = ridge_dir[f[0]][2];

	dir1[0] = ridge_dir[f[1]][0];
	dir1[1] = ridge_dir[f[1]][1];
	dir1[2] = ridge_dir[f[1]][2];

	dir2[0] = ridge_dir[f[2]][0];
	dir2[1] = ridge_dir[f[2]][1];
	dir2[2] = ridge_dir[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;
		double v[3];
		float p[3];
		p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
		p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
		p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		float next_p[3];
		next_p[0] = p[0] + (float)v[0];
		next_p[1] = p[1] + (float)v[1];
		next_p[2] = p[2] + (float)v[2];

		double next_u, next_v;
		int result = barycentricCoor(current.id, next_p, next_u, next_v);
		tmp1[n1].id = current.id;
		tmp1[n1].u = next_u;
		tmp1[n1].v = next_v;

		if(result < 0){
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			int base_v;
			double old_dir[3];
			if(result == -1){
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1])
					base_v = 0;
				else if(face[current.id][1] == f[1])
					base_v = 1;
				else
					base_v = 2;
			}
			else if(result == -2){
				next_u = 0;
				next_v = -c/b;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2])
					base_v = 0;
				else if(face[current.id][1] == f[2])
					base_v = 1;
				else
					base_v = 2;
			}
			else{
				next_u = c/a;
				next_v = 0;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0])
					base_v = 0;
				else if(face[current.id][1] == f[0])
					base_v = 1;
				else
					base_v = 2;
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			t1 = next_u;
			t2 = next_v;
			t0 = 1.0 - t1 - t2;
			
			next_p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
			next_p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
			next_p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
			l += (float)DIST(next_p, p);

			barycentricCoor(current.id, next_p, current.u, current.v);

			f = face[current.id];

			dir0[0] = ridge_dir[f[0]][0];
			dir0[1] = ridge_dir[f[0]][1];
			dir0[2] = ridge_dir[f[0]][2];

			dir1[0] = ridge_dir[f[1]][0];
			dir1[1] = ridge_dir[f[1]][1];
			dir1[2] = ridge_dir[f[1]][2];

			dir2[0] = ridge_dir[f[2]][0];
			dir2[1] = ridge_dir[f[2]][1];
			dir2[2] = ridge_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		else{
			l += step;
			current = tmp1[n1];
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		n1++;
		if(T > (1.0-current.u-current.v)*k_max[face[current.id][0]] 
				+ current.u*k_max[face[current.id][1]] 
				+ current.v*k_max[face[current.id][2]])
				break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -ridge_dir[f[0]][0];
	dir0[1] = -ridge_dir[f[0]][1];
	dir0[2] = -ridge_dir[f[0]][2];

	dir1[0] = -ridge_dir[f[1]][0];
	dir1[1] = -ridge_dir[f[1]][1];
	dir1[2] = -ridge_dir[f[1]][2];

	dir2[0] = -ridge_dir[f[2]][0];
	dir2[1] = -ridge_dir[f[2]][1];
	dir2[2] = -ridge_dir[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;
		double v[3];
		float p[3];
		p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
		p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
		p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		float next_p[3];
		next_p[0] = p[0] + (float)v[0];
		next_p[1] = p[1] + (float)v[1];
		next_p[2] = p[2] + (float)v[2];

		double next_u, next_v;
		int result = barycentricCoor(current.id, next_p, next_u, next_v);
		tmp2[n2].id = current.id;
		tmp2[n2].u = next_u;
		tmp2[n2].v = next_v;

		if(result < 0){
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			int base_v;
			double old_dir[3];
			if(result == -1){
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1])
					base_v = 0;
				else if(face[current.id][1] == f[1])
					base_v = 1;
				else
					base_v = 2;
			}
			else if(result == -2){
				next_u = 0;
				next_v = -c/b;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2])
					base_v = 0;
				else if(face[current.id][1] == f[2])
					base_v = 1;
				else
					base_v = 2;
			}
			else{
				next_u = c/a;
				next_v = 0;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0])
					base_v = 0;
				else if(face[current.id][1] == f[0])
					base_v = 1;
				else
					base_v = 2;
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			t1 = next_u;
			t2 = next_v;
			t0 = 1.0 - t1 - t2;
			
			next_p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
			next_p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
			next_p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
			l += (float)DIST(next_p, p);

			barycentricCoor(current.id, next_p, current.u, current.v);

			f = face[current.id];

			dir0[0] = ridge_dir[f[0]][0];
			dir0[1] = ridge_dir[f[0]][1];
			dir0[2] = ridge_dir[f[0]][2];

			dir1[0] = ridge_dir[f[1]][0];
			dir1[1] = ridge_dir[f[1]][1];
			dir1[2] = ridge_dir[f[1]][2];

			dir2[0] = ridge_dir[f[2]][0];
			dir2[1] = ridge_dir[f[2]][1];
			dir2[2] = ridge_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		else{
			l += step;
			current = tmp2[n2];
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		n2++;
		if(T > (1.0-current.u-current.v)*k_max[face[current.id][0]] 
				+ current.u*k_max[face[current.id][1]] 
				+ current.v*k_max[face[current.id][2]])
				break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}

void MeshData::traceRavineDir(tri_point seed, float length, tri_point *&line, int &point_N, float step, double T)
{
	int max = 200;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = ravine_dir[f[0]][0];
	dir0[1] = ravine_dir[f[0]][1];
	dir0[2] = ravine_dir[f[0]][2];

	dir1[0] = ravine_dir[f[1]][0];
	dir1[1] = ravine_dir[f[1]][1];
	dir1[2] = ravine_dir[f[1]][2];

	dir2[0] = ravine_dir[f[2]][0];
	dir2[1] = ravine_dir[f[2]][1];
	dir2[2] = ravine_dir[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;
		double v[3];
		float p[3];
		p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
		p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
		p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		float next_p[3];
		next_p[0] = p[0] + (float)v[0];
		next_p[1] = p[1] + (float)v[1];
		next_p[2] = p[2] + (float)v[2];

		double next_u, next_v;
		int result = barycentricCoor(current.id, next_p, next_u, next_v);
		tmp1[n1].id = current.id;
		tmp1[n1].u = next_u;
		tmp1[n1].v = next_v;

		if(result < 0){
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			int base_v;
			double old_dir[3];
			if(result == -1){
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1])
					base_v = 0;
				else if(face[current.id][1] == f[1])
					base_v = 1;
				else
					base_v = 2;
			}
			else if(result == -2){
				next_u = 0;
				next_v = -c/b;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2])
					base_v = 0;
				else if(face[current.id][1] == f[2])
					base_v = 1;
				else
					base_v = 2;
			}
			else{
				next_u = c/a;
				next_v = 0;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0])
					base_v = 0;
				else if(face[current.id][1] == f[0])
					base_v = 1;
				else
					base_v = 2;
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			t1 = next_u;
			t2 = next_v;
			t0 = 1.0 - t1 - t2;
			
			next_p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
			next_p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
			next_p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
			l += (float)DIST(next_p, p);

			barycentricCoor(current.id, next_p, current.u, current.v);

			f = face[current.id];

			dir0[0] = ravine_dir[f[0]][0];
			dir0[1] = ravine_dir[f[0]][1];
			dir0[2] = ravine_dir[f[0]][2];

			dir1[0] = ravine_dir[f[1]][0];
			dir1[1] = ravine_dir[f[1]][1];
			dir1[2] = ravine_dir[f[1]][2];

			dir2[0] = ravine_dir[f[2]][0];
			dir2[1] = ravine_dir[f[2]][1];
			dir2[2] = ravine_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		else{
			l += step;
			current = tmp1[n1];
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		n1++;
		if(T < (1.0-current.u-current.v)*k_min[face[current.id][0]] 
				+ current.u*k_min[face[current.id][1]] 
				+ current.v*k_min[face[current.id][2]])
				break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -ravine_dir[f[0]][0];
	dir0[1] = -ravine_dir[f[0]][1];
	dir0[2] = -ravine_dir[f[0]][2];

	dir1[0] = -ravine_dir[f[1]][0];
	dir1[1] = -ravine_dir[f[1]][1];
	dir1[2] = -ravine_dir[f[1]][2];

	dir2[0] = -ravine_dir[f[2]][0];
	dir2[1] = -ravine_dir[f[2]][1];
	dir2[2] = -ravine_dir[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;
		double v[3];
		float p[3];
		p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
		p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
		p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		float next_p[3];
		next_p[0] = p[0] + (float)v[0];
		next_p[1] = p[1] + (float)v[1];
		next_p[2] = p[2] + (float)v[2];

		double next_u, next_v;
		int result = barycentricCoor(current.id, next_p, next_u, next_v);
		tmp2[n2].id = current.id;
		tmp2[n2].u = next_u;
		tmp2[n2].v = next_v;

		if(result < 0){
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			int base_v;
			double old_dir[3];
			if(result == -1){
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1])
					base_v = 0;
				else if(face[current.id][1] == f[1])
					base_v = 1;
				else
					base_v = 2;
			}
			else if(result == -2){
				next_u = 0;
				next_v = -c/b;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2])
					base_v = 0;
				else if(face[current.id][1] == f[2])
					base_v = 1;
				else
					base_v = 2;
			}
			else{
				next_u = c/a;
				next_v = 0;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0])
					base_v = 0;
				else if(face[current.id][1] == f[0])
					base_v = 1;
				else
					base_v = 2;
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			t1 = next_u;
			t2 = next_v;
			t0 = 1.0 - t1 - t2;
			
			next_p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
			next_p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
			next_p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
			l += (float)DIST(next_p, p);

			barycentricCoor(current.id, next_p, current.u, current.v);

			f = face[current.id];

			dir0[0] = ravine_dir[f[0]][0];
			dir0[1] = ravine_dir[f[0]][1];
			dir0[2] = ravine_dir[f[0]][2];

			dir1[0] = ravine_dir[f[1]][0];
			dir1[1] = ravine_dir[f[1]][1];
			dir1[2] = ravine_dir[f[1]][2];

			dir2[0] = ravine_dir[f[2]][0];
			dir2[1] = ravine_dir[f[2]][1];
			dir2[2] = ravine_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		else{
			l += step;
			current = tmp2[n2];
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		n2++;
		if(T < (1.0-current.u-current.v)*k_min[face[current.id][0]] 
				+ current.u*k_min[face[current.id][1]] 
				+ current.v*k_min[face[current.id][2]])
				break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}

bool MeshData::barycentricVector(int id, double t[], double &u, double &v)
{
	float v1[3], v2[3];
	VEC(v1, vertex[face[id][0]], vertex[face[id][1]]);
	VEC(v2, vertex[face[id][0]], vertex[face[id][2]]);
	double det01 = v1[0]*v2[1] - v2[0]*v1[1];
	double det12 = v1[1]*v2[2] - v2[1]*v1[2];
	double det20 = v1[2]*v2[0] - v2[2]*v1[0];

	if(fabs(det01) > fabs(det12)){
		if(fabs(det01) > fabs(det20)){
			//max is |det01|
			if((float)det01 == 0)
				return false;
			u = (v2[1]*t[0] - v2[0]*t[1])/det01;
			v = (-v1[1]*t[0] + v1[0]*t[1])/det01;
		}
		else{
			//max is |det20|
			if((float)det20 == 0)
				return false;
			u = (v2[0]*t[2] - v2[2]*t[0])/det20;
			v = (-v1[0]*t[2] + v1[2]*t[0])/det20;
		}
	}
	else if(fabs(det12) > fabs(det20)){
		//max is |det12|
		if((float)det12 == 0)
			return false;
		u = (v2[2]*t[1] - v2[1]*t[2])/det12;
		v = (-v1[2]*t[1] + v1[1]*t[2])/det12;
	}
	else{
		//max is |det20|
		if((float)det20 == 0)
			return false;
		u = (v2[0]*t[2] - v2[2]*t[0])/det20;
		v = (-v1[0]*t[2] + v1[2]*t[0])/det20;
	}
	return true;
}

void MeshData::traceRidgeDir2(tri_point seed, float length, tri_point *&line, int &point_N, float step, double T)
{
	int max = 5000;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = ridge_dir[f[0]][0];
	dir0[1] = ridge_dir[f[0]][1];
	dir0[2] = ridge_dir[f[0]][2];

	dir1[0] = ridge_dir[f[1]][0];
	dir1[1] = ridge_dir[f[1]][1];
	dir1[2] = ridge_dir[f[1]][2];

	dir2[0] = ridge_dir[f[2]][0];
	dir2[1] = ridge_dir[f[2]][1];
	dir2[2] = ridge_dir[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp1[n1].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;
			l += step;
			current = tmp1[n1];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = ridge_dir[f[0]][0];
			dir0[1] = ridge_dir[f[0]][1];
			dir0[2] = ridge_dir[f[0]][2];

			dir1[0] = ridge_dir[f[1]][0];
			dir1[1] = ridge_dir[f[1]][1];
			dir1[2] = ridge_dir[f[1]][2];

			dir2[0] = ridge_dir[f[2]][0];
			dir2[1] = ridge_dir[f[2]][1];
			dir2[2] = ridge_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n1++;
		else
			break;
		double max_in = (1.0-current.u-current.v)*k_max[face[current.id][0]] 
							+ current.u*k_max[face[current.id][1]] 
							+ current.v*k_max[face[current.id][2]];
		double min_in = (1.0-current.u-current.v)*k_min[face[current.id][0]] 
							+ current.u*k_min[face[current.id][1]] 
							+ current.v*k_min[face[current.id][2]];
		if(fabs(max_in) < fabs(min_in))
			break;
		if(T > max_in)
			break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -ridge_dir[f[0]][0];
	dir0[1] = -ridge_dir[f[0]][1];
	dir0[2] = -ridge_dir[f[0]][2];

	dir1[0] = -ridge_dir[f[1]][0];
	dir1[1] = -ridge_dir[f[1]][1];
	dir1[2] = -ridge_dir[f[1]][2];

	dir2[0] = -ridge_dir[f[2]][0];
	dir2[1] = -ridge_dir[f[2]][1];
	dir2[2] = -ridge_dir[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp2[n2].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;
			l += step;
			current = tmp2[n2];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = ridge_dir[f[0]][0];
			dir0[1] = ridge_dir[f[0]][1];
			dir0[2] = ridge_dir[f[0]][2];

			dir1[0] = ridge_dir[f[1]][0];
			dir1[1] = ridge_dir[f[1]][1];
			dir1[2] = ridge_dir[f[1]][2];

			dir2[0] = ridge_dir[f[2]][0];
			dir2[1] = ridge_dir[f[2]][1];
			dir2[2] = ridge_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n2++;
		else
			break;
		double max_in = (1.0-current.u-current.v)*k_max[face[current.id][0]] 
							+ current.u*k_max[face[current.id][1]] 
							+ current.v*k_max[face[current.id][2]];
		double min_in = (1.0-current.u-current.v)*k_min[face[current.id][0]] 
							+ current.u*k_min[face[current.id][1]] 
							+ current.v*k_min[face[current.id][2]];
		if(fabs(max_in) < fabs(min_in))
			break;
		if(T > max_in)
			break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}

void MeshData::traceRavineDir2(tri_point seed, float length, tri_point *&line, int &point_N, float step, double T)
{
	int max = 1000;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = ravine_dir[f[0]][0];
	dir0[1] = ravine_dir[f[0]][1];
	dir0[2] = ravine_dir[f[0]][2];

	dir1[0] = ravine_dir[f[1]][0];
	dir1[1] = ravine_dir[f[1]][1];
	dir1[2] = ravine_dir[f[1]][2];

	dir2[0] = ravine_dir[f[2]][0];
	dir2[1] = ravine_dir[f[2]][1];
	dir2[2] = ravine_dir[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp1[n1].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;
			l += step;
			current = tmp1[n1];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = ravine_dir[f[0]][0];
			dir0[1] = ravine_dir[f[0]][1];
			dir0[2] = ravine_dir[f[0]][2];

			dir1[0] = ravine_dir[f[1]][0];
			dir1[1] = ravine_dir[f[1]][1];
			dir1[2] = ravine_dir[f[1]][2];

			dir2[0] = ravine_dir[f[2]][0];
			dir2[1] = ravine_dir[f[2]][1];
			dir2[2] = ravine_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n1++;
		else
			break;
		double max_in = (1.0-current.u-current.v)*k_max[face[current.id][0]] 
							+ current.u*k_max[face[current.id][1]] 
							+ current.v*k_max[face[current.id][2]];
		double min_in = (1.0-current.u-current.v)*k_min[face[current.id][0]] 
							+ current.u*k_min[face[current.id][1]] 
							+ current.v*k_min[face[current.id][2]];
		if(fabs(max_in) > fabs(min_in))
			break;
		if(T < min_in)
			break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -ravine_dir[f[0]][0];
	dir0[1] = -ravine_dir[f[0]][1];
	dir0[2] = -ravine_dir[f[0]][2];

	dir1[0] = -ravine_dir[f[1]][0];
	dir1[1] = -ravine_dir[f[1]][1];
	dir1[2] = -ravine_dir[f[1]][2];

	dir2[0] = -ravine_dir[f[2]][0];
	dir2[1] = -ravine_dir[f[2]][1];
	dir2[2] = -ravine_dir[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp2[n2].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;
			l += step;
			current = tmp2[n2];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = ravine_dir[f[0]][0];
			dir0[1] = ravine_dir[f[0]][1];
			dir0[2] = ravine_dir[f[0]][2];

			dir1[0] = ravine_dir[f[1]][0];
			dir1[1] = ravine_dir[f[1]][1];
			dir1[2] = ravine_dir[f[1]][2];

			dir2[0] = ravine_dir[f[2]][0];
			dir2[1] = ravine_dir[f[2]][1];
			dir2[2] = ravine_dir[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n2++;
		else
			break;
		double max_in = (1.0-current.u-current.v)*k_max[face[current.id][0]] 
							+ current.u*k_max[face[current.id][1]] 
							+ current.v*k_max[face[current.id][2]];
		double min_in = (1.0-current.u-current.v)*k_min[face[current.id][0]] 
							+ current.u*k_min[face[current.id][1]] 
							+ current.v*k_min[face[current.id][2]];
		if(fabs(max_in) > fabs(min_in))
			break;
		if(T < min_in)
			break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}


void MeshData::worldcoor(tri_point bary, float p[])
{
	double t1 = bary.u;
	double t2 = bary.v;
	double t0 = 1.0 - t1 - t2;
	int* f = face[bary.id];

	p[0] = (float)(t0*vertex[f[0]][0] + t1*vertex[f[1]][0] + t2*vertex[f[2]][0]);
	p[1] = (float)(t0*vertex[f[0]][1] + t1*vertex[f[1]][1] + t2*vertex[f[2]][1]);
	p[2] = (float)(t0*vertex[f[0]][2] + t1*vertex[f[1]][2] + t2*vertex[f[2]][2]);
}

void MeshData::computePrincipal2()
{
	if(k_max != NULL){
		delete[] k_max;
	}
	k_max = new double[vertex_N];

	if(k_min != NULL){
		delete[] k_min;
	}
	k_min = new double[vertex_N];

	/*
	if(t_max != NULL){
		delete[] t_max;
	}
	t_max = new double[vertex_N][3];

	if(t_min != NULL){
		delete[] t_min;
	}
	t_min = new double[vertex_N][3];
	*/

	for(int i=0; i<vertex_N; i++){
		double wk_max = 0;
		double wk_min = 0;
		double wt_max = 0;
		double wt_min = 0;
		double k_max_tmp = 0;
		double k_min_tmp = 0;
		int v_N = degree_v[i];
		int f_N = degree_f[i];
		int *nei_v = vertex_link_v[i];
		int *nei_f = vertex_link_f[i];
		for(int j=0; j<v_N; j++){
			if(isBound[nei_v[j]])
				continue;
			double w = DIST(vertex[i], vertex[nei_v[j]]);

			int f1 = nei_f[(j+f_N-1)%f_N];
			int f2 = nei_f[j];
			double k;
			if(face_link_E[f1][0] == f2)
				k = k_edge[f1][0];
			else if(face_link_E[f1][1] == f2)
				k = k_edge[f1][1];
			else
				k = k_edge[f1][2];
			if(k > 0){
				wk_max += w;
				k_max_tmp += w*k;
			}
			else if(k < 0){
				wk_min += w;
				k_min_tmp += w*k;
			}
		}
		if(wk_max != 0)
			k_max_tmp /= wk_max;
		if(wk_min != 0)
			k_min_tmp /= wk_min;
		k_max[i] = k_max_tmp;
		k_min[i] = k_min_tmp;
	}
}

void MeshData::traceTmax(tri_point seed, float length, tri_point *&line, int &point_N, float step)
{
	int max = 200;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = t_max[f[0]][0];
	dir0[1] = t_max[f[0]][1];
	dir0[2] = t_max[f[0]][2];

	dir1[0] = t_max[f[1]][0];
	dir1[1] = t_max[f[1]][1];
	dir1[2] = t_max[f[1]][2];

	dir2[0] = t_max[f[2]][0];
	dir2[1] = t_max[f[2]][1];
	dir2[2] = t_max[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp1[n1].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;
			l += step;
			current = tmp1[n1];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_max[f[0]][0];
			dir0[1] = t_max[f[0]][1];
			dir0[2] = t_max[f[0]][2];

			dir1[0] = t_max[f[1]][0];
			dir1[1] = t_max[f[1]][1];
			dir1[2] = t_max[f[1]][2];

			dir2[0] = t_max[f[2]][0];
			dir2[1] = t_max[f[2]][1];
			dir2[2] = t_max[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n1++;
		else
			break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -t_max[f[0]][0];
	dir0[1] = -t_max[f[0]][1];
	dir0[2] = -t_max[f[0]][2];

	dir1[0] = -t_max[f[1]][0];
	dir1[1] = -t_max[f[1]][1];
	dir1[2] = -t_max[f[1]][2];

	dir2[0] = -t_max[f[2]][0];
	dir2[1] = -t_max[f[2]][1];
	dir2[2] = -t_max[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp2[n2].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;
			l += step;
			current = tmp2[n2];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_max[f[0]][0];
			dir0[1] = t_max[f[0]][1];
			dir0[2] = t_max[f[0]][2];

			dir1[0] = t_max[f[1]][0];
			dir1[1] = t_max[f[1]][1];
			dir1[2] = t_max[f[1]][2];

			dir2[0] = t_max[f[2]][0];
			dir2[1] = t_max[f[2]][1];
			dir2[2] = t_max[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n2++;
		else
			break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}

void MeshData::traceTmax2(tri_point seed, float length, tri_point *&line, int &point_N, float step, float T)
{
	int max = 200;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = t_max[f[0]][0];
	dir0[1] = t_max[f[0]][1];
	dir0[2] = t_max[f[0]][2];

	dir1[0] = t_max[f[1]][0];
	dir1[1] = t_max[f[1]][1];
	dir1[2] = t_max[f[1]][2];

	dir2[0] = t_max[f[2]][0];
	dir2[1] = t_max[f[2]][1];
	dir2[2] = t_max[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp1[n1].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;
			l += step;
			current = tmp1[n1];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_max[f[0]][0];
			dir0[1] = t_max[f[0]][1];
			dir0[2] = t_max[f[0]][2];

			dir1[0] = t_max[f[1]][0];
			dir1[1] = t_max[f[1]][1];
			dir1[2] = t_max[f[1]][2];

			dir2[0] = t_max[f[2]][0];
			dir2[1] = t_max[f[2]][1];
			dir2[2] = t_max[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}

			double v[3];
			paralellTrans(v, dir1, f[1], f[0]);
			if(DOT(dir0, v) < T)
				break;
			paralellTrans(v, dir2, f[2], f[1]);
			if(DOT(dir1, v) < T)
				break;
			paralellTrans(v, dir0, f[0], f[2]);
			if(DOT(dir2, v) < T)
				break;
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n1++;
		else
			break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -t_max[f[0]][0];
	dir0[1] = -t_max[f[0]][1];
	dir0[2] = -t_max[f[0]][2];

	dir1[0] = -t_max[f[1]][0];
	dir1[1] = -t_max[f[1]][1];
	dir1[2] = -t_max[f[1]][2];

	dir2[0] = -t_max[f[2]][0];
	dir2[1] = -t_max[f[2]][1];
	dir2[2] = -t_max[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp2[n2].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;
			l += step;
			current = tmp2[n2];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_max[f[0]][0];
			dir0[1] = t_max[f[0]][1];
			dir0[2] = t_max[f[0]][2];

			dir1[0] = t_max[f[1]][0];
			dir1[1] = t_max[f[1]][1];
			dir1[2] = t_max[f[1]][2];

			dir2[0] = t_max[f[2]][0];
			dir2[1] = t_max[f[2]][1];
			dir2[2] = t_max[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}

			double v[3];
			paralellTrans(v, dir1, f[1], f[0]);
			if(DOT(dir0, v) < T)
				break;
			paralellTrans(v, dir2, f[2], f[1]);
			if(DOT(dir1, v) < T)
				break;
			paralellTrans(v, dir0, f[0], f[2]);
			if(DOT(dir2, v) < T)
				break;
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n2++;
		else
			break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}

void MeshData::traceTmin(tri_point seed, float length, tri_point *&line, int &point_N, float step)
{
	int max = 200;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = t_min[f[0]][0];
	dir0[1] = t_min[f[0]][1];
	dir0[2] = t_min[f[0]][2];

	dir1[0] = t_min[f[1]][0];
	dir1[1] = t_min[f[1]][1];
	dir1[2] = t_min[f[1]][2];

	dir2[0] = t_min[f[2]][0];
	dir2[1] = t_min[f[2]][1];
	dir2[2] = t_min[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp1[n1].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;
			l += step;
			current = tmp1[n1];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_min[f[0]][0];
			dir0[1] = t_min[f[0]][1];
			dir0[2] = t_min[f[0]][2];

			dir1[0] = t_min[f[1]][0];
			dir1[1] = t_min[f[1]][1];
			dir1[2] = t_min[f[1]][2];

			dir2[0] = t_min[f[2]][0];
			dir2[1] = t_min[f[2]][1];
			dir2[2] = t_min[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n1++;
		else
			break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -t_min[f[0]][0];
	dir0[1] = -t_min[f[0]][1];
	dir0[2] = -t_min[f[0]][2];

	dir1[0] = -t_min[f[1]][0];
	dir1[1] = -t_min[f[1]][1];
	dir1[2] = -t_min[f[1]][2];

	dir2[0] = -t_min[f[2]][0];
	dir2[1] = -t_min[f[2]][1];
	dir2[2] = -t_min[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp2[n2].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;
			l += step;
			current = tmp2[n2];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_min[f[0]][0];
			dir0[1] = t_min[f[0]][1];
			dir0[2] = t_min[f[0]][2];

			dir1[0] = t_min[f[1]][0];
			dir1[1] = t_min[f[1]][1];
			dir1[2] = t_min[f[1]][2];

			dir2[0] = t_min[f[2]][0];
			dir2[1] = t_min[f[2]][1];
			dir2[2] = t_min[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n2++;
		else
			break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}
	
void MeshData::traceTmin2(tri_point seed, float length, tri_point *&line, int &point_N, float step, float T)
{
	int max = 200;
	int n1 = 0;
	int n2 = 0;
	tri_point *tmp1 = new tri_point[max];
	tri_point *tmp2 = new tri_point[max];
	tri_point current;
	int *f;
	double dir0[3], dir1[3], dir2[3];

	//trace forward
	current = seed;

	f = face[current.id];

	dir0[0] = t_min[f[0]][0];
	dir0[1] = t_min[f[0]][1];
	dir0[2] = t_min[f[0]][2];

	dir1[0] = t_min[f[1]][0];
	dir1[1] = t_min[f[1]][1];
	dir1[2] = t_min[f[1]][2];

	dir2[0] = t_min[f[2]][0];
	dir2[1] = t_min[f[2]][1];
	dir2[2] = t_min[f[2]][2];

	double l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n1 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp1[n1].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;
			l += step;
			current = tmp1[n1];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp1[n1].u = next_u;
			tmp1[n1].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_min[f[0]][0];
			dir0[1] = t_min[f[0]][1];
			dir0[2] = t_min[f[0]][2];

			dir1[0] = t_min[f[1]][0];
			dir1[1] = t_min[f[1]][1];
			dir1[2] = t_min[f[1]][2];

			dir2[0] = t_min[f[2]][0];
			dir2[1] = t_min[f[2]][1];
			dir2[2] = t_min[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}

			double v[3];
			paralellTrans(v, dir1, f[1], f[0]);
			if(DOT(dir0, v) < T)
				break;
			paralellTrans(v, dir2, f[2], f[1]);
			if(DOT(dir1, v) < T)
				break;
			paralellTrans(v, dir0, f[0], f[2]);
			if(DOT(dir2, v) < T)
				break;

		}
		if(n1 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp1[n1-2], p1);
			worldcoor(tmp1[n1-1], p2);
			worldcoor(tmp1[n1], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n1++;
		else
			break;
	}

	//trace backward
	current = seed;

	f = face[current.id];

	dir0[0] = -t_min[f[0]][0];
	dir0[1] = -t_min[f[0]][1];
	dir0[2] = -t_min[f[0]][2];

	dir1[0] = -t_min[f[1]][0];
	dir1[1] = -t_min[f[1]][1];
	dir1[2] = -t_min[f[1]][2];

	dir2[0] = -t_min[f[2]][0];
	dir2[1] = -t_min[f[2]][1];
	dir2[2] = -t_min[f[2]][2];

	l = 0;

	if(DOT(dir0, dir1) < 0){
		dir1[0] = - dir1[0];
		dir1[1] = - dir1[1];
		dir1[2] = - dir1[2];
	}
	if(DOT(dir0, dir2) < 0){
		dir2[0] = - dir2[0];
		dir2[1] = - dir2[1];
		dir2[2] = - dir2[2];
	}
	
	while(l < 0.5*length && n2 < max){
		double t1 = current.u;
		double t2 = current.v;
		double t0 = 1.0 - t1 - t2;

		double v[3];
		v[0] = t0*dir0[0] + t1*dir1[0] + t2*dir2[0];
		v[1] = t0*dir0[1] + t1*dir1[1] + t2*dir2[1];
		v[2] = t0*dir0[2] + t1*dir1[2] + t2*dir2[2];
		float* n = normal_f[current.id];
		double dot = DOT(v, n);
		v[0] = v[0] - dot*n[0];
		v[1] = v[1] - dot*n[1];
		v[2] = v[2] - dot*n[2];
		double len = LENGTH(v);
		if((float)len == 0)
			break;
		v[0] = v[0]*step/len;
		v[1] = v[1]*step/len;
		v[2] = v[2]*step/len;
		
		double du, dv;
		if(!barycentricVector(current.id, v, du, dv))
			break;
		double next_u = current.u + du;
		double next_v = current.v + dv;

		tmp2[n2].id = current.id;
		//The next position is inside the same triangle.
		if(next_u + next_v <= 1.0 && next_u >= 0 && next_v >= 0){
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;
			l += step;
			current = tmp2[n2];
		}
		//Outside the triangle.
		else{
			double a = next_v - t2;
			double b = next_u - t1;
			double c = a*t1 - b*t2;

			double s0 = -c;
			double s1 = a*(1.0-t1) + b*t2;
			double s2 = -a*t1 - b*(1.0-t2);

			int result;
			if(s1*s2 > 0){
				if(next_v > 0)
					result = 1;
				else
					result = 2;
			}
			else if(s2*s0 > 0){
				if(next_v > 0)
					result = 3;
				else
					result = 2;
			}
			else if(s0*s1 > 0){
				if(next_u > 0)
					result = 3;
				else
					result = 1;
			}
			else{
				break;
			}

			int base_v;
			double old_dir[3];
			if(result == 1){
				next_u = 0;
				if(b == 0)
					break;
				next_v = -c/b;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][2];
				old_dir[0] = dir2[0];
				old_dir[1] = dir2[1];
				old_dir[2] = dir2[2];
				if(face[current.id][0] == f[2]){
					current.u = 0;
					current.v = (1.0 - next_v);
					base_v = 0;
				}
				else if(face[current.id][1] == f[2]){
					current.u = next_v;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = 1.0 - next_v;
					current.v = next_v;
					base_v = 2;
				}
			}
			else if(result == 2){
				if(a == 0)
					break;
				next_u = c/a;
				next_v = 0;
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;

				current.id = face_link_E[current.id][0];
				old_dir[0] = dir0[0];
				old_dir[1] = dir0[1];
				old_dir[2] = dir0[2];
				if(face[current.id][0] == f[0]){
					current.u = 0;
					current.v = next_u;
					base_v = 0;
				}
				else if(face[current.id][1] == f[0]){
					current.u = 1.0 - next_u;
					current.v = 0;
					base_v = 1;
				}
				else{
					current.u = next_u;
					current.v = 1.0 - next_u;
					base_v = 2;
				}
			}
			else{
				if(a+b == 0)
					break;
				next_u = (b+c)/(a+b);
				next_v = (a-c)/(a+b);
				if(next_u < 0)
					next_u = 0;
				else if(next_u > 1)
					next_u = 1;
				if(next_v < 0)
					next_v = 0;
				else if(next_v > 1)
					next_v = 1;

				current.id = face_link_E[current.id][1];
				old_dir[0] = dir1[0];
				old_dir[1] = dir1[1];
				old_dir[2] = dir1[2];
				if(face[current.id][0] == f[1]){
					base_v = 0;
					current.u = 0;
					current.v = next_v;
				}
				else if(face[current.id][1] == f[1]){
					base_v = 1;
					current.u = next_u;
					current.v = 0;
				}
				else{
					base_v = 2;
					current.u = next_v;
					current.v = next_u;
				}
			}
			tmp2[n2].u = next_u;
			tmp2[n2].v = next_v;

			if(current.id < 0)
				break;

			float ud = (float)(next_u - t1);
			float vd = (float)(next_v - t2);
			float v1[3], v2[3], v12[3];
			VEC(v1, vertex[f[0]], vertex[f[1]]);
			v1[0] *= ud;
			v1[1] *= ud;
			v1[2] *= ud;
			VEC(v2, vertex[f[2]], vertex[f[2]]);
			v2[0] *= vd;
			v2[1] *= vd;
			v2[2] *= vd;
			VEC(v12, v1, v2);
			l += (float)LENGTH(v12);

			f = face[current.id];

			dir0[0] = t_min[f[0]][0];
			dir0[1] = t_min[f[0]][1];
			dir0[2] = t_min[f[0]][2];

			dir1[0] = t_min[f[1]][0];
			dir1[1] = t_min[f[1]][1];
			dir1[2] = t_min[f[1]][2];

			dir2[0] = t_min[f[2]][0];
			dir2[1] = t_min[f[2]][1];
			dir2[2] = t_min[f[2]][2];

			if(base_v == 0){
				if(DOT(old_dir, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir0, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir0, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
			}
			else if(base_v == 1){
				if(DOT(old_dir, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
				if(DOT(dir1, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir1, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
			}
			else{
				if(DOT(old_dir, dir2) < 0){
					dir2[0] = - dir2[0];
					dir2[1] = - dir2[1];
					dir2[2] = - dir2[2];
				}
				if(DOT(dir2, dir0) < 0){
					dir0[0] = - dir0[0];
					dir0[1] = - dir0[1];
					dir0[2] = - dir0[2];
				}
				if(DOT(dir2, dir1) < 0){
					dir1[0] = - dir1[0];
					dir1[1] = - dir1[1];
					dir1[2] = - dir1[2];
				}
			}

			double v[3];
			paralellTrans(v, dir1, f[1], f[0]);
			if(DOT(dir0, v) < T)
				break;
			paralellTrans(v, dir2, f[2], f[1]);
			if(DOT(dir1, v) < T)
				break;
			paralellTrans(v, dir0, f[0], f[2]);
			if(DOT(dir2, v) < T)
				break;
		}
		if(n2 > 1){
			float p1[3], p2[3], p3[3];
			worldcoor(tmp2[n2-2], p1);
			worldcoor(tmp2[n2-1], p2);
			worldcoor(tmp2[n2], p3);
			if((p1[0]-p2[0])*(p3[0]-p2[0]) + (p1[1]-p2[1])*(p3[1]-p2[1]) + (p1[2]-p2[2])*(p3[2]-p2[2]) > 0)
				break;
		}
		if(current.u >= 0 && current.u <= 1 && current.v >= 0 && current.v <= 1)
			n2++;
		else
			break;
	}

	//concatnation
	point_N = n1+n2+1;
	line = new tri_point[n1+n2+1];
	for(int i=0; i<n2; i++){
		line[i].id = tmp2[n2-i-1].id;
		line[i].u = tmp2[n2-i-1].u;
		line[i].v = tmp2[n2-i-1].v;
	}
	line[n2].id = seed.id;
	line[n2].u = seed.u;
	line[n2].v = seed.v;
	for(i=0; i<n1; i++){
		line[i+n2+1].id = tmp1[i].id;
		line[i+n2+1].u = tmp1[i].u;
		line[i+n2+1].v = tmp1[i].v;
	}
	delete tmp1;
	delete tmp2;
}

void MeshData::paralellTrans(double vj[], double vi[], int i, int j)
{
	double cross[3];
	CROSS(cross, normal[i], normal[j]);
	double len = LENGTH(cross);
	double R[3][3];
	if(len < 0.000001){
		R[0][0] = R[1][1] = R[2][2] = 1;
		R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
	}
	else{
		double dot = MeshData::DOT(normal[i], normal[j]);
		if(dot > 1.0)
			dot = 1;
		else if(dot < -1.0)
			dot = -1;
					
			GENERATE_MAT(R, acos(dot), cross);
	}
	MAT_VEC(vj, R, vi);
}

void MeshData::computeHatchTriangle(float T)
{
	is_hatch = new BOOL[face_N];
	for(int i=0; i<face_N; i++){
		is_hatch[i] = true;
		int* f = face[i];
		for(int j=0; j<3; j++){
			int i1 = f[j];
			int i2 = f[(j+1)%3];
			double d = DIST(vertex[i1], vertex[i2]);
			double v[3];
			paralellTrans(v, t_max[i2], i2, i1);
			double dot = MeshData::DOT(t_max[i1], v);
			if(dot > 1)
				dot = 1;
			else if(dot < -1)
				dot = -1;
			if(acos(fabs(dot)) > d*T){
				is_hatch[i] = false;
				break;
			}
		}
	}
}

void MeshData::tangent(double t[], int from, int to)
{
	float v[3];
	VEC(v, vertex[from], vertex[to]);
	float *n = normal[from];
	double dot = DOT(n, v);
	t[0] = v[0] - dot*n[0];
	t[1] = v[1] - dot*n[1];
	t[2] = v[2] - dot*n[2];
	double l = LENGTH(t);
	if(l != 0){
		t[0] /= l;
		t[1] /= l;
		t[2] /= l;
	}
	else
		t[0] = t[1] = t[2] = 0;
}

void MeshData::undoNormalizeScale()
{
	for(int i=0; i<vertex_N; i++){
		vertex[i][0] = vertex[i][0]*normalize_scale + 0.5f*(original_max[0] + original_min[0]);
		vertex[i][1] = vertex[i][1]*normalize_scale + 0.5f*(original_max[1] + original_min[1]);
		vertex[i][2] = vertex[i][2]*normalize_scale + 0.5f*(original_max[2] + original_min[2]);
	}
}

void MeshData::meanCurvatureNormalWithoutA(int index, double *Hn)
{
	Hn[0] = Hn[1] = Hn[2] = 0;
	if(isBound[index] != INNER)
		return;

	double A = 0;
	int n = degree_v[index];
	int* link_v = vertex_link_v[index];
	for(int i=0; i<n; i++){
		int t = link_v[i];
		int s = link_v[(i+1)%n];

		//vectors of edges of triangles
		float PQ1[3], PQ2[3], QQ[3];
		VEC(PQ1, vertex[index], vertex[s]);
		VEC(PQ2, vertex[index], vertex[t]);
		VEC(QQ, vertex[s], vertex[t]);

		//normal vector of triangle
		double n[3];
		CROSS(n, PQ1, PQ2); 

		//area of triangle
		double Ai = AREA(vertex[index], vertex[s], vertex[t]);
		
		if((float)Ai != 0){ //> 0.001){
			A += 0.5*Ai;

			double dot1 = DOT(PQ1, QQ);
			double dot2 = -DOT(PQ2,QQ);

			double cot1 = dot2/Ai;
			double cot2 = dot1/Ai;

			Hn[0] += cot1*PQ1[0] + cot2*PQ2[0];
			Hn[1] += cot1*PQ1[1] + cot2*PQ2[1];
			Hn[2] += cot1*PQ1[2] + cot2*PQ2[2];
		}
	}
	/*
	if((float)A != 0){ //> 0.006){
		Hn[0] /= A;
		Hn[1] /= A;
		Hn[2] /= A;
	}
	else{
		Hn[0] = 0;
		Hn[1] = 0;
		Hn[2] = 0;
	}
	*/
}

void MeshData::computeCenter()
{
	if(center != NULL)
		delete[] center;
	center = new float[face_N][3];
	for(int i=0; i<face_N; i++)
		faceCenter(center[i], i);
}

void MeshData::computePrincipalHeckbert()
{
	if(k_max != NULL){
		delete[] k_max;
	}
	k_max = new double[vertex_N];

	if(k_min != NULL){
		delete[] k_min;
	}
	k_min = new double[vertex_N];

	if(t_max != NULL){
		delete[] t_max;
	}
	t_max = new double[vertex_N][3];

	if(t_min != NULL){
		delete[] t_min;
	}
	t_min = new double[vertex_N][3];

	int *degree = degree_f;
	int **link = vertex_link_f;
	float **A = new float*[4];
	A[1] = new float[4];
	A[2] = new float[4];
	A[3] = new float[4];
	float *w = new float[4];
	float **v = new float*[4];
	v[1] = new float[4];
	v[2] = new float[4];
	v[3] = new float[4];
	for(int i=0; i<vertex_N; i++){
		A[1][1] = A[1][2] = A[1][3] = A[2][1] = A[2][2] = A[2][3] = A[3][1] = A[3][2] = A[3][3] = 0;
		int deg = degree[i];
		int *l = link[i];
		for(int j=0; j<deg; j++){
			int f = l[j];
			float *n = normal_f[f];
			double a = faceArea(f);
			A[1][1] += a*n[0]*n[0];
			A[1][2] += a*n[0]*n[1];
			A[1][3] += a*n[0]*n[2];
			A[2][1] += a*n[1]*n[0];
			A[2][2] += a*n[1]*n[1];
			A[2][3] += a*n[1]*n[2];
			A[3][1] += a*n[2]*n[0];
			A[3][2] += a*n[2]*n[1];
			A[3][3] += a*n[2]*n[2];
		}
		SVD::svdcmp(A, 3, 3, w, v);
		float tmp;
		float tmp_v[3];
		if(w[1] < w[2]){
			tmp = w[2];
			w[2] = w[1];
			w[1] = tmp;
			
			tmp_v[0] = v[2][1];
			tmp_v[1] = v[2][2];
			tmp_v[2] = v[2][3];
			v[2][1] = v[1][1];
			v[2][2] = v[1][2];
			v[2][3] = v[1][3];
			v[1][1] = tmp_v[0];
			v[1][2] = tmp_v[1];
			v[1][3] = tmp_v[2];
		}
		if(w[2] < w[3]){
			tmp = w[3];
			w[3] = w[2];
			w[2] = tmp;
	
			tmp_v[0] = v[3][1];
			tmp_v[1] = v[3][2];
			tmp_v[2] = v[3][3];
			v[3][1] = v[2][1];
			v[3][2] = v[2][2];
			v[3][3] = v[2][3];
			v[2][1] = tmp_v[0];
			v[2][2] = tmp_v[1];
			v[2][3] = tmp_v[2];
		}
		if(w[1] < w[2]){
			tmp = w[2];
			w[2] = w[1];
			w[1] = tmp;

			tmp_v[0] = v[2][1];
			tmp_v[1] = v[2][2];
			tmp_v[2] = v[2][3];
			v[2][1] = v[1][1];
			v[2][2] = v[1][2];
			v[2][3] = v[1][3];
			v[1][1] = tmp_v[0];
			v[1][2] = tmp_v[1];
			v[1][3] = tmp_v[2];
		}
		
		normal[i][0] = v[3][1];
		normal[i][1] = v[3][2];
		normal[i][2] = v[3][3];

		k_max[i] = w[1];
		k_min[i] = w[2];
	}
}

int MeshData::selectTriangle(float s[], float d[])
{
	double min = 1000000000;
	int index = -1;
	for(int i=0; i<face_N; i++){
		float v[3];
		VEC(v, vertex[face[i][0]], s);
		double dot = DOT(v, normal_f[i]);
		double c = DOT(d, normal_f[i]);
		if((float)c == 0)
			continue;
		double dist = -dot/c;
		float p[3];
		p[0] = s[0] + dist*d[0];
		p[1] = s[1] + dist*d[1];
		p[2] = s[2] + dist*d[2];
		bool flag = true;
		for(int j=0; j<3; j++){
			float *p1 = vertex[face[i][j]];
			float *p2 = vertex[face[i][(j+1)%3]];
			float *p3 = vertex[face[i][(j+2)%3]];
			float e[3];
			VEC(e, p1, p2);
			double nor[3];
			CROSS(nor, normal_f[i], e);
			if((nor[0]*(p3[0] - p1[0]) + nor[1]*(p3[1] - p1[1]) + nor[2]*(p3[2] - p1[2]))*
				(nor[0]*(p[0] - p1[0]) + nor[1]*(p[1] - p1[1]) + nor[2]*(p[2] - p1[2])) < 0){
				flag = false;
				break;
			}
		}
		if(!flag)
			continue;
		if(min > dist){
			index = i;
			min = dist;
		}
	}
	if(index < 0)
		return -1;
	if(selected_triangle != NULL)
		delete selected_triangle;
	selected_triangle = new Node;
	selected_triangle->append(-1, index);
	return index;
}

float MeshData::averageOfEdgeLength()
{
	double total = 0;
	int n = 0;
	for(int i=0; i<vertex_N; i++){
		int *l = vertex_link_v[i];
		int deg = degree_v[i];
		for(int j=0; j<deg; j++){
			if(i < l[j]){
				total += DIST(vertex[i], vertex[l[j]]);
				n++;
			}
		}
	}
	return (float)(total/n);
}

BOOL MeshData::ReadModelFile(CString filePath)
{
	CString info, str;
	int vCount = 0, fCount = 0;
	//读取文件
	CStdioFile file;
	if (!file.Open(filePath, CFile::modeRead))return FALSE;
	//读取一行
	file.ReadString(str);
	//获得第一行中的数据，放到i中
	_stscanf(str, _T("%d"), &vCount);
	file.ReadString(str);
	_stscanf(str, _T("%d"), &fCount);
	//info.Format(_T("info:%d %d\n"), vCount, fCount);
	//OutputDebugString(info);
	setVertexCount(vCount);
	setFaceCount(fCount);
	//info.Format(_T("test:%d %d\n"), this->face_N, this->vertex_N);
	//OutputDebugString(info);
	//if (!in){
	//	OutputDebugString(_T("读取文件失败\n"));
	//	//return FALSE;
	//}

	//fscanf(in,"%d",&vertexCount);
	//str.Format(_T("%d"), vertexCount);
	//OutputDebugString(str);
	//fscanf(in,"%d",&faceCount);


	float p[3];
	float x, y, z;
	long i;
	for (i = 0; i < vCount; i++)
	{
		file.ReadString(str);
		//OutputDebugString(str);	
		_stscanf(str, _T("%f %f %f"), &x, &y, &z);
		//info.Format(_T("read v:%f %f %f\n"), x, y, z);
		p[0] = x, p[1] = y, p[2] = z;
		//OutputDebugString(info);
		setVertex(i, p);
	}

	int temp, f[3];
	int fx, fy, fz;
	for (i = 0; i < fCount; i++)
	{
		file.ReadString(str);
		_stscanf(str, _T("%d %d %d %d"), &temp, &fx, &fy, &fz);
		//info.Format(_T("read f:%d %d %d\n"), fx, fy, fz);
		//OutputDebugString(info);
		f[0] = fx, f[1] = fy, f[2] = fz;
		setFace(i, f);
	}
	for (long point = 0; point < this->getVertexCount(); point++){
		this->getVertex(point, p);
		//info.Format(_T("point:%d (%f %f %f)"), point, p[0], p[1], p[2]);
		//OutputDebugString(info);

	}
	computeFaceNormal();
	computeNormal();
	file.Close();
	return TRUE;
}


BOOL MeshData::Smooth(int times, float step)
{
	//Hn is mean curvature normals
	generateVertexLink();
	generateFaceLink();
	//OutputDebugString(_T("smooth:1"));
	double(*Hn)[3] = new double[vertex_N][3];

	for (int i = 0; i<times; i++){
		for (int j = 0; j<vertex_N; j++){
			//Boundary points are fixed
			//OutputDebugString(_T("smooth:2"));
			if (/*isBound[j]*/0){
				//OutputDebugString(_T("smooth:2 i"));
				Hn[j][0] = Hn[j][1] = Hn[j][2] = 0;
			}
			else{
				//smoother->MeanCurvatureFlowImplicit(step);
			}
		}
		//OutputDebugString(_T("smooth:4"));
		for (j = 0; j<vertex_N; j++){
			float *p = this->vertex[j];
			//replace each vertex
			//OutputDebugString(_T("smooth:5"));
			p[0] -= (float)(step*Hn[j][0]);
			p[1] -= (float)(step*Hn[j][1]);
			p[2] -= (float)(step*Hn[j][2]);
		}
	}
	OutputDebugString(_T("smooth:6"));
	delete[] Hn;
	return TRUE;
}

void MeshData::Test()
{
	OutputDebugString(_T("test\n"));
}