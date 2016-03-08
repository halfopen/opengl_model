#if !defined(AFX_NODE_H__0A5F9E9B_304F_43AF_A664_C329DA50A6C0__INCLUDED_)
#define AFX_NODE_H__0A5F9E9B_304F_43AF_A664_C329DA50A6C0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Node  
{
public:
	Node();
	virtual ~Node();

public:
	Node* doesInclude(int v);
	void append(int v, int f);
	int v;
	int f;
	Node* next;
	Node* tail;

};

#endif // !defined(AFX_NODE_H__0A5F9E9B_304F_43AF_A664_C329DA50A6C0__INCLUDED_)
