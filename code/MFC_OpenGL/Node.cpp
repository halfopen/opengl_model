// Node.cpp: Node クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include "Node.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

Node::Node()
{
	next = NULL;
	tail = this;
	v = -1;
	f = -1;
}

Node::~Node()
{
	/*
	if(tail != this){
		Node* current_next;
		for(Node* current=this->next; current!=NULL; current=current_next){
			current_next = current->next;
			delete current;
		}
	}*/
	if(next != NULL)
		delete next;
}

void Node::append(int v, int f)
{
	tail->v = v;
	tail->f = f;
	tail->next = new Node;
	tail = tail->next;
	/*
	for(Node* current = this; current->next != NULL; current=current->next);
	current->v = v;
	current->f = f;
	current->next = new Node;
	*/
}



Node* Node::doesInclude(int v)
{
	for(Node* current = this; current->next!=NULL; current=current->next){
		if(current->v == v)
			return current;
	}
	return NULL;
}
