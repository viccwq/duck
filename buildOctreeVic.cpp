//%20130329 duck3 �޸İ˲��������ڵ��ڰ���n��pointʱ��ֹͣsplit
//%20130514 bunny2 ʹ�����ʣ���ߺ���createOctree(����v��ֵ����Num
#include "mex.h" 
#include <math.h>
#include <iostream>
using namespace std;

const int NUM=20;

int Size;		//�������
double *Point;	//�����Ϣ

int ID=0;
int i=0;

/*
http://blog.csdn.net/timzc/article/details/6060591
2��ʵ��Octree��ԭ��
(1). �趨���ݹ����
(2). �ҳ����������ߴ磬���Դ˳ߴ罨����һ��������
(3). ���򽫵�λԪԪ�ض����ܱ�������û���ӽڵ��������
(4). ��û�дﵽ���ݹ���ȣ��ͽ���ϸ�ְ˵ȷݣ��ٽ�����������װ�ĵ�λԪԪ��ȫ���ֵ�����
����������
(5). �������������������䵽�ĵ�λԪԪ��������Ϊ���Ҹ�����������һ���ģ������������ֹͣ
ϸ�֣���Ϊ���ݿռ�ָ����ۣ�ϸ�ֵĿռ����õ��ķ���ض����٣�����һ����Ŀ��������ô����Ŀ
����һ��������������и�����Ρ�
(6). �ظ�3��ֱ���ﵽ���ݹ���ȡ�


4��BSP Tree��Octree�Ա�
a) BSP Tree�������ָ�Ϊ1���棬��Octree�ָ�Ϊ3���档
b) BSP Treeÿ���ڵ������2���ӽ�㣬��Octree�����8���ӽ��
���BSP Tree�������ڲ��ۼ�Ψ�ĳ����У���Octree��������ά����
*/

//////////////////////////////////////////////////////////////////////////
//����˲����ڵ���
template<class T>
struct OctreeNode
{
	T data; //�ڵ�����
	T xmin,xmax; //�ڵ����꣬������������������
	T ymin,ymax;
	T zmin,zmax;
	int nodedepth;
	OctreeNode <T> *top_left_front,*top_left_back; //�ýڵ�ĸ��ӽ��
	OctreeNode <T> *top_right_front,*top_right_back;
	OctreeNode <T> *bottom_left_front,*bottom_left_back;
	OctreeNode <T> *bottom_right_front,*bottom_right_back;
	OctreeNode <T> *father;
	OctreeNode //�ڵ���
		(T nodeValue = T(),
		T xminValue = T(),T xmaxValue = T(),
		T yminValue = T(),T ymaxValue = T(),
		T zminValue = T(),T zmaxValue = T(),
		int depthInf=int(),
		OctreeNode<T>* top_left_front_Node = NULL,
		OctreeNode<T>* top_left_back_Node = NULL,
		OctreeNode<T>* top_right_front_Node = NULL,
		OctreeNode<T>* top_right_back_Node = NULL,
		OctreeNode<T>* bottom_left_front_Node = NULL,
		OctreeNode<T>* bottom_left_back_Node = NULL,
		OctreeNode<T>* bottom_right_front_Node = NULL,
		OctreeNode<T>* bottom_right_back_Node = NULL,
		OctreeNode<T>* father_Node = NULL )
		:data(nodeValue),
		xmin(xminValue),xmax(xmaxValue),
		ymin(yminValue),ymax(ymaxValue),
		zmin(zminValue),zmax(zmaxValue),
		nodedepth(depthInf),
		top_left_front(top_left_front_Node),
		top_left_back(top_left_back_Node),
		top_right_front(top_right_front_Node),
		top_right_back(top_right_back_Node),
		bottom_left_front(bottom_left_front_Node),
		bottom_left_back(bottom_left_back_Node),
		bottom_right_front(bottom_right_front_Node),
		bottom_right_back(bottom_right_back_Node),
		father(father_Node)
	{}
};
//////////////////////////////////////////////////////////////////////////
//��������
struct LeafList
{
	int data;
	OctreeNode<double> * node;
	LeafList * next;
	LeafList(int iData=int(),OctreeNode<double> * pNode=NULL,LeafList * pNext=NULL)
		:data(iData),node(pNode),next(pNext)
	{}
};

//////////////////////////////////////////////////////////////////////////
OctreeNode<double> * *pointMark=NULL;
OctreeNode<double> * forFree=NULL;
LeafList * MyList=new LeafList();
LeafList * nextList=NULL;
LeafList * MyList2=new LeafList();
LeafList * nextList2=NULL;
LeafList * MyList3=new LeafList();
LeafList * nextList3=NULL;

int LeafCounter=1;
int LeafCounter2=1;
int LeafCounter3=1;

//////////////////////////////////////////////////////////////////////////
//�����˲���
template <class T>
void createOctree(OctreeNode<T> * &root,OctreeNode<T> * pfather,int depth,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
	//cout<<"�����У����Ժ򡭡�"<<endl;
	//maxdepth=maxdepth-1; //ÿ�ݹ�һ�ξͽ����ݹ����-1
	int v=0;
	int de=0;

	/*int V_top_left_front=0,V_top_left_back=0; 
	int V_top_right_front=0,V_top_right_back=0;
	int V_bottom_left_front=0,V_bottom_left_back=0;
	int V_bottom_right_front=0,V_bottom_right_back=0;
	int De_top_left_front=0,De_top_left_back=0; 
	int De_top_right_front=0,De_top_right_back=0;
	int De_bottom_left_front=0,De_bottom_left_back=0;
	int De_bottom_right_front=0,De_bottom_right_back=0;*/


	
	
	if(true)
	{
		root=new OctreeNode<T>();
		if ((ID==1)&&(NULL==int(pfather)))
		{
			for (int i=0;i<Size;i++)
			{
				pointMark[i]=root;

			}
			root->father=NULL;			
			mexPrintf("Start to build the tree...\n");
		}
		else
		{
			root->father=pfather;
		}

		root->data = ID; //Ϊ�ڵ㸳ֵ�����Դ洢�ڵ���Ϣ��������ɼ��ԡ������Ǽ�ʵ�ְ˲������ܣ��򵥸�ֵΪ��
		//mexPrintf("Current ID is:%d\n",ID);
		ID++;
		
		root->nodedepth=depth;
		depth++;

		root->xmin=xmin; //Ϊ�ڵ����긳ֵ
		root->xmax=xmax;
		root->ymin=ymin;
		root->ymax=ymax;
		root->zmin=zmin;
		root->zmax=zmax;
		double xm=(xmax-xmin)/2;//����ڵ��ά���ϵİ�߳�
		double ym=(ymax-ymin)/2;
		double zm=(zmax-zmin)/2;

		//�鿴�������е�ĸ���
		for (int i=0;i<Size;i++)
		{
			if ((Point[i]>=xmin)&&(Point[i]<xmax)&&(Point[i+Size]>=ymin)&&(Point[i+Size]<ymax)&&(Point[i+2*Size]>=zmin)&&(Point[i+2*Size]<zmax))
			{
				v++;
				pointMark[i]=root;
				if (de<(int)Point[i+3*Size])
					de=(int)Point[i+3*Size];
				
			}

		}

		//�鿴�����������е�ĸ���
		//for (int i=0;i<Size;i++)
		//{
		//	
		//	////1
		//	if ((Point[i]>=xmin)&&(Point[i]<xmax-xm)&&(Point[i+Size]>=ymax-ym)&&(Point[i+Size]<=ymax)&&(Point[i+2*Size]>=zmax-zm)&&(Point[i+2*Size]<=zmax))
		//	{
		//		V_top_left_front++;
		//		pointMark[i]=root;
		//		De_top_left_front=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->top_left_front;
		//	} 
		//	////2
		//	else if ((Point[i]>=xmin)&&(Point[i]<xmax-xm)&&(Point[i+Size]>=ymin)&&(Point[i+Size]<ymax-ym)&&(Point[i+2*Size]>=zmax-zm)&&(Point[i+2*Size]<=zmax))
		//	{
		//		V_top_left_back++;
		//		pointMark[i]=root;
		//		De_top_left_back=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->top_left_back;
		//	}
		//	////3
		//	else if ((Point[i]>=xmax-xm)&&(Point[i]<=xmax)&&(Point[i+Size]>=ymax-ym)&&(Point[i+Size]<=ymax)&&(Point[i+2*Size]>=zmax-zm)&&(Point[i+2*Size]<=zmax))
		//	{
		//		V_top_right_front++;
		//		pointMark[i]=root;
		//		De_top_right_front=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->top_right_front;
		//	}
		//	////4
		//	else if ((Point[i]>=xmax-xm)&&(Point[i]<=xmax)&&(Point[i+Size]>=ymin)&&(Point[i+Size]<ymax-ym)&&(Point[i+2*Size]>=zmax-zm)&&(Point[i+2*Size]<=zmax))
		//	{
		//		V_top_right_back++;
		//		pointMark[i]=root;
		//		De_top_right_back=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->top_right_back;
		//	}
		//	////5
		//	else if ((Point[i]>=xmin)&&(Point[i]<xmax-xm)&&(Point[i+Size]>=ymax-ym)&&(Point[i+Size]<=ymax)&&(Point[i+2*Size]>=zmin)&&(Point[i+2*Size]<zmax-zm))
		//	{
		//		V_bottom_left_front++;
		//		pointMark[i]=root;
		//		De_bottom_left_front=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->bottom_left_front;
		//	}
		//	////6
		//	else if ((Point[i]>=xmin)&&(Point[i]<xmax-xm)&&(Point[i+Size]>=ymin)&&(Point[i+Size]<ymax-ym)&&(Point[i+2*Size]>=zmin)&&(Point[i+2*Size]<zmax-zm))
		//	{
		//		V_bottom_left_back++;
		//		pointMark[i]=root;
		//		De_bottom_left_back=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->bottom_left_back;
		//	}
		//	////7
		//	else if ((Point[i]>=xmax-xm)&&(Point[i]<=xmax)&&(Point[i+Size]>=ymax-ym)&&(Point[i+Size]<=ymax)&&(Point[i+2*Size]>=zmin)&&(Point[i+2*Size]<zmax-zm))
		//	{
		//		V_bottom_right_front++;
		//		pointMark[i]=root;
		//		De_bottom_right_front=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->bottom_right_front;
		//	}
		//	////8
		//	else if ((Point[i]>=xmax-xm)&&(Point[i]<=xmax)&&(Point[i+Size]>=ymin)&&(Point[i+Size]<ymax-ym)&&(Point[i+2*Size]>=zmin)&&(Point[i+2*Size]<zmax-zm))
		//	{
		//		V_bottom_right_back++;
		//		pointMark[i]=root;
		//		De_bottom_right_back=Point[i+3*Size];
		//		//continue;
		//		//pointMark[i]=root->bottom_right_back;
		//	}
		//
		//}

		//if ((V_top_left_front>=NUM)||(V_top_left_back>=NUM)||(V_top_right_front>=NUM)||(V_top_right_back>=NUM)||
		//	(V_bottom_left_front>=NUM)||(V_bottom_left_back>=NUM)||(V_bottom_right_front>=NUM)||(V_bottom_right_back>=NUM))
		if ((v>=NUM)||((depth-1)<de))		
		{
			createOctree(root->top_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmax-zm,zmax);	
			createOctree(root->top_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmax-zm,zmax);
			createOctree(root->top_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmax-zm,zmax);
			createOctree(root->top_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmax-zm,zmax);
			createOctree(root->bottom_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmin,zmax-zm);
			createOctree(root->bottom_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmin,zmax-zm);
			createOctree(root->bottom_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmin,zmax-zm);
			createOctree(root->bottom_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmin,zmax-zm);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//�ͷŰ˲���
template <class T>
void freeOctree(OctreeNode<T> * &root)
{
	OctreeNode<double> * father=NULL,*tempNode=NULL;
	father=root->father;


	if ((root->top_left_front!=NULL)||
		(root->top_left_back!=NULL)||
		(root->top_right_front!=NULL)||
		(root->top_right_back!=NULL)||
		(root->bottom_left_front!=NULL)||
		(root->bottom_left_back!=NULL)||
		(root->bottom_right_front!=NULL)||
		(root->bottom_right_back!=NULL))
	{
		if (root->top_left_front!=NULL)		{	//�к���
			freeOctree(root->top_left_front);	root->top_left_front=NULL;	}
		if(root->top_left_back!=NULL)		{
			freeOctree(root->top_left_back);	root->top_left_back=NULL;	}
		if(root->top_right_front!=NULL)		{
			freeOctree(root->top_right_front);	root->top_right_front=NULL;	}
		if(root->top_right_back!=NULL)		{
			freeOctree(root->top_right_back);	root->top_right_back=NULL;	}	
		if(root->bottom_left_front!=NULL)	{
			freeOctree(root->bottom_left_front);root->bottom_left_front=NULL;	}
		if(root->bottom_left_back!=NULL)	{
			freeOctree(root->bottom_left_back);	root->bottom_left_back=NULL;	}
		if(root->bottom_right_front!=NULL)	{
			freeOctree(root->bottom_right_front);root->bottom_right_front=NULL;	}
		if(root->bottom_right_back!=NULL)	{
			freeOctree(root->bottom_right_back);root->bottom_right_back=NULL;	}

		tempNode=root;	
		delete tempNode;
		//mexPrintf("deleting the %5d th (ID) node\n",--ID);

	}
	else
	{
		if(root==father->top_left_front)
		{
			tempNode=root;	
			//father->top_left_front=NULL;	
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->top_left_back)
		{
			tempNode=root;
			//father->top_left_back=NULL;		
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->top_right_front)
		{
			tempNode=root;	
			//father->top_right_front=NULL;	
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->top_right_back)
		{
			tempNode=root;	
			//father->top_right_back=NULL;	
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->bottom_left_front)
		{
			tempNode=root;	
			//father->bottom_left_front=NULL;	
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->bottom_left_back)
		{	
			tempNode=root;	
			//father->bottom_left_back=NULL;	
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->bottom_right_front)
		{
			tempNode=root;	
			//father->bottom_right_front=NULL;
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}
		else if(root==father->bottom_right_back)
		{
			tempNode=root;	
			//father->bottom_right_back=NULL;	
			delete tempNode;	//mexPrintf("deleting the %5d th (ID) node\n",--ID);
		}			
		
	}	
}
//////////////////////////////////////////////////////////////////////////
//�����ڵ�1
template <class T>
void createOctree1(OctreeNode<T> * &root,OctreeNode<T> * pfather,int depth,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
	root=new OctreeNode<T>();

	root->father=pfather;


	root->data = ID; //Ϊ�ڵ㸳ֵ�����Դ洢�ڵ���Ϣ��������ɼ��ԡ������Ǽ�ʵ�ְ˲������ܣ��򵥸�ֵΪ��
	//mexPrintf("Splitting...Current ID is:%d\n",ID);
	ID++;

	root->nodedepth=depth+1;

	root->xmin=xmin; //Ϊ�ڵ����긳ֵ
	root->xmax=xmax;
	root->ymin=ymin;
	root->ymax=ymax;
	root->zmin=zmin;
	root->zmax=zmax;
	//double xm=(xmax-xmin)/2;//����ڵ��ά���ϵİ�߳�
	//double ym=(ymax-ymin)/2;
	//double zm=(ymax-ymin)/2;
	//if(depth-1<expectDepth-1)
	//{
	//	createOctree1(root->top_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmax-zm,zmax);	
	//	createOctree1(root->top_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmax-zm,zmax);
	//	createOctree1(root->top_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmax-zm,zmax);
	//	createOctree1(root->top_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmax-zm,zmax);
	//	createOctree1(root->bottom_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmin,zmax-zm);
	//	createOctree1(root->bottom_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmin,zmax-zm);
	//	createOctree1(root->bottom_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmin,zmax-zm);
	//	createOctree1(root->bottom_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmin,zmax-zm);
	//}
	//����Ҷ����ӵ��б�
	LeafList * localList=new LeafList();
	localList->data=LeafCounter;
	localList->next=NULL;
	localList->node=root;
	//mexPrintf("LeafCounter is the %-6dth leaf\n",LeafCounter);
	LeafCounter+=1;
	nextList->next=localList;
	nextList=localList;

}
//////////////////////////////////////////////////////////////////////////
//�����ڵ�2
template <class T>
void createOctree2(OctreeNode<T> * &root,OctreeNode<T> * pfather,int depth,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
	root=new OctreeNode<T>();

	root->father=pfather;


	root->data = ID; //Ϊ�ڵ㸳ֵ�����Դ洢�ڵ���Ϣ��������ɼ��ԡ������Ǽ�ʵ�ְ˲������ܣ��򵥸�ֵΪ��
	//mexPrintf("Splitting...Current ID is:%d\n",ID);
	ID++;

	root->nodedepth=depth+1;

	root->xmin=xmin; //Ϊ�ڵ����긳ֵ
	root->xmax=xmax;
	root->ymin=ymin;
	root->ymax=ymax;
	root->zmin=zmin;
	root->zmax=zmax;
	//����Ҷ����ӵ��б�2
	LeafList * localList=new LeafList();
	localList->data=LeafCounter2;
	localList->next=NULL;
	localList->node=root;
	//mexPrintf("LeafCounter2 is the %-6dth leaf\n",LeafCounter2);
	LeafCounter2+=1;
	nextList2->next=localList;
	nextList2=localList;

}
//////////////////////////////////////////////////////////////////////////
//�����ڵ�3
template <class T>
void createOctree3(OctreeNode<T> * &root,OctreeNode<T> * pfather,int depth,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
	root=new OctreeNode<T>();

	root->father=pfather;


	root->data = ID; //Ϊ�ڵ㸳ֵ�����Դ洢�ڵ���Ϣ��������ɼ��ԡ������Ǽ�ʵ�ְ˲������ܣ��򵥸�ֵΪ��
	//mexPrintf("Splitting...Current ID is:%d\n",ID);
	ID++;

	root->nodedepth=depth+1;

	root->xmin=xmin; //Ϊ�ڵ����긳ֵ
	root->xmax=xmax;
	root->ymin=ymin;
	root->ymax=ymax;
	root->zmin=zmin;
	root->zmax=zmax;
	//����Ҷ����ӵ��б�2
	LeafList * localList=new LeafList();
	localList->data=LeafCounter3;
	localList->next=NULL;
	localList->node=root;
	//mexPrintf("LeafCounter3 is the %-6dth leaf\n",LeafCounter3);
	LeafCounter3+=1;
	nextList3->next=localList;
	nextList3=localList;
}
//////////////////////////////////////////////////////////////////////////
//ϸ����Ƚ�С�Ľڵ�
//����t ��ʾϸ�ֽڵ������б�
template <class T>
void splitNode(OctreeNode<T> * &root,int t)	
{
	int depth=root->nodedepth;


	double xmin=root->xmin;
	double xmax=root->xmax;
	double ymin=root->ymin;
	double ymax=root->ymax;
	double zmin=root->zmin;
	double zmax=root->zmax;
	double xm=(xmax-xmin)/2;//����ڵ��ά���ϵİ�߳�
	double ym=(ymax-ymin)/2;
	double zm=(ymax-ymin)/2;

	if (t==1)
	{
		if (LeafCounter>1)	//�Ƿ��һ��ϸ��
		{
			LeafCounter=LeafCounter-1;
		}
		createOctree1(root->top_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmax-zm,zmax);	
		createOctree1(root->top_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmax-zm,zmax);
		createOctree1(root->top_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmax-zm,zmax);
		createOctree1(root->top_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmax-zm,zmax);
		createOctree1(root->bottom_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmin,zmax-zm);
		createOctree1(root->bottom_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmin,zmax-zm);
		createOctree1(root->bottom_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmin,zmax-zm);
		createOctree1(root->bottom_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmin,zmax-zm);
	} 
	else if(t==2)
	{
		if (LeafCounter2>1)	//�Ƿ��һ��ϸ��
		{
			LeafCounter2=LeafCounter2-1;
		}
		createOctree2(root->top_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmax-zm,zmax);	
		createOctree2(root->top_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmax-zm,zmax);
		createOctree2(root->top_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmax-zm,zmax);
		createOctree2(root->top_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmax-zm,zmax);
		createOctree2(root->bottom_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmin,zmax-zm);
		createOctree2(root->bottom_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmin,zmax-zm);
		createOctree2(root->bottom_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmin,zmax-zm);
		createOctree2(root->bottom_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmin,zmax-zm);
	}
	else if (t==3)
	{		
		if (LeafCounter3>1)	//�Ƿ��һ��ϸ��
		{
			LeafCounter3=LeafCounter3-1;
		}
		createOctree3(root->top_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmax-zm,zmax);	
		createOctree3(root->top_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmax-zm,zmax);
		createOctree3(root->top_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmax-zm,zmax);
		createOctree3(root->top_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmax-zm,zmax);
		createOctree3(root->bottom_left_front,root,depth,xmin,xmax-xm,ymax-ym,ymax,zmin,zmax-zm);
		createOctree3(root->bottom_left_back,root,depth,xmin,xmax-xm,ymin,ymax-ym,zmin,zmax-zm);
		createOctree3(root->bottom_right_front,root,depth,xmax-xm,xmax,ymax-ym,ymax,zmin,zmax-zm);
		createOctree3(root->bottom_right_back,root,depth,xmax-xm,xmax,ymin,ymax-ym,zmin,zmax-zm);
	}
}
//////////////////////////////////////////////////////////////////////////
//Ҷ�Ӵ��������б�
//����findLeaf����֮ǰ��Ҫ��ʼ�����ָ��MyList nextList
template <class T>
void findLeaf( OctreeNode<T> * & p)
{
	
	if((p!=NULL)&&(p->top_left_back!=NULL))
	{		
		findLeaf(p->top_left_front);
		findLeaf(p->top_left_back);
		findLeaf(p->top_right_front);
		findLeaf(p->top_right_back);
		findLeaf(p->bottom_left_front);
		findLeaf(p->bottom_left_back);
		findLeaf(p->bottom_right_front);
		findLeaf(p->bottom_right_back);		
	}
	else
	{
		LeafList * localList=new LeafList();
		localList->data=LeafCounter;
		localList->next=NULL;
		localList->node=p;
		//mexPrintf("LeafCounter is the %-6dth leaf\n",LeafCounter);
		LeafCounter+=1;
		nextList->next=localList;
		nextList=localList;
	}
}
//////////////////////////////////////////////////////////////////////////
//���ͬ��ȵĽڵ��Ҷ��
template <class T>
OctreeNode<double> * findLeftNeighbor( OctreeNode<T> * p)
{
	OctreeNode<double> * father=p->father,* tempNode=NULL ;
	if (p->nodedepth>0)
	{
		if (p==father->top_right_front)
		{
			return father->top_left_front;
		} 
		else if(p==father->top_right_back)
		{
			return father->top_left_back;
		}
		else if(p==father->bottom_right_front)
		{
			return father->bottom_left_front;
		}
		else if(p==father->bottom_right_back)
		{
			return father->bottom_left_back;
		}
		else //Ҫ�ҵĶ�����p�ĺ���
		{
			tempNode=findLeftNeighbor(father);
			if ((tempNode!=NULL)&&(tempNode->top_left_back!=NULL))
			{
				if (p==father->top_left_front)
				{
					return tempNode->top_right_front;
				} 
				else if(p==father->top_left_back)
				{
					return tempNode->top_right_back;
				}
				else if(p==father->bottom_left_front)
				{
					return tempNode->bottom_right_front;
				}
				else if(p==father->bottom_left_back)
				{
					return tempNode->bottom_right_back;
				}
			} 
			else
			{
				return NULL;
			}
		}
	} 
	else
	{
		return NULL;		
	}
}
//////////////////////////////////////////////////////////////////////////
//�ұ�ͬ��ȵĽڵ��Ҷ��
template <class T>
OctreeNode<double> * findRightNeighbor( OctreeNode<T> * p)
{
	OctreeNode<double> * father=p->father,* tempNode=NULL ;
	if (p->nodedepth>0)
	{
		if (p==father->top_left_front)
		{
			return father->top_right_front;
		} 
		else if(p==father->top_left_back)
		{
			return father->top_right_back;
		}
		else if(p==father->bottom_left_front)
		{
			return father->bottom_right_front;
		}
		else if(p==father->bottom_left_back)
		{
			return father->bottom_right_back;
		}
		else //Ҫ�ҵĶ�����p�ĺ���
		{
			tempNode=findRightNeighbor(father);
			if ((tempNode!=NULL)&&(tempNode->top_left_back!=NULL))
			{
				if (p==father->top_right_front)
				{
					return tempNode->top_left_front;
				} 
				else if(p==father->top_right_back)
				{
					return tempNode->top_left_back;
				}
				else if(p==father->bottom_right_front)
				{
					return tempNode->bottom_left_front;
				}
				else if(p==father->bottom_right_back)
				{
					return tempNode->bottom_left_back;
				}
			} 
			else
			{
				return NULL;
			}
		}
	} 
	else
	{
		return NULL;		
	}
}
//////////////////////////////////////////////////////////////////////////
//�ϱ�ͬ��ȵĽڵ��Ҷ��
template <class T>
OctreeNode<double> * findTopNeighbor( OctreeNode<T> * p)
{
	OctreeNode<double> * father=p->father,* tempNode=NULL ;
	if (p->nodedepth>0)
	{
		if (p==father->bottom_left_front)
		{
			return father->top_left_front;
		} 
		else if(p==father->bottom_left_back)
		{
			return father->top_left_back;
		}
		else if(p==father->bottom_right_front)
		{
			return father->top_right_front;
		}
		else if(p==father->bottom_right_back)
		{
			return father->top_right_back;
		}
		else //Ҫ�ҵĶ�����p�ĺ���
		{
			tempNode=findTopNeighbor(father);
			if ((tempNode!=NULL)&&(tempNode->top_left_back!=NULL))
			{
				if (p==father->top_left_front)
				{
					return tempNode->bottom_left_front;
				} 
				else if(p==father->top_left_back)
				{
					return tempNode->bottom_left_back;
				}
				else if(p==father->top_right_front)
				{
					return tempNode->bottom_right_front;
				}
				else if(p==father->top_right_back)
				{
					return tempNode->bottom_right_back;
				}
			} 
			else
			{
				return NULL;
			}
		}
	} 
	else
	{
		return NULL;		
	}
}
//////////////////////////////////////////////////////////////////////////
//�±�ͬ��ȵĽڵ��Ҷ��
template <class T>
OctreeNode<double> * findBottomNeighbor( OctreeNode<T> * p)
{
	OctreeNode<double> * father=p->father,* tempNode=NULL ;
	if (p->nodedepth>0)
	{
		if (p==father->top_left_front)
		{
			return father->bottom_left_front;
		} 
		else if(p==father->top_left_back)
		{
			return father->bottom_left_back;
		}
		else if(p==father->top_right_front)
		{
			return father->bottom_right_front;
		}
		else if(p==father->top_right_back)
		{
			return father->bottom_right_back;
		}
		else //Ҫ�ҵĶ�����p�ĺ���
		{
			tempNode=findBottomNeighbor(father);
			if ((tempNode!=NULL)&&(tempNode->top_left_back!=NULL))
			{
				if (p==father->bottom_left_front)
				{
					return tempNode->top_left_front;
				} 
				else if(p==father->bottom_left_back)
				{
					return tempNode->top_left_back;
				}
				else if(p==father->bottom_right_front)
				{
					return tempNode->top_right_front;
				}
				else if(p==father->bottom_right_back)
				{
					return tempNode->top_right_back;
				}
			} 
			else
			{
				return NULL;
			}
		}
	} 
	else
	{
		return NULL;		
	}
}
//////////////////////////////////////////////////////////////////////////
//ǰ��ͬ��ȵĽڵ��Ҷ��
template <class T>
OctreeNode<double> * findFrontNeighbor( OctreeNode<T> * p)
{
	OctreeNode<double> * father=p->father,* tempNode=NULL ;
	if (p->nodedepth>0)
	{
		if (p==father->top_left_back)
		{
			return father->top_left_front;
		} 
		else if(p==father->top_right_back)
		{
			return father->top_right_front;
		}
		else if(p==father->bottom_left_back)
		{
			return father->bottom_left_front;
		}
		else if(p==father->bottom_right_back)
		{
			return father->bottom_right_front;
		}
		else //Ҫ�ҵĶ�����p�ĺ���
		{
			tempNode=findFrontNeighbor(father);
			if ((tempNode!=NULL)&&(tempNode->top_left_back!=NULL))
			{
				if (p==father->top_left_front)
				{
					return tempNode->top_left_back;
				} 
				else if(p==father->top_right_front)
				{
					return tempNode->top_right_back;
				}
				else if(p==father->bottom_left_front)
				{
					return tempNode->bottom_left_back;
				}
				else if(p==father->bottom_right_front)
				{
					return tempNode->bottom_right_back;
				}
			} 
			else
			{
				return NULL;
			}
		}
	} 
	else
	{
		return NULL;		
	}
}
//////////////////////////////////////////////////////////////////////////
//���ͬ��ȵĽڵ��Ҷ��
template <class T>
OctreeNode<double> * findBackNeighbor( OctreeNode<T> * p)
{
	OctreeNode<double> * father=p->father,* tempNode=NULL ;
	if (p->nodedepth>0)
	{
		if (p==father->top_left_front)
		{
			return father->top_left_back;
		} 
		else if(p==father->top_right_front)
		{
			return father->top_right_back;
		}
		else if(p==father->bottom_left_front)
		{
			return father->bottom_left_back;
		}
		else if(p==father->bottom_right_front)
		{
			return father->bottom_right_back;
		}
		else //Ҫ�ҵĶ�����p�ĺ���
		{
			tempNode=findBackNeighbor(father);
			if ((tempNode!=NULL)&&(tempNode->top_left_back!=NULL))
			{
				if (p==father->top_left_back)
				{
					return tempNode->top_left_front;
				} 
				else if(p==father->top_right_back)
				{
					return tempNode->top_right_front;
				}
				else if(p==father->bottom_left_back)
				{
					return tempNode->bottom_left_front;
				}
				else if(p==father->bottom_right_back)
				{
					return tempNode->bottom_right_front;
				}
			} 
			else
			{
				return NULL;
			}
		}
	} 
	else
	{
		return NULL;		
	}
}
//////////////////////////////////////////////////////////////////////////
//��ȡ�ڵ㣬����ȡ��Ƚ�С��Ҷ�ӷ���
LeafList* getLeaf(LeafList * &MyList)
{
	LeafList *previousLeaf=NULL,*currentLeaf=NULL,*nextLeaf=NULL,*tempLeaf=NULL,*outputLeaf=NULL;
	int c=0,cn=0,minDepth=0;
	if (MyList->next!=NULL)
	{
		currentLeaf=MyList->next;
		if (currentLeaf->next!=NULL)
		{
			c=currentLeaf->node->nodedepth;

			nextLeaf=currentLeaf->next;
			if (nextLeaf->next!=NULL)
			{
				//Ԥ�ȸ�ֵ
				previousLeaf=MyList;
				outputLeaf=currentLeaf;//outputLeaf=MyList->next;
				//ѡ�������С�Ķ���
				while(nextLeaf->next!=NULL)
				{
					cn=nextLeaf->node->nodedepth;
					if (cn<c)
					{
						previousLeaf=currentLeaf;
						outputLeaf=nextLeaf;
						c=cn;
					/*	if (previousLeaf==NULL)
						{
							mexPrintf("");
						}*/
					} 
					
					currentLeaf=nextLeaf;
					nextLeaf=nextLeaf->next;
					
				}				
				//���һ��������ָҶ�ӵ����ֵ��Ҫ�Ƚ�
				if (nextLeaf->node->nodedepth<c)
				{
					currentLeaf->next=NULL;
					return nextLeaf;
				} 
				else
				{
					previousLeaf->next=outputLeaf->next;
					return outputLeaf;
				}
			} 
			else
			{
				//nextLeaf->nextָ��NULL���б�ֻ������������currentLeaf=MyList->next��nextLeaf
				if (nextLeaf->node->nodedepth<c)
				{
					currentLeaf->next=NULL;
					return nextLeaf;
				} 
				else
				{
					MyList->next=nextLeaf;
					return currentLeaf;
				}
			}

			//while (currentLeaf->next->next!=NULL)
			//{
			//	if (currentLeaf->next->node->nodedepth<c)
			//	{
			//		c=currentLeaf->next->node->nodedepth;
			//		

			//	}
			//	else
			//	{

			//	}
			//} 	
			//return outputLeaf;
		}
		else
		{
			//currentLeaf->nextָ��NULL���б�ֻ��һ������currentLeaf=MyList->next
			MyList->next=NULL;
			return currentLeaf;
		}

	} 
	else
	{
		//MyList->nextָ��NULL���б�Ϊ��
		return NULL;
	}	
}
//////////////////////////////////////////////////////////////////////////
//�˲���ƽ����
template <class T>
void balanceTree( OctreeNode<T> * & p)
{
	mexPrintf("Start to balance the tree...\n");
	LeafList * tempList=NULL;
	OctreeNode<double>* neighbor[6];
	bool flag;

	//////////////////////////////////////////////////////////////////////////
	//�����б�1
	flag=true;
	for (int i=0;i<6;i++)
	{
		neighbor[i]=NULL;
	}
	tempList=getLeaf(MyList);
	//
	mexPrintf("It is spliting the 1th Leaf list...\n");
	if (tempList!=NULL)
	{
		while (tempList!=NULL)
		{
			neighbor[0]=findLeftNeighbor(tempList->node);
			neighbor[1]=findRightNeighbor(tempList->node);
			neighbor[2]=findTopNeighbor(tempList->node);
			neighbor[3]=findBottomNeighbor(tempList->node);
			neighbor[4]=findFrontNeighbor(tempList->node);
			neighbor[5]=findBackNeighbor(tempList->node);

			if ((neighbor[0]!=NULL)||(neighbor[1]!=NULL)||(neighbor[2]!=NULL)||(neighbor[3]!=NULL)||(neighbor[4]!=NULL)||(neighbor[5]!=NULL))
			{
				//�жϺ�ͬ����ھ��ڽӵĺ�����û�к���
				if ((flag)&&(neighbor[0]!=NULL)&&(neighbor[0]->top_left_back!=NULL)&&
					((neighbor[0]->top_right_front->top_left_back!=NULL)||
					(neighbor[0]->top_right_back->top_left_back!=NULL)||
					(neighbor[0]->bottom_right_front->top_left_back!=NULL)||
					(neighbor[0]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,1);
					delete tempList;
					flag=false;
				} 
				else if((flag)&&(neighbor[1]!=NULL)&&(neighbor[1]->top_left_back!=NULL)&&
					((neighbor[1]->top_left_front->top_left_back!=NULL)||
					(neighbor[1]->top_left_back->top_left_back!=NULL)||
					(neighbor[1]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[1]->bottom_left_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,1);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[2]!=NULL)&&(neighbor[2]->top_left_back!=NULL)&&
					((neighbor[2]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[2]->bottom_left_back->top_left_back!=NULL)||
					(neighbor[2]->bottom_right_front->top_left_back!=NULL)||
					(neighbor[2]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,1);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[3]!=NULL)&&(neighbor[3]->top_left_back!=NULL)&&
					((neighbor[3]->top_left_front->top_left_back!=NULL)||
					(neighbor[3]->top_left_back->top_left_back!=NULL)||
					(neighbor[3]->top_right_front->top_left_back!=NULL)||
					(neighbor[3]->top_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,1);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[4]!=NULL)&&(neighbor[4]->top_left_back!=NULL)&&
					((neighbor[4]->top_left_back->top_left_back!=NULL)||
					(neighbor[4]->top_right_back->top_left_back!=NULL)||
					(neighbor[4]->bottom_left_back->top_left_back!=NULL)||
					(neighbor[4]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,1);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[5]!=NULL)&&(neighbor[5]->top_left_back!=NULL)&&
					((neighbor[5]->top_left_front->top_left_back!=NULL)||
					(neighbor[5]->top_right_front->top_left_back!=NULL)||
					(neighbor[5]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[5]->bottom_right_front->top_left_back!=NULL)))
				{
					splitNode(tempList->node,1);
					delete tempList;
					flag=false;
				}

				if(flag)
				{
					//δ�����Ҷ�����ֵ�б�2
					tempList->next=NULL;
					nextList2->next=tempList;
					nextList2=tempList;
				}
			} 
			else
			{
				//δ�����Ҷ�����ֵ�б�2
				tempList->next=NULL;
				nextList2->next=tempList;
				nextList2=tempList;
			}
			//
			flag=true;
			for (int i=0;i<6;i++)
			{
				neighbor[i]=NULL;
			}
			tempList=getLeaf(MyList);
			//
		}
		
	}
	//////////////////////////////////////////////////////////////////////////
	//�����б�2
	flag=true;
	for (int i=0;i<6;i++)
	{
		neighbor[i]=NULL;
	}
	tempList=getLeaf(MyList2);
	//
	mexPrintf("It is spliting the 2th Leaf list...\n");
	if (tempList!=NULL)
	{
		while (tempList!=NULL)
		{
			neighbor[0]=findLeftNeighbor(tempList->node);
			neighbor[1]=findRightNeighbor(tempList->node);
			neighbor[2]=findTopNeighbor(tempList->node);
			neighbor[3]=findBottomNeighbor(tempList->node);
			neighbor[4]=findFrontNeighbor(tempList->node);
			neighbor[5]=findBackNeighbor(tempList->node);

			if ((neighbor[0]!=NULL)||(neighbor[1]!=NULL)||(neighbor[2]!=NULL)||(neighbor[3]!=NULL)||(neighbor[4]!=NULL)||(neighbor[5]!=NULL))
			{
				//�жϺ�ͬ����ھ��ڽӵĺ�����û�к���
				if ((flag)&&(neighbor[0]!=NULL)&&(neighbor[0]->top_left_back!=NULL)&&
					((neighbor[0]->top_right_front->top_left_back!=NULL)||
					(neighbor[0]->top_right_back->top_left_back!=NULL)||
					(neighbor[0]->bottom_right_front->top_left_back!=NULL)||
					(neighbor[0]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,2);
					delete tempList;
					flag=false;
				} 
				else if((flag)&&(neighbor[1]!=NULL)&&(neighbor[1]->top_left_back!=NULL)&&
					((neighbor[1]->top_left_front->top_left_back!=NULL)||
					(neighbor[1]->top_left_back->top_left_back!=NULL)||
					(neighbor[1]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[1]->bottom_left_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,2);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[2]!=NULL)&&(neighbor[2]->top_left_back!=NULL)&&
					((neighbor[2]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[2]->bottom_left_back->top_left_back!=NULL)||
					(neighbor[2]->bottom_right_front->top_left_back!=NULL)||
					(neighbor[2]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,2);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[3]!=NULL)&&(neighbor[3]->top_left_back!=NULL)&&
					((neighbor[3]->top_left_front->top_left_back!=NULL)||
					(neighbor[3]->top_left_back->top_left_back!=NULL)||
					(neighbor[3]->top_right_front->top_left_back!=NULL)||
					(neighbor[3]->top_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,2);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[4]!=NULL)&&(neighbor[4]->top_left_back!=NULL)&&
					((neighbor[4]->top_left_back->top_left_back!=NULL)||
					(neighbor[4]->top_right_back->top_left_back!=NULL)||
					(neighbor[4]->bottom_left_back->top_left_back!=NULL)||
					(neighbor[4]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,2);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[5]!=NULL)&&(neighbor[5]->top_left_back!=NULL)&&
					((neighbor[5]->top_left_front->top_left_back!=NULL)||
					(neighbor[5]->top_right_front->top_left_back!=NULL)||
					(neighbor[5]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[5]->bottom_right_front->top_left_back!=NULL)))
				{
					splitNode(tempList->node,2);
					delete tempList;
					flag=false;
				}
				
				if(flag)
				{
					//δ�����Ҷ�����ֵ�б�3
					tempList->next=NULL;
					nextList3->next=tempList;
					nextList3=tempList;
				}
			} 
			else
			{
				//δ�����Ҷ�����ֵ�б�3
				tempList->next=NULL;
				nextList3->next=tempList;
				nextList3=tempList;
			}
			//
			flag=true;
			for (int i=0;i<6;i++)
			{
				neighbor[i]=NULL;
			}
			tempList=getLeaf(MyList2);	
			//
		}
	}
	/////////////////////////////////////////////////////////////////////////
	//�����б�3
	flag=true;
	for (int i=0;i<6;i++)
	{
		neighbor[i]=NULL;
	}
	tempList=getLeaf(MyList3);
	//
	mexPrintf("It is spliting the 3th Leaf list...\n");
	if (tempList!=NULL)
	{
		while (tempList!=NULL)
		{
			neighbor[0]=findLeftNeighbor(tempList->node);
			neighbor[1]=findRightNeighbor(tempList->node);
			neighbor[2]=findTopNeighbor(tempList->node);
			neighbor[3]=findBottomNeighbor(tempList->node);
			neighbor[4]=findFrontNeighbor(tempList->node);
			neighbor[5]=findBackNeighbor(tempList->node);

			if ((neighbor[0]!=NULL)||(neighbor[1]!=NULL)||(neighbor[2]!=NULL)||(neighbor[3]!=NULL)||(neighbor[4]!=NULL)||(neighbor[5]!=NULL))
			{
				//�жϺ�ͬ����ھ��ڽӵĺ�����û�к���
				if ((flag)&&(neighbor[0]!=NULL)&&(neighbor[0]->top_left_back!=NULL)&&
					((neighbor[0]->top_right_front->top_left_back!=NULL)||
					(neighbor[0]->top_right_back->top_left_back!=NULL)||
					(neighbor[0]->bottom_right_front->top_left_back!=NULL)||
					(neighbor[0]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,3);
					delete tempList;
					flag=false;
				} 
				else if((flag)&&(neighbor[1]!=NULL)&&(neighbor[1]->top_left_back!=NULL)&&
					((neighbor[1]->top_left_front->top_left_back!=NULL)||
					(neighbor[1]->top_left_back->top_left_back!=NULL)||
					(neighbor[1]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[1]->bottom_left_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,3);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[2]!=NULL)&&(neighbor[2]->top_left_back!=NULL)&&
					((neighbor[2]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[2]->bottom_left_back->top_left_back!=NULL)||
					(neighbor[2]->bottom_right_front->top_left_back!=NULL)||
					(neighbor[2]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,3);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[3]!=NULL)&&(neighbor[3]->top_left_back!=NULL)&&
					((neighbor[3]->top_left_front->top_left_back!=NULL)||
					(neighbor[3]->top_left_back->top_left_back!=NULL)||
					(neighbor[3]->top_right_front->top_left_back!=NULL)||
					(neighbor[3]->top_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,3);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[4]!=NULL)&&(neighbor[4]->top_left_back!=NULL)&&
					((neighbor[4]->top_left_back->top_left_back!=NULL)||
					(neighbor[4]->top_right_back->top_left_back!=NULL)||
					(neighbor[4]->bottom_left_back->top_left_back!=NULL)||
					(neighbor[4]->bottom_right_back->top_left_back!=NULL)))
				{
					splitNode(tempList->node,3);
					delete tempList;
					flag=false;
				}
				else if((flag)&&(neighbor[5]!=NULL)&&(neighbor[5]->top_left_back!=NULL)&&
					((neighbor[5]->top_left_front->top_left_back!=NULL)||
					(neighbor[5]->top_right_front->top_left_back!=NULL)||
					(neighbor[5]->bottom_left_front->top_left_back!=NULL)||
					(neighbor[5]->bottom_right_front->top_left_back!=NULL)))
				{
					splitNode(tempList->node,3);
					delete tempList;
					flag=false;
				}
				
				if (flag)
				{
					//δ�����Ҷ��ֱ��ɾ��
					delete tempList;
				}
			} 
			else
			{
				//δ�����Ҷ��ֱ��ɾ��
				delete tempList;
			}
			//
			flag=true;
			for (int i=0;i<6;i++)
			{
				neighbor[i]=NULL;
			}
			tempList=getLeaf(MyList3);	
			//
		}
	}


}

//////////////////////////////////////////////////////////////////////////
//��������˲���
template <class T>
void preOrder( OctreeNode<T> * & p)
{
	if(p)
	{
		//cout<<i<<".��ǰ�ڵ��ֵΪ��"<<p->data<<"\n����Ϊ��";
		//cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		//cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		//cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		//i+=1;
		//cout<<endl;
		preOrder(p->top_left_front);
		preOrder(p->top_left_back);
		preOrder(p->top_right_front);
		preOrder(p->top_right_back);
		preOrder(p->bottom_left_front);
		preOrder(p->bottom_left_back);
		preOrder(p->bottom_right_front);
		preOrder(p->bottom_right_back);
		//cout<<endl;
	}
	else
	{
		/*leafCount++;*/
	}
}
//////////////////////////////////////////////////////////////////////////
//��˲��������
template<class T>
int depth(OctreeNode<T> *& p)
{
	if(p == NULL)
		return -1;
	int h = depth(p->top_left_front);
	return h+1;
}
//////////////////////////////////////////////////////////////////////////
//���㵥λ���ȣ�Ϊ���ҵ���׼��
int cal(int num)
{
	int result=1;
	if(1==num)
		result=1;
	else
	{
		for(int i=1;i<num;i++)
			result=2*result;
	}
	return result;
}
//////////////////////////////////////////////////////////////////////////
//���ҵ�
//int depth=0;

int times=0;
static double xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,zmax=0;
int tmaxdepth=0;
double txm=1,tym=1,tzm=1;
template<class T>
void find(OctreeNode<T> *& p,double x,double y,double z)
{
	double xm=(p->xmax-p->xmin)/2;
	double ym=(p->ymax-p->ymin)/2;
	double zm=(p->ymax-p->ymin)/2;
	times++;
	if(x>xmax || x<xmin || y>ymax || y<ymin || z>zmax || z<zmin)
	{
		cout<<"�õ㲻�ڳ����У�"<<endl;
		return;
	}
	if(x<=p->xmin+txm && x>=p->xmax-txm && y<=p->ymin+tym && y>=p->ymax-tym && z<=p->zmin+tzm && z>=p->zmax-tzm )
	{
		cout<<endl<<"�ҵ��õ㣡"<<"�õ�λ��"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<"�ڵ��ڣ�"<<endl;
		cout<<"������"<<times<<"�εݹ飡"<<endl;
	}
	else if(x<(p->xmax-xm) && y<(p->ymax-ym) && z<(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->bottom_left_back,x,y,z);
	}
	else if(x<(p->xmax-xm) && y<(p->ymax-ym) && z>(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->top_left_back,x,y,z);
	}
	else if(x>(p->xmax-xm) && y<(p->ymax-ym) && z<(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->bottom_right_back,x,y,z);
	}
	else if(x>(p->xmax-xm) && y<(p->ymax-ym) && z>(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->top_right_back,x,y,z);
	}
	else if(x<(p->xmax-xm) && y>(p->ymax-ym) && z<(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->bottom_left_front,x,y,z);
	}
	else if(x<(p->xmax-xm) && y>(p->ymax-ym) && z>(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->top_left_front,x,y,z);
	}
	else if(x>(p->xmax-xm) && y>(p->ymax-ym) && z<(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->bottom_right_front,x,y,z);
	}
	else if(x>(p->xmax-xm) && y>(p->ymax-ym) && z>(p->zmax-zm))
	{
		cout<<"��ǰ�����ڵ����꣺"<<endl;
		cout<<" xmin: "<<p->xmin<<" xmax: "<<p->xmax;
		cout<<" ymin: "<<p->ymin<<" ymax: "<<p->ymax;
		cout<<" zmin: "<<p->zmin<<" zmax: "<<p->zmax;
		cout<<endl;
		find(p->top_right_front,x,y,z);
	}
}
//////////////////////////////////////////////////////////////////////////
//main����


void process(double *range, double *outData)
{
	
	OctreeNode<double> * rootNode = NULL;

	int choiced = 0;
	bool flag=true;
	while(flag)
	{
		//system("cls");
		//cout<<"��ѡ�������\n";
		//cout<<"1.�����˲��� 2.��������˲���\n";
		//cout<<"3.�鿴����� 4.���ҽڵ�   \n";
		//cout<<"0.�˳�\n\n";
		//cin>>choiced;

		choiced=1;

		if(choiced == 0)
			outData[0]=0;
		else if(choiced == 1)
		{
			//system("cls");
			//cout<<"���������ݹ���ȣ�"<<endl;
			//cin>>maxdepth;
			int depth=0;
			
			//cout<<"��������������꣬˳�����£�xmin,xmax,ymin,ymax,zmin,zmax"<<endl;
			//cin>>xmin>>xmax>>ymin>>ymax>>zmin>>zmax;
			xmin=range[0];
			xmax=range[1];
			ymin=range[2];
			ymax=range[3];
			zmin=range[4];
			zmax=range[5];

			if(depth>=0 || xmax>xmin || ymax>ymin || zmax>zmin || xmin>0 || ymin>0 ||zmin>0)
			{
				/*tmaxdepth=cal(maxdepth);
				txm=(xmax-xmin)/tmaxdepth;
				tym=(ymax-ymin)/tmaxdepth;
				tzm=(zmax-zmin)/tmaxdepth;*/
				createOctree(rootNode,rootNode,depth,xmin,xmax,ymin,ymax,zmin,zmax);
				findLeaf(rootNode);//Ҷ�Ӵ��������б�
				mexPrintf("�����ɵĽڵ���Ϊ %-6d nodes\n",ID-1);
				mexPrintf("�����ɵ�Ҷ����Ϊ %-6d leaves\n",LeafCounter-1);
				balanceTree(rootNode);
				mexPrintf("ƽ�⻯֮���Ҷ��\n");
				
				LeafCounter=1;
				nextList=MyList;
				findLeaf(rootNode);
				mexPrintf("�����ɵĽڵ���Ϊ %-6d nodes\n",ID-1);
				mexPrintf("�����ɵ�Ҷ����Ϊ %-6d leaves\n",LeafCounter-1);
				forFree=rootNode;
				//OctreeNode<double> * leaf[leafCount];
				//preOrder(rootNode);//�洢Ҷ�ӵ�ָ��


			}
			else
			{
				//cout<<"�������";
				mexPrintf("�������\n");
				outData[0]=0;
			}

		}
		else if(choiced == 2)
		{
			system("cls");
			cout<<"��������˲��������\n";
			i=1;
			preOrder(rootNode);
			cout<<endl;
			system("pause");
		}
		else if(choiced == 3)
		{
			system("cls");
			int dep = depth(rootNode);
			cout<<"�˰˲��������Ϊ"<<dep+1<<endl;
			system("pause");
		}
		else if(choiced == 4)
		{
			system("cls");
			cout<<"��������ϣ�����ҵĵ�����꣬˳�����£�x,y,z\n";
			double x,y,z;
			cin>>x>>y>>z;
			times=0;
			cout<<endl<<"��ʼ��Ѱ�õ㡭��"<<endl;
			find(rootNode,x,y,z);
			system("pause");
		}
		else
		{
			system("cls");
			cout<<"\n\n����ѡ��\n";
			system("pause");
		}
		flag=false;
	}
	
}

/*
int _tmain(int argc, _TCHAR* argv[])
{
	OctreeNode<double> * rootNode = NULL;
	int choiced = 0;
	while(true)
	{
		system("cls");
		cout<<"��ѡ�������\n";
		cout<<"1.�����˲��� 2.��������˲���\n";
		cout<<"3.�鿴����� 4.���ҽڵ�   \n";
		cout<<"0.�˳�\n\n";
		cin>>choiced;
		if(choiced == 0)
			return 0;
		else if(choiced == 1)
		{
			system("cls");
			cout<<"���������ݹ���ȣ�"<<endl;
			cin>>maxdepth;
			cout<<"��������������꣬˳�����£�xmin,xmax,ymin,ymax,zmin,zmax"<<endl;
			cin>>xmin>>xmax>>ymin>>ymax>>zmin>>zmax;
			if(maxdepth>=0 || xmax>xmin || ymax>ymin || zmax>zmin || xmin>0 || ymin>0 ||zmin>0)
			{
				tmaxdepth=cal(maxdepth);
				txm=(xmax-xmin)/tmaxdepth;
				tym=(ymax-ymin)/tmaxdepth;
				tzm=(zmax-zmin)/tmaxdepth;
				createOctree(rootNode,maxdepth,xmin,xmax,ymin,ymax,zmin,zmax);
			}
			else
			{
				cout<<"�������";
				return 0;
			}
		}
		else if(choiced == 2)
		{
			system("cls");
			cout<<"��������˲��������\n";
			i=1;
			preOrder(rootNode);
			cout<<endl;
			system("pause");
		}
		else if(choiced == 3)
		{
			system("cls");
			int dep = depth(rootNode);
			cout<<"�˰˲��������Ϊ"<<dep+1<<endl;
			system("pause");
		}
		else if(choiced == 4)
		{
			system("cls");
			cout<<"��������ϣ�����ҵĵ�����꣬˳�����£�x,y,z\n";
			double x,y,z;
			cin>>x>>y>>z;
			times=0;
			cout<<endl<<"��ʼ��Ѱ�õ㡭��"<<endl;
			find(rootNode,x,y,z);
			system("pause");
		}
		else
		{
			system("cls");
			cout<<"\n\n����ѡ��\n";
			system("pause");
		}
	}


	return 0;
}

*/


void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
	//�Ӻ�������

	//const EPS=1.0E-10;
	//const double EPS=1.0E-6;    

	//ȫ�ֱ�����ʼ��
	ID=1;
	i=0;

	MyList->next=NULL;	MyList->node=NULL;	nextList=MyList;
	MyList2->next=NULL;	MyList2->node=NULL;	nextList2=MyList2;
	MyList3->next=NULL;	MyList3->node=NULL;	nextList3=MyList3;
	
	LeafCounter=1;
	LeafCounter2=1;
	LeafCounter3=1;

    double *range,*point;//*threshold; 
	double *outData,*leafInf; 
    int M0,N0,M1,N1; 
    
    //�쳣���� 
    if(nrhs!=2)
        mexErrMsgTxt("USAGE: b=reverse(a)\n"); 
    if(!mxIsDouble(prhs[0])) 
        mexErrMsgTxt("the Input Matrix must be double!\n"); 


    range=mxGetPr(prhs[0]);		//��ķ�Χ 1*6 [xmin xmax ymin ymax zmin zmax]
    point=mxGetPr(prhs[1]);		//����ĵ㼯 3*n :����Ĺ̶���
	M0=mxGetM(prhs[0]); 
	N0=mxGetN(prhs[0]); 
	M1=mxGetM(prhs[1]); 
	N1=mxGetN(prhs[1]); 
	if (N1==3)
		mexErrMsgTxt("Make sure that you have specified the depth of each node in the 4th column!\n"); 	
	
	Point=point;				//ȫ�ֱ�����ֵ
	Size=M1;					//ȫ�ֱ�����ֵ
	pointMark=new OctreeNode<double> * [Size];

    mexPrintf("Range:%d * %d\n",M0,N0);
    mexPrintf("FixPoints:%d * %d\n",M1,N1);

    plhs[0]=mxCreateDoubleMatrix(Size,8,mxREAL); 
    outData=mxGetPr(plhs[0]);

    process(range,outData);

	//�������
	for (int i=0;i<Size;i++)
	{
		if ((pointMark[i]->top_left_back)==NULL)
		{
			outData[i+Size*0]=pointMark[i]->xmin;
			outData[i+Size*1]=pointMark[i]->xmax;
			outData[i+Size*2]=pointMark[i]->ymin;
			outData[i+Size*3]=pointMark[i]->ymax;
			outData[i+Size*4]=pointMark[i]->zmin;
			outData[i+Size*5]=pointMark[i]->zmax;
			outData[i+Size*6]=pointMark[i]->nodedepth;
			outData[i+Size*7]=0;
		}
		else
		{
			outData[i+Size*0]=pointMark[i]->xmin;
			outData[i+Size*1]=pointMark[i]->xmax;
			outData[i+Size*2]=pointMark[i]->ymin;
			outData[i+Size*3]=pointMark[i]->ymax;
			outData[i+Size*4]=pointMark[i]->zmin;
			outData[i+Size*5]=pointMark[i]->zmax;
			outData[i+Size*6]=pointMark[i]->nodedepth;
			outData[i+Size*7]=1;
		}
	}

	delete [ ] pointMark;

	plhs[1]=mxCreateDoubleMatrix(LeafCounter-1,8,mxREAL);
	leafInf=mxGetPr(plhs[1]);

	LeafList * tempLeaf=NULL;
	for (int i=0;i<LeafCounter-1;i++)
	{
		tempLeaf=MyList->next;
		if (tempLeaf==NULL)
		{
			break;
		}
		leafInf[i+(LeafCounter-1)*0]=tempLeaf->node->xmin;		
		leafInf[i+(LeafCounter-1)*1]=tempLeaf->node->xmax;
		leafInf[i+(LeafCounter-1)*2]=tempLeaf->node->ymin;
		leafInf[i+(LeafCounter-1)*3]=tempLeaf->node->ymax;
		leafInf[i+(LeafCounter-1)*4]=tempLeaf->node->zmin;
		leafInf[i+(LeafCounter-1)*5]=tempLeaf->node->zmax;
		leafInf[i+(LeafCounter-1)*6]=tempLeaf->node->nodedepth;
		leafInf[i+(LeafCounter-1)*7]=tempLeaf->node->data;
		MyList->next=tempLeaf->next;
		delete tempLeaf;
	}
	
	

	//delete MyList;		
	//delete MyList2;
	//delete MyList3;
	freeOctree(forFree);
} 

