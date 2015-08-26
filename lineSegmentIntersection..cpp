/*
@Author: Pranay Yadav, 12CS30025
@Description: Computing the Intersection Points and Convex Hull (Chan's Algorithm) resulting from those Intersection Points, from given n- line segments in a 2-D Plane
*/
#include <string>
#include <queue>
#include <iostream>
#include <stdlib.h>
#include <set>
#include <vector>
#include <algorithm>
#define RIGHT_TURN -1  // CW
#define LEFT_TURN 1  // CCW
#define COLLINEAR 0  // COLLINEAR
using namespace std;  
struct point {
	double x;
	double y;
	int segID;
	int seg1ID;
	int seg2ID;
}p0;

typedef struct segment{
	double x1;
	double x2;
	double y1;
	double y2;
	int ID;
	double x;
	double y;
};

/*
	Comparator class for Points to be put in Priority Queue insertion of points based on y coordinates
*/
class CompareEventPoint {
public:
	bool operator()(point& t1, point& t2) // Returns true if t1 is earlier than t2
	{
		if(t1.y<t2.y)return true;
		return false;
	}
};

/*
	Comparator Class for the AVL Tree based Set STL which uses the orientation of each line segment rather than x-coordinates 
*/
struct classComparator {
	bool operator() (const segment& lhs,const segment& rhs) 
	{
		if(rhs.x==lhs.x && rhs.y==lhs.y)
	   		return ((rhs.x1 - lhs.x2) * (rhs.y2 - lhs.y2) - (rhs.y1 - lhs.y2) *     (rhs.x2 - lhs.x2))<0;
	 	else
	   		return ((rhs.x1 - lhs.x1) * (rhs.y2 - lhs.y1) - (rhs.y1 - lhs.y1) * (rhs.x2 - lhs.x1))<0;
	}
};


/*
	Returns square of the distance between the two point class objects
	@param p1: Object of class point aka first point
	@param p2: Object of class point aka second point
*/
int dist(point p1, point p2)
{
	return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}

/*
	Returns orientation of the line joining points p and q and line joining points q and r
	Returns -1 : CW orientation
			+1 : CCW orientation
			0 : Collinear 
	@param p: Object of class point aka first point
	@param q: Object of class point aka second point
	@param r: Object of class point aka third point
*/
int orientation(point p, point q, point r)
{
	int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
	if (val == 0) return 0;  // Collinear
	return (val > 0)? -1: 1; // CW: -1 or CCW: 1
}

/*
	Predicate function used while sorting the points using qsort() inbuilt function in C++
	@param p: Object of class point aka first point
	@param p: Object of class point aka second point
*/
int compare(const void* vp1, const void* vp2)
{
	point* p1 = (point*)vp1;
	point* p2 = (point*)vp2;
	int orient = orientation(p0,*p1,*p2);
	if (orient == 0)
		return (dist(p0, *p2) >= dist(p0, *p1))? -1 : 1;
	return (orient == 1)? -1: 1;
}

/*
	Returns the index of the point to which the tangent is drawn from point p.
	Uses a modified Binary Search Algorithm to yield tangent in O(log n) complexity
	@param v: vector of objects of class points representing the hull aka the vector of hull points
	@param p: Object of class point from where tangent needs to be drawn
*/
int tangent(vector<point> v,point p){
	int l=0;
	int r= v.size();
	int l_before = orientation(p, v[0], v[v.size()-1]);
	int l_after = orientation(p, v[0], v[(l + 1) % v.size()]);
	while (l < r){
		int c = ((l + r)>>1);
		int c_before = orientation(p, v[c], v[(c - 1) % v.size()]);
		int c_after = orientation(p, v[c], v[(c + 1) % v.size()]);
		int c_side = orientation(p, v[l], v[c]);
		if (c_before != RIGHT_TURN and c_after != RIGHT_TURN)
			return c;
		else if ((c_side == LEFT_TURN) and (l_after == RIGHT_TURN or l_before == l_after) or (c_side == RIGHT_TURN and c_before == RIGHT_TURN))
			r = c;
		else
			l = c + 1 ;
		l_before = -c_after; 
		l_after = orientation(p, v[l], v[(l + 1) % v.size()]);
	}
	return l;
}

/*
	Returns the pair of integers representing the Hull # and the point in that Hull which is the extreme amongst all given Hull Points
	@param hulls: Vector containing the hull points for various hulls stored as individual vectors.
*/
pair<int,int> extreme_hullpt_pair(vector<vector<point> >& hulls){
	int h= 0,p= 0;
	for (int i=0; i<hulls.size(); ++i){
		int min_index=0, min_y = hulls[i][0].y;
		for(int j=1; j< hulls[i].size(); ++j){
			if(hulls[i][j].y < min_y){
				min_y=hulls[i][j].y;
				min_index=j;
			}
		}
		if(hulls[i][min_index].y < hulls[h][p].y){
			h=i;
			p=min_index;
		}   
	}
	return make_pair(h,p);
}

/*
	Returns the pair of integers representing the Hull # and the point in that Hull to which the point lpoint will be joined
	@param hulls: Vector containing the hull points for various hulls stored as individual vectors.
	@param lpoint: Pair of the Hull # and the leftmost extreme point contained in that hull, amongst all the obtained hulls
*/
pair<int,int> next_hullpt_pair(vector<vector<point> >& hulls, pair<int,int> lpoint){
	point p = hulls[lpoint.first][lpoint.second];
	pair<int,int> next = make_pair(lpoint.first, (lpoint.second + 1) % hulls[lpoint.first].size());
	for (int h=0; h< hulls.size(); h++){
		if(h != lpoint.first){
			int s= tangent(hulls[h],p);
			point q= hulls[next.first][next.second];
			point r= hulls[h][s];
			int t= orientation(p,q,r);
			if( t== RIGHT_TURN || (t==COLLINEAR) && dist(p,r)>dist(p,q))
				next = make_pair(h,s);
		}
	}
	return next;
}    

/*
	Constraint to find the outermost boundary of the points by checking if the points lie to the left otherwise adding the given point p 
	Returns the Hull Points
	@param v: Vector of all the points 
	@param p: New point p which will be checked to be in the Hull Points or not 
*/
bool isUnequal(point p, point q){
	return (!(p.x==q.x && p.y==q.y && p.segID==q.segID && p.seg1ID==q.seg1ID && p.seg2ID==q.seg2ID));
}	

vector<point> keep_left (vector<point>& v,point p){
	while(v.size()>1 && orientation(v[v.size()-2],v[v.size()-1],p) != LEFT_TURN)
		v.pop_back();
	if(!v.size() || isUnequal(v[v.size()-1],p))
		v.push_back(p);
	return v;
}

/*
	Graham Scan algorithm to find convex hull from the given set of points
	@param points: List of the given points in the cluster (as obtained by Chan's Algorithm grouping)
	Returns the Hull Points in a vector
*/
vector<point> GrahamScan(vector<point>& points){
	if(points.size()<=1)
		return points;
	qsort(&points[0], points.size(), sizeof(point), compare);
	vector<point> lower_hull;
	for(int i=0; i<points.size(); ++i)
		lower_hull = keep_left(lower_hull,points[i]);
	reverse(points.begin(),points.end());
	vector<point> upper_hull;
	for(int i=0; i<points.size(); ++i)
		upper_hull = keep_left(upper_hull,points[i]);
	for(int i=1;i<upper_hull.size();++i)
		lower_hull.push_back(upper_hull[i]);
	return lower_hull;   
}

/*
	Implementation of Chan's Algorithm to compute Convex Hull in O(nlogh) complexity 
*/
vector<point> chansalgorithm(vector<point> v){
	for(int t=0; t< v.size(); ++t){
		for(int m=1; m< (1<<(1<<t)); ++m){
			vector<vector<point> > hulls;
			for(int i=0;i<v.size();i=i+m){
				vector<point> chunk;
				if(v.begin()+i+m <= v.end())
					chunk.assign(v.begin()+i,v.begin()+i+m);
				else
					chunk.assign(v.begin()+i,v.end());        	
				hulls.push_back(GrahamScan(chunk));
			}
			/*cout<<"\nM (Chunk Size): "<<m<<"\n";
			for(int i=0;i<hulls.size();++i){
				cout<<"Convex Hull for Hull #"<<i<<" (Obtained using Graham Scan!!)\n";
				for(int j=0; j<hulls[i].size();++j)
					cout<<hulls[i][j]<<" ";
				cout<<"\n";
			}*/
			vector<pair<int,int> > hull;
			hull.push_back(extreme_hullpt_pair(hulls));
			for(int i=0; i<m; ++i){
				pair<int,int> p= next_hullpt_pair(hulls,hull[hull.size()-1]);
				vector<point> output;
				if(p==hull[0]){
					for(int j=0; j<hull.size();++j){
						output.push_back(hulls[hull[j].first][hull[j].second]);
					}
					return output;
				}
				hull.push_back(p);
			}
		}
	}
}


/*

Earlier implemented classical AVL Tree completely. but it was giving segmentation faults for particular inputs. I have left the code snippets for the
AVLTree implementation however I shall be using Set STL.

class AVLNode{
public:
	segID seg;
	AVLNode* left;
	AVLNode* right;
	int height;
	AVLNode(segID s){
		seg=s;
		left=0;
		right=0;
		height=0;
	}
	bool operator< (const AVLNode* rhs){
		return (seg < rhs->seg);
	}
	bool operator== (const AVLNode* rhs){
		return (seg == rhs->seg && height== rhs->height && left==rhs->left && right==rhs->right);
	}
	bool operator!= (const AVLNode* rhs){
		return (!(seg == rhs->seg && height== rhs->height && left==rhs->left && right==rhs->right));
	}
} *root;

int height(AVLNode* N)
{
    if (N == NULL)
        return 0;
    return N->height;
}

AVLNode* rightRotate(AVLNode* y)
{
	AVLNode* x = y->left;
	AVLNode* T2 = x->right;
	x->right = y;
	y->left = T2;
	y->height = max(height(y->left), height(y->right))+1;
	x->height = max(height(x->left), height(x->right))+1;
	return x;
}
 
AVLNode* leftRotate(AVLNode *x)
{
	AVLNode *y = x->right;
	AVLNode *T2 = y->left;
	y->left = x;
	x->right = T2;
	x->height = max(height(x->left), height(x->right))+1;
	y->height = max(height(y->left), height(y->right))+1;
	return y;
}
 
int getBalance(AVLNode *N)
{
	if (N == NULL)
		return 0;
	return height(N->left) - height(N->right);
}

AVLNode* minValueNode(AVLNode* node)
{
	AVLNode* current = node;
	while (current->left != NULL)
		current = current->left;
	return current;
}

AVLNode* maxValueNode(AVLNode* node)
{
		AVLNode* n = node;
		while(n->right)
			n=n->right;
		return n;
}

pair<AVLNode*,AVLNode*> findNode(AVLNode* root, segID key){
	AVLNode* curr = root, *prev = 0;
	while (curr){
		prev = curr;
		if (key < curr->seg)
			curr = curr->left;
		else if (curr->seg < key)
			curr = curr->right;
		else
			return make_pair(curr, prev);
	}
	return make_pair((AVLNode*)0,prev);
}

AVLNode * inOrderSuccessor(AVLNode *root, AVLNode *n)
{
	if( n->right != NULL )
		return minValueNode(n->right);
	AVLNode *succ = NULL;
	while (root != NULL)
	{
		if (n->seg < root->seg)
		{
			succ = root;
			root = root->left;
		}
		else if (root->seg < n->seg)
			root = root->right;
		else
		   break;
	}
	return succ;
}

AVLNode* inOrderPredecessor(AVLNode* root, AVLNode* x)
{
	if( x->left != NULL )
		return maxValueNode(x->left);
	AVLNode *predecessor = NULL;
	while (root != NULL)
	{
		if (root->seg < x->seg)
		{
			predecessor = root;
			root = root->right;
		}
		else if (x->seg < root->seg)
			root = root->left;
		else
		   break;
	}
	return predecessor;
}

AVLNode* lower_bound(AVLNode* root, segID key){
	if (!root) 
		return 0;
	pair<AVLNode*, AVLNode*> result = findNode(root,key);
	if ((result.first)!=0)
		return result.first;
	if ((result.second)->seg < key)
		return inOrderSuccessor(root,result.second);
	return result.second;
}

AVLNode* upper_bound(AVLNode* root, segID key){
	if (!root) 
		return 0;
	pair<AVLNode*, AVLNode*> result = findNode(root,key);
	if (result.first != (AVLNode*)0)
		return inOrderSuccessor(root,result.first);
	if (result.second->seg < key)
		return inOrderSuccessor(root,result.second);
	return result.second;
}

AVLNode* insert(AVLNode* node, segID key)
{
	if (node == NULL)
		return new AVLNode(key);
	if (key < node->seg)
		node->left  = insert(node->left, key);
	else
		node->right = insert(node->right, key);
	node->height = max(height(node->left), height(node->right)) + 1;
	int balance = getBalance(node);
	// Left Left Case
	if (balance > 1 && key < node->left->seg)
		return rightRotate(node);
	// Right Right Case
	if (balance < -1 && node->right->seg < key)
		return leftRotate(node);
	// Left Right Case
	if (balance > 1 && node->left->seg < key)
	{
		node->left =  leftRotate(node->left);
		return rightRotate(node);
	}
	// Right Left Case
	if (balance < -1 && key < node->right->seg)
	{
		node->right = rightRotate(node->right);
		return leftRotate(node);
	}
return node;
}

AVLNode* erase(AVLNode* root, segID key)
{
	if (root == NULL)
		return root;
	if ( key < root->seg )
		root->left = erase(root->left, key);
	else if( root->seg < key )
		root->right = erase(root->right, key);
	else
	{
		if( (root->left == NULL) || (root->right == NULL) )
		{
			AVLNode *temp = root->left ? root->left : root->right;
			if(temp == NULL)
			{
				temp = root;
				root = NULL;
			}
			else // One child case
			 *root = *temp; // Copy the contents of the non-empty child
			delete temp;
		}
		else
		{
			AVLNode* temp = minValueNode(root->right);
			root->seg = temp->seg;
			root->right = erase(root->right, temp->seg);
		}
	}
	if (root == NULL)
	  return root;
	root->height = max(height(root->left), height(root->right)) + 1;
	int balance = getBalance(root);
	if (balance > 1 && getBalance(root->left) >= 0)
		return rightRotate(root);
	if (balance > 1 && getBalance(root->left) < 0)
	{
		root->left =  leftRotate(root->left);
		return rightRotate(root);
	}
	if (balance < -1 && getBalance(root->right) <= 0)
		return leftRotate(root);
	if (balance < -1 && getBalance(root->right) > 0)
	{
		root->right = rightRotate(root->right);
		return leftRotate(root);
	}
	return root;
}*/

/*
	Checks whether 2 line segments intersect or not and returns the intersection point if yes.
	Uses the method told in lecture i.e. using s and t line equation parameters
*/
point isIntersection(segment a,segment b)
{
	point p;
	p.segID=-1;
	double num=(a.y1-b.y1)*(b.x2-b.x1)-(a.x1-b.x1)*(b.y2-b.y1);
	double den=(a.y1-a.y2)*(b.x2-b.x1)-(a.x1-a.x2)*(b.y2-b.y1);
	double s=num/den;
	if(s>0 && s<1)
	{
		double t=((1-s)*a.x1+s*a.x2-b.x1)/(b.x2-b.x1);
	  	if(t>0 && t<1)
	  	{
	    	p.segID=0;
	    	if(a.x1>b.x1)
	    	{
	      		p.seg1ID=a.ID;
	      		p.seg2ID=b.ID;
	    	}
	    	else
	    	{
	      		p.seg1ID=b.ID;
	      		p.seg2ID=a.ID;    
	    	}
	    	p.x=(1-s)*a.x1+s*a.x2;
	    	p.y=(1-s)*a.y1+s*a.y2;
	    	return p;
	  	}
	}
	return p;
}

/*void printTree(AVLNode* root){
	if(root == (AVLNode*)0)
		return;
	printTree(root->left);
	cout<<"\nNode("<<root->seg.leftend<<" "<<root->seg.rightend<<")\t";
	printTree(root->right);
}*/

/*
	Actual Code where Line Intersections are computed using a priority queue and array of segments taken from user
	Returns list of intersection points out of the given line segments
*/
vector<point> computeIntersectionPoints(priority_queue<point, vector<point>, CompareEventPoint>& pq, segment* segVector)
{
	//AVLNode *lower_it = 0,*upper_it=0;
	set<segment,classComparator> tree;
	set<segment,classComparator>::iterator lower_it,upper_it;
	int segid,iteration=0;
	point temp_p,p;
	double Y_coord_SweepLine;
	vector<point> out;
	while (!pq.empty()) {
		temp_p=pq.top();
		pq.pop(); 
		segid=temp_p.segID;
		Y_coord_SweepLine=temp_p.y;
		iteration++;
		//Left End point event  
		if(segid>0)
		{
		/*	
			Earlier snippet using AVLTrees

			root = insert(root,segVector[segid-1]);
			lower_it= lower_bound (root,segVector[segid-1]);   
			if(lower_it != root)                                                 
			{    
				lower_it= inOrderPredecessor(root,lower_it);
				p=isIntersection(segVector[segid-1],segVector[lower_it->seg.ID-1]);
				if(p.segID==0 && temp_p.y>p.y){
					pq.push(p);
				}
			}
			upper_it=upper_bound (root,segVector[segid-1]); 
			if(upper_it!= ((AVLNode*)0))
			{
				p=isIntersection(segVector[segid-1],segVector[upper_it->seg.ID-1]);
				if(p.segID==0 && temp_p.y>p.y)
					pq.push(p);
			}
		*/	
			tree.insert(segVector[segid-1]);
			lower_it=tree.lower_bound (segVector[segid-1]);   
			if(lower_it!=tree.begin())                                                  //       ^
			{
				lower_it--;
			    p=isIntersection(segVector[segid-1],segVector[lower_it->ID-1]);
			    if(p.segID==0 && temp_p.y>p.y)
			    	pq.push(p);
			}
			upper_it=tree.upper_bound (segVector[segid-1]); 
   			if(upper_it!=tree.end())
   			{
   				p=isIntersection(segVector[segid-1],segVector[upper_it->ID-1]);
     			if(p.segID==0 && temp_p.y>p.y)
       				pq.push(p);
    		}
 		}

 		// Right End Point event
 		else if(segid<0)
   		{
   			segid=-segid;
     		/*			
			Classical AVLTree code snippet

     		if(lower_it!=root && upper_it!= ((AVLNode*)0))
			{
				lower_it=lower_bound (root,segVector[segid-1]);
				lower_it= inOrderPredecessor(root,lower_it);
				upper_it= upper_bound (root,segVector[segid-1]); 
			}
			root = erase(root,segVector[segid-1]);
			if(lower_it!=root && upper_it != ((AVLNode*)0))
			{
				p=isIntersection(segVector[lower_it->seg.ID-1],segVector[upper_it->seg.ID-1]);
				if(p.segID==0 && temp_p.y>p.y)
					pq.push(p);
			}*/
     		if(lower_it!=tree.begin() && upper_it!=tree.end())
     		{
       			lower_it=tree.lower_bound (segVector[segid-1]);
       			lower_it--;
       			upper_it=tree.upper_bound (segVector[segid-1]); 
     		}
     		tree.erase(segVector[segid-1]); 
     		if(lower_it!=tree.begin() && upper_it!=tree.end())
     		{
       			p=isIntersection(segVector[lower_it->ID-1],segVector[upper_it->ID-1]);
       			if(p.segID==0 && temp_p.y>p.y)
         			pq.push(p);
       		}

   		}
 		// Case of intersection Point
 		else
 		{
			int t1=temp_p.seg1ID;
		   	int t2=temp_p.seg2ID;
		   	segment s1=segVector[t1-1];
		   	segment s2=segVector[t2-1];
/*			root = erase(root,s1); 
			root = erase(root,s2);  */
		   	tree.erase(s1); 
		   	tree.erase(s2); 
		   	segVector[t1-1].x=temp_p.x;
		   	segVector[t1-1].x1=temp_p.x;
		   	segVector[t2-1].x=temp_p.x;
		   	segVector[t2-1].x1=temp_p.x;
		   	segVector[t1-1].y=temp_p.y;
		   	segVector[t1-1].y1=temp_p.y;
		   	segVector[t2-1].y=temp_p.y;
		   	segVector[t2-1].y1=temp_p.y;
		   	s1=segVector[t1-1];
		   	s2=segVector[t2-1];
		   	/*
		   	root = insert(root,s1);
			root = insert(root,s2);

			lower_it=lower_bound(root,s1);
			if(lower_it!=root)
			{
				lower_it= inOrderPredecessor(root,lower_it);
				p=isIntersection(s1,segVector[lower_it->seg.ID-1]);
				if(p.segID==0 && temp_p.y>p.y)
					pq.push(p);
			}
			upper_it=upper_bound(root,s2); 
			if(upper_it!= ((AVLNode*)0))
			{
				upper_it= inOrderSuccessor(root,upper_it);
				p=isIntersection(s2,segVector[upper_it->seg.ID-1]);
				if(p.segID==0 && temp_p.y>p.y)
					pq.push(p);
			}*/
		   	tree.insert(s1);
		   	tree.insert(s2);
		   	lower_it=tree.lower_bound (s1);
		   	if(lower_it!=tree.begin())
		    {
		    	lower_it--;
		    	p=isIntersection(s1,segVector[lower_it->ID-1]);
		    	if(p.segID==0 && temp_p.y>p.y)
		        	pq.push(p);
		    }
		    upper_it=tree.upper_bound (s2); 
		    if(upper_it!=tree.end())
		    {
		        p=isIntersection(s2,segVector[upper_it->ID-1]);
			    if(p.segID==0 && temp_p.y>p.y)
		        	pq.push(p);
		    }
			out.push_back(temp_p);
		}
 	}
 	return out;
}

int main()
{
	int N=0,i=0;
	double x1=0,x2=0,y1=0,y2=0;
	priority_queue<point, vector<point>, CompareEventPoint> pq;
	cin>>N;
	segment segVector[N];
	point p;
	for(i=0;i<N;i++)
	{
		cin>>p.x>>p.y;
		p.segID=i+1;
		pq.push(p);
		segVector[i].ID=i+1;
		segVector[i].x1=p.x;
		segVector[i].y1=p.y;
		cin>>p.x>>p.y;
		p.segID=-(i+1);
		pq.push(p);
		segVector[i].x2=p.x;
		segVector[i].y2=p.y;
	}
	vector<point> output= computeIntersectionPoints(pq,segVector);
	cout<<"\n*********************** INTERSECTION POINTS ********************************\n";
	for (int i=0; i< output.size(); ++i){
		cout<<"("<<output[i].x<<","<<output[i].y<<") ";
	}
	cout<<endl;
	vector<point> hull = chansalgorithm(output);
	cout<<"\n**************************** CONVEX HULL *********************************\n";
	for (int i=0; i< hull.size(); ++i){
		cout<<"("<<hull[i].x<<","<<hull[i].y<<") ";
	}
	cout<<endl;
	return 0;
}

