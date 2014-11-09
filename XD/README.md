#DataStructure/SplaySequence.cpp
``` cpp
struct SplaySequence{
  static const int MAXN = 201111;
  struct Node{
    Node *ch[2], *par;
    int s, val, flip;
    Node(){ s = 0; val = -1; }
    int cmp(int k)const{
      int d = k-ch[0]->s;
      if(d==1)return -1;
      else if(d<=0)return 0;
      else if(d>1)return 1;
    }
    void maintain(){ s = ch[0]->s+ch[1]->s+1; }
    void pushDown(){
      if(flip){
        flip=0;
        swap(ch[0], ch[1]);
        ch[0]->reverse();
        ch[1]->reverse();
      }
    }
    Node * & parentPointer(){
      if(par->ch[0]==this)return par->ch[0];
      else return par->ch[1];
    }
    void reverse(){flip^=1;}
  };
  int n;
  Node seq[MAXN], *null, *root;
  SplaySequence(){ null = new Node(); }
  Node* newNode(int val){
    Node &o = seq[n];
    o.s = 1;
    o.flip = 0;
    o.val = val;
    o.ch[0] = null;
    o.ch[1] = null;
    return &seq[n++];
  }
  void rotateR(Node * &o, int d){
    Node *a = o->ch[d^1];
    a->ch[d]->par = o;
    o->ch[d^1] = a->ch[d];
    a->ch[d] = o;
    a->par = o->par;
    o->par = a;
    o->maintain();
    a->maintain();
    o = a;
  }
  void rotate(Node *o, int d){
    o->pushDown();
    if(o==root){
        rotateR(root, d);
    }else{
      int d2 = (o==o->par->ch[1]);
      rotateR(o->par->ch[d2], d);
    }
  }
  void splay (Node *o, Node* &rt){
    o->pushDown();
    while (o != rt){
      o->par->pushDown();
      int d = (o == o->par->ch[1]);
      if (o->par == rt){
        rotate(o->par, d ^ 1);
      }
      else{
        o->par->par->pushDown();
        if (o->par == o->par->par->ch[d]){
          rotate(o->par->par, d ^ 1);
          rotate(o->par, d ^ 1);
        }else{
          rotate(o->par, d ^ 1);
          rotate(o->par, d);
        }
      }
    }
  }
  void splay(Node *o){
    splay(o, root);
  }
  Node * splayK(Node * rt, int k){
    Node *v = rt;
    Node * &rtR = (rt==root?root:rt->parentPointer());
    while(1){
      v->pushDown();
      int d = v->cmp(k);
      if(d==-1){
        splay(v, rtR);
        return v;
      }
      if(d==1)k-= v->ch[0]->s+1;
      v = v->ch[d];
    }
  }
  Node * merge(Node * left, Node * right){// root doesn't make sense
    if(left==null)return right;
    root = left;
    left = splayK(left, left->s);
    left->ch[1] = right;
    right->par = left;
    left->maintain();
    root = null;
    return left;
  }
  void split(int k, Node * &left, Node * &right){// split root
    Node *o = splayK(root, k);

    root = o->par = null;
    left = o;
    right = o->ch[1];
    o->ch[1] = null;
    left->maintain();
  }
  void push_back(Node *o){
    Node *v = root;
    if (root == null){
      root = o;
      o->par = null;
    }
    else{
      while (v->ch[1] != null) v = v->ch[1];
      v->ch[1] = o;
      o->par = v;
      splay(o);
    }
  }
  Node* build(int sz) {
    if(!sz) return null;
    Node* L = build(sz/2);
    Node* o = newNode(n);
    Node* R = build(sz - sz/2 - 1);
    o->ch[0] = L;
    o->ch[1] = R;
    L->par = o;
    R->par = o;
    o->maintain();
    return o;
  }
  void clear(){
    n=0;
    root = null;
  }
  int size(){
    return root->s;
  }
  void sweep (Node *x){
    stack<Node*> stk;
    Node *v = x;
    while(v!=root){
      stk.push(v);
      v = v->par;
    }
    stk.push(root);
    while(!stk.empty()){
      stk.top()->pushDown();
      stk.pop();
    }
  }
  void init(int sz){
    n=0;
    root = build(sz);
  }
};
/*
   10
   2 3 4 5 6 7 8 9 10 1

*/

```
#DataStructure/SegmentTree(kth).cpp
``` cpp
struct pa{
    int val, i, od, sum;
    pa(int _val=0, int _i=0, int _od=0){
        val = _val;
        i = _i;
        od = _od;
    }
};
bool cmp1(const pa &a, const pa &b){
    if(a.val!=b.val)return a.val<b.val;
    return a.i<b.i;
}
bool cmp2(const pa &a, const pa &b){
    return a.i<b.i;
}
struct segment_tree{
    struct segment{
        vector<pa> num;
    }sg[401111];
    int n;
    void init(int _n){
        n=_n;
    }
    void init(vector<int> &ar){
        n = ar.size();
        vector<pa> a(n);
        for(int i=0; i<401111; i++)sg[i].num.clear();
        for(int i=0; i<n; i++){
            a[i].val=ar[i];
            a[i].i=i;
        }
        sort(a.begin(), a.end(), cmp1);
        for(int i=0; i<n; i++)a[i].od=i;
        sort(a.begin(), a.end(), cmp2);
        sg[0].num = a;
        build(0,n,0);
    }
    void build(int l, int r, int k){
        static bool clr[401111];
        vector<pa> &num = sg[k].num;
        if(num.size()<=1)return;
        for(int i=0; i<num.size(); i++){
            if(num[i].od < (l+r)/2){
                sg[k*2+1].num.push_back(num[i]);
                clr[num[i].i]=1;
            }else{
                sg[k*2+2].num.push_back(num[i]);
                clr[num[i].i]=0;
            }
        }
        for(int i=0; i<num.size(); i++){
            num[i].sum = clr[num[i].i];
            if(i)num[i].sum+= num[i-1].sum;
        }
        build(l, (l+r)/2, k*2+1);
        build((l+r)/2, r, k*2+2);
    }
    int kth_element(int l, int r, int k){
        return kth_element(l,r,k,0);
    }
    int kth_element(int l, int r, int k, int sid){
        vector<pa> &num = sg[sid].num;
        if(num.size()<=1)return num[0].val;
        int n_left = num[r-1].sum - (l?num[l-1].sum:0);
        int mid=num.size()/2, n=num.size();
        if(k<n_left){
            int nl = (l?num[l-1].sum:0) ;
            int nr = num[r-1].sum ;
            return kth_element(nl, nr, k, sid*2+1);
        }else{
            int nl = (l-(l?num[l-1].sum:0)) ;
            int nr = (r-num[r-1].sum) ;
            return kth_element(nl, nr, k-n_left, sid*2+2);
        }
    }
}st;

```
#Geometry/geometry.cpp
``` cpp
const double EPS = 1e-8, PI = 3.14159265358979323846;
int dcmp(double x){
        if(x<-EPS)return -1;
        else if(EPS)return 1;
        else return 0;
}
bool between(double x, double a, double b){// closed
        if(b < a)swap(a,b);
        return 0<=dcmp(x-a) && 0<=dcmp(b-x);
}
struct Vector;
struct Line;
struct Circle;
typedef Vector Point;
typedef Line Seg;
struct Vector{
    double x, y;
    Vector(double _x=0, double _y=0){
        this->x = _x;
        this->y = _y;
    }
    Vector operator + (const Vector &that)const{
        return Vector(x+that.x, y+that.y);
    }
    Vector operator - (const Vector &that)const{
        return Vector(x-that.x, y-that.y);
    }
    Vector operator * (double that)const{
        return Vector(x*that, y*that);
    }
    Vector operator / (double that)const{
        return Vector(x/that, y/that);
    }
    double operator ^ (const Vector &that)const{// dot
        return x*that.x + y*that.y;
    }
    double operator * (const Vector &that)const{// cross
        return x*that.y - y*that.x;
    }
    bool operator <(const Vector &that)const{
        return dcmp(x-that.x)<0 || (dcmp(x-that.x)==0 && dcmp(y-that.y)<0);
    }
    bool operator == (const Vector &that)const{
        return dcmp(x-that.x)==0 && dcmp(y-that.y)==0;
    }
    double length()const{
        return sqrt(x*x+y*y);
    }
    double angle()const{//
        return atan2(y,x);
    }
    double interAngle(const Vector &that)const{
        return acos(((*this)^that)/this->length()/that.length());
    }
    double disToPoint(const Point &that)const{
        return (that-*this).length();
    }
    double disToLine(const Line &that)const;
    double disToSeg(const Seg &that)const;
    bool onLine(const Line &that)const;
    bool onSeg(const Seg &that)const;
    Vector projectionToLine(const Line & that)const;
    Vector rotate(double rad)const{
        return Vector(x*cos(rad)-y*sin(rad), x*sin(rad)+y*cos(rad));
    }
    Vector unit()const{
        return (*this)/this->length();
    }
    Vector normal()const{
        return Vector(-y,x).unit();
    }
};
struct Line{
    Point o;
    Vector v;
    Line(){}
    Line(Point _o,Vector _v){
        this->o = _o;
        this->v = _v;
    }
    Vector e()const{
        return o+v;
    }
    bool operator < (const Line &that)const{
        if( dcmp(v*that.v) == 0 )return dcmp(v*(that.o-o)) > 0;
        return dcmp(v*that.v) > 0;
    }
    bool operator == (const Line &that)const{
        return !(*this<that) && !(that<*this);
    }
    Point intersectionWithLine(const Line &that)const;
    vector<Point> intersectionWithCircle(const Circle &that)const;
};
struct Circle{
    Point c;
    double r;
    Circle(Point _c, double _r):c(_c), r(_r){}
    Point point(double ang)const{
        return Point(c.x + r*cos(ang), c.y+r*sin(ang));
    }
    bool operator == (const Circle &that)const{
        return c==that.c && dcmp(r-that.r)==0;
    }
    vector<Point> intersectionWithLine(const Line &that)const;
    vector<Point> intersectionWithCircle(const Circle &that)const;
};
double Point::disToLine(const Line &that)const{
    return fabs((that.v * (*this-that.o)) / that.v.length());
}
double Point::disToSeg(const Seg &that)const{
    double crossa = that.v ^ (*this-that.o);
    double crossb = that.v ^ (*this-that.e());
    if(dcmp(crossa*crossb) < 0)
        return this->disToLine(that);
    else
        return min(this->disToPoint(that.o), this->disToPoint(that.e()));
}
bool Point::onLine(const Line &that)const{
    return dcmp(this->disToLine(that))==0;
}
bool Point::onSeg(const Seg &that)const{
    return dcmp(this->disToSeg(that))==0;
}
Point Point::projectionToLine(const Line &that)const{
    const Point &o = that.o;
    const Vector &v = that.v;
    return o + v*( (v^(*this-o))/(v^v));
}
Point Line::intersectionWithLine(const Line &that)const{
    Vector u = o-that.o;
    double t = (that.v*u) / (v*that.v);
    return o+v*t;
}
vector<Point> Line::intersectionWithCircle(const Circle &that)const{
    vector<Point> res;
    Point c = that.c, p;
    double d = c.disToLine(*this), r = that.r;
    if(dcmp(r-d)>=0){
        p = c.projectionToLine(*this);
        double l = sqrt(r*r-d*d);
        if(dcmp(r-d)==0){
            res.push_back(p);
        }else{
            res.push_back(p+v.unit()*l);
            res.push_back(p-v.unit()*l);
        }
    }
    return res;
}
vector<Point> Circle::intersectionWithLine(const Line &that)const{
    return that.intersectionWithCircle(*this);
}
vector<Point> Circle::intersectionWithCircle(const Circle &that)const{
    vector<Point> res;
    double d = (c-that.c).length();
    if(dcmp(d)==0);
    else if(dcmp(r+that.r-d)<0);
    else if(dcmp(fabs(r-that.r)-d)>0);
    else{
        double ang = (that.c-c).angle();
        double da = acos((r*r+d*d-that.r*that.r)/(2*r*d));
        Point p1 = point(ang-da), p2 = point(ang+da);
        res.push_back(p1);
        if(!(p1==p2))
            res.push_back(p2);
    }
    return res;
}
double triangleArea(const Point &a, const Point &b, const Point &c){
    Vector va = a-c, vb = b-c;
    return abs(va*vb)/2;
}
vector<Point> halfplane(vector<Line> l)
{
    sort(l.begin(),l.end());
    l.resize(int(unique(l.begin(),l.end())-l.begin()));
    vector<Line> q;
    q.resize(int(l.size()));
    int L=0,R=0;
    for( int i = 0; i < int(l.size()); i++ ){
        while( R-L >= 2 ){
            Point p = q[R-1].intersectionWithLine(q[R-2]);
            if( dcmp((l[i].v)*(p-l[i].o)) <= 0 )R--;
        }
        while( R-L >= 2 ){
            Point p = q[L].intersectionWithLine(q[L+1]);
            if( dcmp((l[i].v)*(p-l[i].o)) <= 0 )L++;
        }
        q[R++] = l[i];
    }
    while( R-L >= 2 ){
        Point p = q[R-1].intersectionWithLine(q[R-2]);
        if( dcmp((q[L].v)*(p-q[L].o)) <= 0 )R--;
    }
    vector<Point> res;
    if( R-L < 3 )return res;
    for( int i = L; i < R-1; i++ )
        res.push_back(q[i].intersectionWithLine(q[i+1]));
    res.push_back(q[L].intersectionWithLine(q[R-1]));
    return res;
}

```
#Graph/BCC.cpp
``` cpp
// Knights of the Round Table
typedef vector<vector<int> > Graph;
Graph G, sG;
int n, m;
struct Pair{
	int a, b;
	Pair(int _a=0, int _b=0): a(_a), b(_b){}
	bool operator !=(const Pair &that){
		return a!=that.a || b!=that.b;
	}
};
bool has[1111], vst[1111];
int dep[1111], back[1111], clr[1111];
stack<Pair> stk;

bool findOddCycle(Graph &G, int v, int c, int clr[1111]){
	if(clr[v]!=-1){
		if(clr[v]==c)return 0;
		else return 1;
	}
	clr[v]=c;
	for(int i=0; i<G[v].size(); i++){
		int u=G[v][i];
		if(findOddCycle(G, u, c^1, clr))return 1;
	}
	return 0;
}
void addEdge(Graph &G, const Pair &pr){
	G[pr.a].push_back(pr.b);
	G[pr.b].push_back(pr.a);
}
void dfs(int v, int p, int d){
	vst[v]=1;
	dep[v] = back[v] = d;
	for(int i=0; i<G[v].size(); i++){
		int u=G[v][i];
		if(u==p|| ( vst[u] && dep[u]>=dep[v] ))continue;
		Pair e(v,u);
		stk.push(e);
		if(vst[u]){// back edge
			back[v] = min(back[v], dep[u]);
			continue;
		}
		dfs(u, v, d+1);
		back[v] = min(back[v], back[u]);
		if(back[u]>=dep[v]){
			Pair pr;
			sG.clear(), sG.resize(n);
			while(!stk.empty() && stk.top()!=e){
				addEdge(sG, stk.top());
				stk.pop();
			}
			addEdge(sG, stk.top());
			stk.pop();
			memset(clr, -1, sizeof(clr));
			if(findOddCycle(sG, v, 0, clr)){
				for(int j=0; j<n; j++)
					if(sG[j].size()>0)has[j]=1;
			}
		}
	}
}
```
#Graph/MinCut.cpp
``` cpp
int V, G[MAX_V][MAX_V], w[MAX_V], contract[MAX_V];
bool vst[MAX_V];

void MAS(int &s, int &t, int &c, int cmbt){//maximum_adjacency_search
    memset(w,0,sizeof(w));
    memset(vst,0,sizeof(vst));

    for(int cnt=0; cnt<V-cmbt; cnt++){
        int x=0, maxw=0;
        for(int v=0; v<V; v++)
            if(contract[v]==-1 && !vst[v] && maxw<w[v]){
                x=v; maxw=w[v];
            }
        vst[x]=1;
        if(cnt==V-cmbt-2)s=x;
        if(cnt==V-cmbt-1){t=x; c=w[t];}

        for(int v=0; v<V; v++)
            if(contract[v]==-1 && !vst[v])
                w[v]+= G[x][v];
    }
}
int mincut(){//Stoer-Wagner O(V^3)
    memset(contract,-1,sizeof(contract));
    int s, t, c, res=INF;
    for(int cnt=0; cnt<V-1; cnt++){
        MAS(s, t, c, cnt);
        res = min(res, c);
        //t is combined to s
        contract[t]=s;
        for(int v=0; v<V; v++)
            if(contract[v]==t) contract[v]=s;
        for(int v=0; v<V; v++)
            if(contract[v]==-1){
                G[s][v]+= G[t][v];
                G[v][s]+= G[v][t];
            }
    }
    return res;
}
```
#Graph/GeneralMatching.cpp
``` cpp
int inQ[MAXN], frm[MAXN], typ[MAXN], wife[MAXN], fa[MAXN], vit[MAXN];
deque<int> Q;
vector<int> E[MAXN];
int rt(int u){ if(fa[u]==u) return u; return fa[u] = rt(fa[u]);}
void unionRt(int u, int v){ if(rt(u)!=rt(v)) fa[rt(u)] = rt(v); }
int lca(int u, int v){
  memset(vit, 0, sizeof(vit));
  while( u ) vit[ rt(u) ] = 1, u = frm[wife[rt(u)]];
  while( v ){
    if(vit[rt(v)]) return rt(v);
    v = frm[wife[rt(v)]];
  }
}
void contractPath(int u, int anc){
  for(int wf, fm; u != anc; u = fm){
    wf = wife[u], fm = frm[wf];
    if(rt(fm) != anc) frm[fm] = wf;
    if(typ[wf]==2) Q.pb(wf), typ[wf] = inQ[wf] = 1;
    if(typ[fm]==2) Q.pb(fm), typ[fm] = inQ[fm] = 1;
    unionRt(u, wf), unionRt(wf, fm);
  }
}
void contract(int u, int v){
  int anc = lca(u, v);
  if(rt(u)!=anc) frm[u] = v;
  if(rt(v)!=anc) frm[v] = u;
  contractPath(u, anc);
  contractPath(v, anc);
}
void path(int st){
  Q.clear();
  for(int i=0;i<=n;i++) fa[i] = i, inQ[i] = typ[i] = frm[i] = 0;
  typ[st] = inQ[st] = 1; Q.pb(st);
  while(sz(Q)){
    int u = Q.front(); Q.pop_front(); inQ[u] = 0;
    for(int i=0;i<sz(E[u]);i++){
      int v = E[u][i];
      if(v != wife[u] && rt(u)!=rt(v) && typ[v]!=2){
        if(typ[v]==1) contract(u, v);
        else if(!wife[v]){
          frm[v] = u;
          for(int tv = v, tu = frm[tv], chW; tu && tv; tv = chW, tu = frm[tv]){
            chW = wife[tu];
            wife[tu] = tv, wife[tv] = tu;
            wife[chW] = 0;
          }
          return ;
        }else{
          frm[v] = u, typ[v] = 2;
          inQ[wife[v]] = typ[wife[v]] = 1; Q.pb(wife[v]);
        }
      }
    }
  }
}
void Edmond(){
  for(int i=1;i<=n;i++) E[i].clear(), wife[i] = 0;
  for(int i=1;i<=n;i++) if(!wife[i]) path(i);
}

```
#Graph/MincostFlow.cpp
``` cpp
struct edge{
	int to, cap, rev, cost;
	edge(int _to, int _cap, int _rev, int _cost): to(_to), cap(_cap), rev(_rev), cost(_cost) {}
};
vector<edge> G[MAXV]
void add_edge(int v, int u, int cap, int cost){
	G[v].push_back(edge(u, cap,   G[u].size(),  cost));
	G[u].push_back(edge(v,   0, G[v].size()-1, -cost));
}
int dis[MAXV], prev[MAXV], pree[MAXV];
bool inq[MAXV];

void sssp(int s){
	queue<int> q;
	fill(dis, dis+MAXV, INF);
	memset(inq,0,sizeof(inq));
	dis[s]=0;
	q.push(s), inq[s]=1;
	while(!q.empty()){
		int v=q.front();
		q.pop(), inq[v]=0;
		for(int i=0; i<G[v].size(); i++){
			edge &e=G[v][i];
			if(0<e.cap && dis[v]+e.cost<dis[e.to]){
				dis[e.to] = dis[v]+e.cost;
				prev[e.to] = v;
				pree[e.to] = i;
				if(!inq[e.to])
					q.push(e.to), inq[e.to]=1;
			}
		}
	}
}
int mincost_flow(int s, int t, int flow){
	int cost=0;
	while(0<flow){
		sssp(s);
		if(dis[t]==INF)return -1;
		int f=flow;
		for(int v=t; v!=s; v=prev[v]){
			edge &e=G[prev[v]][pree[v]];
			f = min(f, e.cap);
		}
		for(int v=t; v!=s; v=prev[v]){
			edge &e  = G[prev[v]][pree[v]];
			edge &re = G[  e.to ][ e.rev ];
			e.cap -= f;
			re.cap+= f;
		}
		cost+= f*dis[t];
		flow-= f;
	}
	return cost;
}

```
#Graph/HeavyLight.cpp
``` cpp
const int V = 100000;
vector<int> adj[V];
int parent[V], heavy[V];
int depth[V], size[V]; 
int chain[V], head[V]; 
 
void dfs(int i){
    size[i] = 1;
    for (int k=0; k<adj[i].size(); ++k){
        int j = adj[i][k];
        if (j == parent[i]) continue;
        parent[j] = i;
        depth[j] = depth[i] + 1;
        DFS(j);
        size[i] += size[j];
        if (heavy[i] == -1 || size[j] > size[heavy[i]])
            heavy[i] = j;
    }
}
void heavylight(int N){
    memset(heavy, -1, sizeof(heavy));
    parent[0] = -1;
    depth[0] = 0;
    DFS(0);
    int c = 0;  // 鏈的編號
    for (int i=0; i<N; ++i)
        if (parent[i] == -1 || heavy[parent[i]] != i){
            // i點是鏈的開頭祖先
            for (int k = i; k != -1; k = heavy[k])
                chain[k] = c, head[k] = i;
            c++;
        }
}
int lca(int i, int j){
    // 深的上升
    while (chain[i] != chain[j])
        if (depth[head[i]] > depth[head[j]]) i = parent[head[i]];
        else j = parent[head[j]];
    return depth[i] < depth[j] ? i : j;
}
```
#Graph/SCC.cpp
``` cpp
#define MAX_V 10000
int V;
vector<int> G[MAX_V];
vector<int> rG[MAX_V];
vector<int> vs;
bool used[MAX_V];
int cmp[MAX_V];
void add_edge(int a, int b){
    G[a].push_back(b);
    rG[b].push_back(a);
}
void dfs(int v){
    used[v] =true;
    for(int i=0; i<int(G[v].size()); i++){
        if(!used[G[v][i]])
        dfs(G[v][i]);
    }
    vs.push_back(v);
}
void rdfs(int v, int k){
    used[v] = true;
    cmp[v] = k;
    for(int i=0; i<int(rG[v].size()); i++)
        if(!used[rG[v][i]])
            rdfs(rG[v][i], k);
}
int scc(){
    memset(used, 0, sizeof(used));
    vs.clear();
    for(int v=0; v<V; v++)
        if(!used[v])
            dfs(v);
    memset(used, 0, sizeof(used));
    int k=0;
    for(int i=int(vs.size())-1; 0<=i; i--)
        if(!used[vs[i]])
            rdfs(vs[i], k++);
    return k;
}
```
#Graph/MinimumDirectedSpanningTree.cpp
``` cpp
map<string, int> lib;
string name[110];
int courN, m;
int edge[200000][3];
int minIn[110], frm[110], leader[110], del[110], use[110];
void Init(){lib.clear();courN = 1;m = 0;}
int getID(string &str){
  if(lib.count(str)!=1) {
    name[courN] = str;
    lib[str] = courN++;
  }
  return lib[str];
}
void add_edge(int u, int v, int cst){
  edge[m][0] = u;edge[m][1] = v;edge[m][2] = cst;
  m++;
}
int main()
{
  int tn,n;
  int i,j;
  string str;
  int a,b,t;
  int tid, id;
  int costNode, costEdge;
  int inCirV;
  name[0] = string("Source");

  cin >> tn;
  while(tn--){
    cin >> n;
    Init();
    for(i=0;i<n;i++){
      cin >> str;
      scanf("%d%d",&a,&b);
      id = getID(str);
      add_edge(0, id, a);
      while(b--){
        cin >> str;
        scanf("%d",&a);
        tid = getID(str);
        add_edge(tid, id, a);
      }
    }
    cin >> t;
    while(t--){
      cin >> str;
      id = getID(str);
      add_edge(0, id, 0);
    }
    n++;
    memset(del, 0, sizeof(del));
    costNode = 0;
    while(1){
      for(i=0;i<n;i++) minIn[i] = INF, frm[i] = -1;
      for(i=0;i<m;i++){
        int u = edge[i][0];
        int v = edge[i][1];
        int cst = edge[i][2];
        if(u==v)continue;
        if(cst < minIn[v])
          minIn[v] = cst, frm[v] = u;
      }
      bool circleFlag = false;
      costEdge = 0;
      memset(use,-1,sizeof(use));
      memset(leader,-1,sizeof(leader));
      for(i=0;i<n;i++){
        if(del[i])continue;
        if(frm[i]==-1 && i!=0)return INF;
        if(i)costEdge += minIn[i];
        for(j=i; j!=-1 && use[j]==-1 ;j = frm[j])use[j] = i;
        if(j!=-1 && use[j]==i){
          circleFlag = true; inCirV = j;
          do{
            del[j] = 1;
            leader[j] = inCirV;
            costNode += minIn[j];
            j = frm[j];
          }while(j != inCirV);
          del[inCirV] = 0;
        }
      }
      if(!circleFlag)break;
      for(i=0;i<m;i++){
        int &u = edge[i][0];
        int &v = edge[i][1];
        int &cst = edge[i][2];
        if(u==v)continue;
        if(leader[v]!=-1)cst -= minIn[v];
        if(leader[u]!=-1)u = leader[u];
        if(leader[v]!=-1)v = leader[v];
        if(u == v) u = v = -1;
      }
    }
    if(debug)cout << "test:" << (costEdge+costNode) << "\n";
    cout << int((costNode+costEdge+4)/5) << "\n";
  }
  return 0;
}

```
#Graph/MaxFlow.cpp
``` cpp
int lvl[MAX_V],iter[MAX_V],V;
struct edge{int to, cap, f, rev;};
vector<edge> G[MAX_V];

void add_edge(int v, int u, int cap){
    G[v].push_back(edge{u, cap, 0,G[u].size()});
    G[u].push_back(edge{v, 0, 0,G[v].size()-1});
}
void BFS(int s){
    fill(lvl, lvl+V, INF);
    queue<int> q;
    q.push(s);
    while(!q.empty()){
        int v=q.front(); q.pop();
        for(int i=0; i<G[v].size(); i++){
            edge e=G[v][i];
            if(0<(e.cap-e.f) && lvl[e.to]==INF){
                lvl[e.to]=lvl[v]+1;
                q.push(e.to);
            }
        }
    }
}
int DFS(int v, int t, int f){
    if(v==t)return f;
    for(int &i=iter[v]; i<G[v].size(); i++){
        edge &e=G[v][i];
        if(0<(e.cap-e.f) && lvl[v]<lvl[e.to]){
            int d = DFS(e.to, t, min(f,e.cap-e.f));
            if(0<d){
                e.f+= d;
                G[e.to][e.rev].f-=d;
                return d;
            }
        }
    }
    return 0;
}
int maxFlow(int s, int t){//Dinic
    int flow=0,f;
    while(1){
        BFS(s);
        if(lvl[t]==INF)return flow;
        memset(iter,0,sizeof(iter));
        while( 0<(f=DFS(s,t,INF)) ) flow+= f;
    }
}

```
#Graph/MaximumClique.cpp
``` cpp
int n;
int G[MAXN][MAXN];
int dp[MAXN], best;

bool dfs(int id[],int top,int cnt)
{
    if(!top){
        if(best<cnt){
            best = cnt;
            return true;
        }
        return false;
    }
    int a[MAXN];
    for(int i=0;i<top;i++){
        if(cnt+top-i<=best) return false;
        if(cnt+dp[id[i]]<=best) return false;
        int k = 0;
        for(int j=i+1;j<top;j++)
            if(G[id[i]][id[j]])
                a[k++] = id[j];
        if(dfs(a,k,cnt+1)) return true;
    }
    return false;
}

int solve(void)
{
    int id[MAXN];
    best = 0;
    for(int i=n-1;i>=0;i--){
        int top = 0;
        for(int j=i+1;j<n;j++)
            if(G[i][j])
                id[top++] = j;
        dfs(id,top,1);
        dp[i] = best;
    }
    return best;
}
```
#Number/Number.cpp
``` cpp
LL powMod(LL a, LL p, LL m){
	if(p==0)return 1%m;
	LL res=powMod(a, p/2, m);
	res = res*res%m;
	if(p%2==1)
		res = res*a%m;
	return res;
}
LL gcd(LL a, LL b){
	if(b==0)return a;
	return gcd(b, a%b);
}
// |x|+|y| min
LL extGcd(LL a, LL b, LL &x, LL &y){
	LL res;
	if(b==0){
		x=1;
		y=0;
		res = a;
	}else{
		LL res = extGcd(b, a%b, y, x);
		y-= x*(a/b);
	}
	return res;
}
LL invMod(LL a, LL m){
	LL x, y;
	extGcd(a, m, x, y);
	return x;
}
int eulerPhi(int n){
	int m = sqrt(n+0.5);
	int res=n;
	for(int i=2; i<=m; i++){
		if(n%i==0){
			res = res*(i-1)/i;
			while(res%i==0)res/=i;
		}
	}
	if(n>1) res = res*(n-1)/n;
}
vector<int> phiTable(int n){
	vector<int> phi(n, 0);
	phi[1] = 1;
	for(int i=2; i<=n; i++){
		for(int j=i; j<=n; j++){
			if(!phi[j])phi[j] = j;
			phi[j] = phi[j]*(i-1)/i;
		}
	}
	return phi;
}
// a[i]*x = b[i](mod m[i])
// return pair(x, M)
typedef vector<int> vI;
pair<int,int> linearCongruence(const vI &a, const vI &b, const vI &m){
	// initially x=0(mod 1)
	int x=0, M=1;
	for(int i=0; i<a.size(); i++){
		int a_ = a[i]*M, b_=b[i]-a[i]*x, d=gcd(m[i],a_);
		if(b_%d!=0)return make_pair(0,-1);
		int t = b_ / d * invMod(a_/d, m[i]/d) %(m[i]/d);
		x = x+M*t;
		M*= m[i]/d;
	}
	return make_pair(x%M, M);
}
// x = a[i] (mod m[i]), a coprime each other
LL china(const vector<int> &a, const vector<int> &m){
	LL M=1, d, res=0, x, y;
	for(int i=0; i<m.size(); i++)M*= m[i];
	for(int i=0; i<a.size(); i++){
		LL w = M/m[i];
		extGcd(m[i], w, x, y);
		res = (res + w*y*a[i])%M;
	}
	return res;
}
// a^x = b (mod p), n is prime
int logMod(int a, int b, int p){
	int m=sqrt(p+0.5), v, e=1;
	v = invMod(powMod(a, m, p), p);
	map<int,int> mp;
	mp[1] = 0;
	for(int i=1; i<m ;i++){
		e = (e*a)%p;
		if(!mp.count(e))mp[e] = i;
	}
	for(int i=0; i<m; i++){
		if(mp.count(b))return i*m + mp[b];
		b = (b*v)%p;
	}
	return -1;
}
// n! (mod p) n! = a*p^e (p is prime)
// based on Wilson's theorem
const MAX_P = 101;
int fact[MAX_P];
int factMod(int n, int p, int &e){
	e = 0;
	if(n==0)return 1;
	int res = factMod(n/p, p, e);
	e+= n/p;
	if(n / p  %2 != 0)return res*(p-fact[n%p])%p;
	return res*fact[n%p]%p;
}
int combMod(int n, int m, int p){
	if(n<0 || k<0 || n<k)return 0;
	int e1, e2, e3;
	int a1 = factMod(n, p, e1), a2 = factMod(k, p, e2), a3 = factMod(n-k, p, e3);
	if(e1 > e2+e2)return 0;
	return a1* invMod(a2*a3%p, p) %p;
}
// Miller-Rabin
// fake prime test always correct when n<2^32
// if n < 4,759,123,141, it is enough to test a = 2, 7, and 61;
// if n < 3,825,123,056,546,413,051, it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, and 23.

bool MRTest(LL x, LL n){
	LL y = n-1;
	while(y%2==0)y/=2;
	LL z = powMod(x, y, n);
	if(z==1)return 1;
	while(y<n-1 && z!=1 && z!=n-1)
		z=(z*z)%n, y*=2;
	return z==n-1;
}
bool isPrime(LL n){
	if(n==2 || n==7 || n==61)return 1;
	if(n==1 || n%2==0) return 0;
	return MRTest(2, n) && MRTest(7, n) && MRTest(61, n);
}
//integration
double simpson(double (*f)(double), double a, double b){
	return (f(a)+4*f((a+b)/2)+f(b))*(b-a)/6;
}
double adaptiveSimpson(double (*f)(double), double a, double b, double eps){
	double c=(a+b)/2;
	double A=simpson(f,a,b), L=simpson(f,a,c), R=simpson(f,c,b);
	if(fabs(L+R-A)<=15*eps) return L+R+(L+R-A)/15;
	return adaptiveSimpson(f,a,c,eps/2) + adaptiveSimpson(f,c,b,eps/2);
}
//moebius function
map<int,int> moebius(int n){
	map<int,int> res;
	vector<int> primes;
	for(int i=2; i*i<=n; i++){
		if(n%i==0){
			primes.push_back(i);
			while(n%i==0)n/=i;
		}
	}
	if(n!=1) primes.push_back(n);
	int m = primes.size();
	for(int i=0; i<(1<<m); i++){
		int mu=1, d=1;
		for(int j=0; j<m; j++){
			if(i>>j &1){
				mu*= -1;
				d*= primes[j];
			}
		}
		res[d] = mu;
	}
	return res;
}

```
#String/ACautomata.cpp
``` cpp
const int SIGMA_SIZE = 26;
const int MAXNODE = 11000;
const int MAXS = 150 + 10;

struct AhoCorasickAutomata {
  int ch[MAXNODE][SIGMA_SIZE];
  int fail[MAXNODE];
  int val[MAXNODE];  
  int last[MAXNODE]; 
  int cnt[MAXS];
  int sz;

  void init() {
    sz = 1;
    memset(ch[0], 0, sizeof(ch[0]));
    memset(cnt, 0, sizeof(cnt));
  }

  int idx(char c) {return c-'a';}

  void insert(char *s, int v) {
    int u = 0, n = strlen(s);
    for(int i = 0; i < n; i++) {
      int c = idx(s[i]);
      if(!ch[u][c]) {
        memset(ch[sz], 0, sizeof(ch[sz]));
        val[sz] = 0;
        ch[u][c] = sz++;
      }
      u = ch[u][c];
    }
    val[u] = v;
  }

  void print(int j) {
    if(j) {
      printf("%d\n",val[j]);
      print(last[j]);
    }
  }

  int find(char* T) {
    int n = strlen(T);
    int j = 0;
    for(int i = 0; i < n; i++) {
      int c = idx(T[i]);
      while(j && !ch[j][c]) j = fail[j]; 
      j = ch[j][c];
      if(val[j]) print(j);
      else if(last[j]) print(last[j]); 
    }
  }

  void getFail() {
    queue<int> q;
    fail[0] = 0;
    for(int c = 0; c < SIGMA_SIZE; c++) {
      int u = ch[0][c];
      if(u) { fail[u] = 0; q.push(u); last[u] = 0; }
    }
    while(!q.empty()) {
      int r = q.front(); q.pop();
      for(int c = 0; c < SIGMA_SIZE; c++) {
        int u = ch[r][c];
        if(!u) continue;
        q.push(u);
        int v = fail[r];
        while(v && !ch[v][c]) v = fail[v];
        fail[u] = ch[v][c];
        last[u] = val[fail[u]] ? fail[u] : last[fail[u]];
      }
    }
  }
};


```
#String/Manacher.cpp
``` cpp
void manacher(int n, const char s[], int p[]){
  for (int i=0, j=0, k=0; i <= 2*(n-1); i++){
    int l = i<k ? min(p[j+j-i], (k-i)/2) : 0;
    int a = i/2-l, b = (i+1)/2+l;
    
    while(0 <= a && b < n && s[a] == s[b]){
      --a;
      ++b;
      ++l;
    }
    
    p[i] = l;
    
    if(k < 2 * b){
      j = i;
      k = 2 * b;
    }
  }
}
```
#String/MorrisPratt.cpp
``` cpp
int mp(int tpl[], int n, int pat[], int m){
	if(n<m)return 0;
	int res=0;
	fail[0]=-1;
	for(int i=1,j=-1; i<m; i++){
		while(0<=j && pat[j+1]!=pat[i]) j=fail[j];
		if(pat[i]==pat[j+1])j++;
		fail[i] = j;
	}
	for(int i=0,j=-1; i<n; i++){
		while(0<=j && tpl[i]!=pat[j+1]) j=fail[j];
		if(pat[j+1]==tpl[i])j++;
		if(m-1<=j){
			res++;
			j=fail[j];
		}
	}
	return res;
}
```
#String/SA+LCPA.cpp
``` cpp
#define MAX_L 1111
int v[MAX_L], prev[MAX_L], sa[MAX_L], rsa[MAX_L], lcpa[MAX_L];
bool cmp(int a, int b){ return v[a]<v[b]; }
void cptSAandLCPA(char s[], int sa[], int lcpa[]){
    int len=strlen(s);
    //SA
    for(int i=0; i<len; i++)sa[i]=i;
    for(int p=0,sh=1; sh<len; p++,sh*=2){
        if(p==0)
            for(int i=0; i<len; i++)v[i]=s[i];
        else{
            for(int i=0; i<len; i++){
                v[i]= prev[i]*MAX_L;
                if(i+sh/2<len)v[i]+= prev[i+sh/2];
            }
        }
        sort(sa, sa+len, cmp);
        for(int i=0; i<len; i++){
            prev[sa[i]] = i;
            if(0<i) if(v[sa[i-1]]==v[sa[i]]) prev[sa[i]] = prev[sa[i-1]];
        }
    }
    for(int i=0; i<len; i++)rsa[sa[i]]=i;
    // LCPA
    for(int lcp=0,i=0; i<len;i++){
        if(rsa[i]==0)
            lcpa[rsa[i]]=0;
        else{
            if(lcp>0)lcp--;
            int j=sa[rsa[i]-1];// behind i in SA order
            while(s[i+lcp] == s[j+lcp])lcp++;
            lcpa[rsa[i]]=lcp;
        }
    }
}

```
#String/Z.cpp
``` cpp
void makeZ( int z[], char s[], int l ){
	z[0] = l;
	int L=0, R=0, i, x;
	for( i = 1 ; i < l ; i++ ){
		if( R < i || z[i-L] >= R-i+1 ){
			R < i ? x = i : x = R+1;
			while( x < l && s[x] == s[x-i] ) x++;
			z[i] = x-i;
			if( i < x ){ L = i; R = x-1; }
		}
		else z[i] = z[i-L];
	}
}

```
