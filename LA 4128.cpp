#include <iostream>
#include <cstdio>
#include <vector>
#include <queue>
#include <cstring>
using namespace std;

const int dx[]={1,0,-1,0}, dy[]={0,1,0,-1};
const int INF = 1011110000;

int n, m, sx,sy,ex,ey,vir[111][111], hor[111][111];

struct edge{
	int to, cost;
	edge(int _to=0, int _cost=0):to(_to), cost(_cost){}
};

vector<edge> G[111*111*5];

void add_edge(int v, int u, int cost){// v->u
	G[v].push_back(edge(u, cost));
	// printf("%d -> %d : %d\n",v,u,cost);
}
int idx(int x, int y, int d){
	return x*m*5 + y*5 + d;
}
bool inside(int x, int y){
	return 0<=x && x<n && 0<=y && y<m;
}
int edgeCost(int x, int y, int d){
	int nx=x+dx[d], ny=y+dy[d];
	if(!inside(nx,ny))return 0;
	if(d==0)return vir[x][y];
	if(d==1)return hor[x][y];
	if(d==2)return vir[x-1][y];
	if(d==3)return hor[x][y-1];
	return 0;
}

int dis[111*111*5];
bool inq[111*111*5];
bool input(){
	scanf("%d%d%d%d%d%d",&n,&m,&sx,&sy,&ex,&ey);
	if(n+m+sx+sy+ex+ey==0)return 0;
	for(int i=0; i<n-1; i++){
		for(int j=0; j<m-1; j++)scanf("%d",&hor[i][j]);// ->
		for(int j=0; j<m; j++)scanf("%d",&vir[i][j]);// \|/
	}
	for(int j=0; j<m-1; j++)scanf("%d",&hor[n-1][j]);// ->
	sx--,sy--,ex--,ey--;
	return 1;
}

void solve(){
	for(int i=0; i<111*111*5; i++)G[i].clear();
	for(int x=0; x<n; x++){
		for(int y=0; y<m; y++){
			for(int d=0; d<4; d++){
				int nx=x+dx[d], ny=y+dy[d], cost = edgeCost(x,y,d);
				if(0<cost){
					add_edge(idx(x,y,d), idx(nx,ny,d), cost);// straight
					add_edge(idx(x,y,d), idx(nx,ny,4), 2*cost);// before turn
				}
			}
			for(int d=0; d<4; d++){
				int nx=x+dx[d], ny=y+dy[d], cost = edgeCost(x,y,d);
				if(0<cost){
					add_edge(idx(x,y,4), idx(nx,ny,d), 2*cost);// after turn
					add_edge(idx(x,y,4), idx(nx,ny,4), 2*cost);// turn after turn
				}
			}
		}
	}
	// to be finished~
}
int main(){
	while(input())solve();
}

/*
4 4 1 1 4 4 
 10  10  10 
9  0  0  10 
  0  0  0 
9  0  0  10 
  9  0  0 
0  9  0  10 
  0  9  9 

2 2 1 1 2 2 0 1 1 0 

0 0 0 0 0 0
*/