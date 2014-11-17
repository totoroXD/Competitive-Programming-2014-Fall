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
