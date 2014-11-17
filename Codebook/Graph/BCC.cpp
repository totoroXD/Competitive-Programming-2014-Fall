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