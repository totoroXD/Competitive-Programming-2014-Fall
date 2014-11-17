// POJ 1741
struct edge {int to, cost; };
vector<edge> G[MAX_N], dis;
int minval, fctd;
bool vst[MAX_N];


int find_w(int v, int par, int V)
{
	int w=1, val=0;
	for(int i=0; i<G[v].size(); i++){
		edge &e = G[v][i];
		if(!vst[e.to] && e.to!=par){
			int wson = find_w(e.to, v, V);
			w+= wson;
			val = max(wson, val);
		}
	}
	val = max(V-w, val);
	if(val<minval)
	{
		minval = val;
		fctd = v;
	}
	return w;
}
int find_ctd(int root, int V)
{
	minval = INF;
	find_w(root, root, V);
	return fctd;
}
void DFS(int v, int par, int d)		// find distance and size of a tree
{
	dis.push_back((edge){v, d});
	for(int i=0; i<G[v].size(); i++){
		edge &e = G[v][i];
		if(!vst[e.to] && e.to!=par)
			DFS(e.to, v, d+e.cost);
	}
}
bool cmp(edge a, edge b){ return a.cost<b.cost; }

int find_pair(int root, int k)
{
	int res=0;
	dis.clear();
	DFS(root, root, 0);
	sort(dis.begin(), dis.end(), cmp);
	for(int i=0, j=dis.size()-1; i<j; i++)
	{
		while(i<j && dis[j].cost+dis[i].cost>k) j--;
		res+= j-i;
	}
	return res;
}
int solve(int root, int k)
{
	dis.clear();
	DFS(root, root, 0);
	int res=0, V=dis.size();
	root = find_ctd(root, V);
	res+= find_pair(root, k);

	vst[root] = 1;
	for(int i=0; i<G[root].size(); i++)
	{
		edge &e = G[root][i];
		if(!vst[e.to])
		{
			res+= solve(e.to, k);
			if(k-2*e.cost>=0)
			res-= find_pair(e.to, k-2*e.cost);
		}
	}
	vst[root] = 0;
	return res;
}