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