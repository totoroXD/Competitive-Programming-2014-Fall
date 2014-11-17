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