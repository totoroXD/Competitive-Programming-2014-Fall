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
