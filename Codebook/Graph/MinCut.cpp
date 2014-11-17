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