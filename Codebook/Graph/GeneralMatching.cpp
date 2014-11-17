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
