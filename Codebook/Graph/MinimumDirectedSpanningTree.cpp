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
