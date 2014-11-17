#define squ(x) ((x)*(x))
int ind;
struct pt{
  int id;
  ll xx[2];
  pt(){}
  pt(int _id, ll a, ll b){ id = _id; xx[0] = a; xx[1] = b; }
  pt(ll a, ll b){ xx[0] = a; xx[1] = b; }
  bool operator < (const pt &a) const { return xx[ind] < a.xx[ind]; }
  bool operator ==(const pt &a) const{ return xx[0]==a.xx[0] && xx[1] == a.xx[1];}
  pt operator -(const pt &a){ return pt(xx[0]-a.xx[0], xx[1]-a.xx[1]);}
};
ll dis(const pt &a, const pt &b){ return squ(a.xx[0]-b.xx[0])+squ(a.xx[1]-b.xx[1]); }
ll crs(const pt &a, const pt &b){ return a.xx[0]*b.xx[1] - a.xx[1]*b.xx[0]; }
ll dot(const pt &a, const pt &b){ return a.xx[0]*b.xx[0] + a.xx[1]*b.xx[1]; }

pt arr[MAXN];
pt point[MAXN];
priority_queue< pL > ansQ;

struct kd_tree{
  int child[MAXN << 3];
  pt node[MAXN << 3];
  void build(int LB, int RB, int rt=1, int dep=0){
    if(LB > RB) return ;
    child[rt] = RB - LB;
    child[rt<<1] = child[(rt<<1)|1] = -1;
    ind = dep%2;
    int mid = (LB+RB)>>1;
    nth_element(point+LB, point+mid, point+RB+1);
    node[rt] = point[mid];
    build(LB, mid-1, rt<<1, dep+1);
    build(mid+1, RB, (rt<<1)|1, dep+1);
  }
  void Q(const pt &p, int rt=1, int dep=0){
    if(-1 == child[rt]) return ;
    pL now = mp(dis(p, node[rt]), ll(node[rt].id));
    int dim = dep%2, lc = (rt<<1), rc = (rt<<1)|1;
    bool flag = false;
    if(p.xx[dim] >= node[rt].xx[dim]) swap(lc, rc);
    if(child[lc] != -1) Q(p, lc, dep+1);
    if(sz(ansQ) < 2){
      ansQ.push(now);
      flag = true;
    }else{
      if(now < ansQ.top())
          ansQ.pop(), ansQ.push(now);
      if(squ(p.xx[dim]-node[rt].xx[dim]) < ansQ.top().x)
        flag = true;
    }
    if(flag && child[rc]!=-1) Q(p, rc, dep+1);
  }
};

kd_tree tree;
vector< pt > rpt;
int coor(const pt &a){
  if(a.xx[1] >=0 ){
    if(a.xx[0] >=0 )return 0;
    else return 1;
  }else{
    if(a.xx[0] >=0 )return 3;
    else return 2;
  }
}
bool test(int now){
  pt o, p, q;
  int coorP, coorQ;
  int cnt = 0;
  o = arr[now];
  for(int i=0;i<sz(rpt);i++){
    if(rpt[i] == o)
      return true;

    p = rpt[i] - o;
    if(i+1==sz(rpt)) q = rpt[0] - o;
    else q = rpt[i+1] - o;
    coorP = coor(p);
    coorQ = coor(q);
    if(crs(p, q)==0 && dot(p, q) < 0) return true;

    if(coorP==coorQ)continue;
    if((coorP+1)%4==coorQ){
      cnt++;
      continue;
    }
    if((coorQ+1)%4==coorP){
      cnt--;
      continue;
    }
    if(crs(p, q) > 0) cnt += 2;
    else cnt -= 2;
  }
  return (cnt !=0);
}

int main(){
  int tn;
  int n;
  ll a, b;
  int R, qn;
  int bk;
  int i;
  int Sz;
  int ans1, ans2;

  cin >> tn;
  for(int z=1;z<=tn;z++){
    printf("Case #%d:\n", z);
    cin >> n;
    for(i=0;i<n;i++){
      cin >> a >> b;
      arr[i] = pt(i+1, a, b);
    }
    cin >> R;
    for(int r=1;r<=R;r++){
      printf("Region %d\n", r);
      cin >> bk;
      rpt.clear();
      for(i=0;i<bk;i++){
        cin >> a >> b;
        rpt.pb(pt(a,b));
      }
      reverse(all(rpt));

      for(i=Sz=0;i<n;i++)
        if(test(i))
          point[Sz++] = arr[i];
        
      tree.build(0, Sz-1);
      cin >> qn;
      while(qn--){
        cin >> a >> b;
        while(sz(ansQ))ansQ.pop();
        tree.Q(pt(a,b));
        ans2 = ansQ.top().y; ansQ.pop();
        ans1 = ansQ.top().y; ansQ.pop();
        printf("%d %d\n", ans1, ans2);
      }
    }
  }
  return 0;
}
