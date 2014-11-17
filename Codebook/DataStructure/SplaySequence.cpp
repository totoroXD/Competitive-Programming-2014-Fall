struct SplaySequence{
  static const int MAXN = 201111;
  struct Node{
    Node *ch[2], *par;
    int s, val, flip;
    Node(){ s = 0; val = -1; }
    int cmp(int k)const{
      int d = k-ch[0]->s;
      if(d==1)return -1;
      else if(d<=0)return 0;
      else if(d>1)return 1;
    }
    void maintain(){ s = ch[0]->s+ch[1]->s+1; }
    void pushDown(){
      if(flip){
        flip=0;
        swap(ch[0], ch[1]);
        ch[0]->reverse();
        ch[1]->reverse();
      }
    }
    Node * & parentPointer(){
      if(par->ch[0]==this)return par->ch[0];
      else return par->ch[1];
    }
    void reverse(){flip^=1;}
  };
  int n;
  Node seq[MAXN], *null, *root;
  SplaySequence(){ null = new Node(); }
  Node* newNode(int val){
    Node &o = seq[n];
    o.s = 1;
    o.flip = 0;
    o.val = val;
    o.ch[0] = null;
    o.ch[1] = null;
    return &seq[n++];
  }
  void rotateR(Node * &o, int d){
    Node *a = o->ch[d^1];
    a->ch[d]->par = o;
    o->ch[d^1] = a->ch[d];
    a->ch[d] = o;
    a->par = o->par;
    o->par = a;
    o->maintain();
    a->maintain();
    o = a;
  }
  void rotate(Node *o, int d){
    o->pushDown();
    if(o==root){
        rotateR(root, d);
    }else{
      int d2 = (o==o->par->ch[1]);
      rotateR(o->par->ch[d2], d);
    }
  }
  void splay (Node *o, Node* &rt){
    o->pushDown();
    while (o != rt){
      o->par->pushDown();
      int d = (o == o->par->ch[1]);
      if (o->par == rt){
        rotate(o->par, d ^ 1);
      }
      else{
        o->par->par->pushDown();
        if (o->par == o->par->par->ch[d]){
          rotate(o->par->par, d ^ 1);
          rotate(o->par, d ^ 1);
        }else{
          rotate(o->par, d ^ 1);
          rotate(o->par, d);
        }
      }
    }
  }
  void splay(Node *o){
    splay(o, root);
  }
  Node * splayK(Node * rt, int k){
    Node *v = rt;
    Node * &rtR = (rt==root?root:rt->parentPointer());
    while(1){
      v->pushDown();
      int d = v->cmp(k);
      if(d==-1){
        splay(v, rtR);
        return v;
      }
      if(d==1)k-= v->ch[0]->s+1;
      v = v->ch[d];
    }
  }
  Node * merge(Node * left, Node * right){// root doesn't make sense
    if(left==null)return right;
    root = left;
    left = splayK(left, left->s);
    left->ch[1] = right;
    right->par = left;
    left->maintain();
    root = null;
    return left;
  }
  void split(int k, Node * &left, Node * &right){// split root
    Node *o = splayK(root, k);

    root = o->par = null;
    left = o;
    right = o->ch[1];
    o->ch[1] = null;
    left->maintain();
  }
  void push_back(Node *o){
    Node *v = root;
    if (root == null){
      root = o;
      o->par = null;
    }
    else{
      while (v->ch[1] != null) v = v->ch[1];
      v->ch[1] = o;
      o->par = v;
      splay(o);
    }
  }
  Node* build(int sz) {
    if(!sz) return null;
    Node* L = build(sz/2);
    Node* o = newNode(n);
    Node* R = build(sz - sz/2 - 1);
    o->ch[0] = L;
    o->ch[1] = R;
    L->par = o;
    R->par = o;
    o->maintain();
    return o;
  }
  void clear(){
    n=0;
    root = null;
  }
  int size(){
    return root->s;
  }
  void sweep (Node *x){
    stack<Node*> stk;
    Node *v = x;
    while(v!=root){
      stk.push(v);
      v = v->par;
    }
    stk.push(root);
    while(!stk.empty()){
      stk.top()->pushDown();
      stk.pop();
    }
  }
  void init(int sz){
    n=0;
    root = build(sz);
  }
};
/*
   10
   2 3 4 5 6 7 8 9 10 1

*/
