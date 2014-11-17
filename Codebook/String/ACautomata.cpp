const int SIGMA_SIZE = 26;
const int MAXNODE = 11000;
const int MAXS = 150 + 10;

struct AhoCorasickAutomata {
  int ch[MAXNODE][SIGMA_SIZE];
  int fail[MAXNODE];
  int val[MAXNODE];  
  int last[MAXNODE]; 
  int cnt[MAXS];
  int sz;

  void init() {
    sz = 1;
    memset(ch[0], 0, sizeof(ch[0]));
    memset(cnt, 0, sizeof(cnt));
  }

  int idx(char c) {return c-'a';}

  void insert(char *s, int v) {
    int u = 0, n = strlen(s);
    for(int i = 0; i < n; i++) {
      int c = idx(s[i]);
      if(!ch[u][c]) {
        memset(ch[sz], 0, sizeof(ch[sz]));
        val[sz] = 0;
        ch[u][c] = sz++;
      }
      u = ch[u][c];
    }
    val[u] = v;
  }

  void print(int j) {
    if(j) {
      printf("%d\n",val[j]);
      print(last[j]);
    }
  }

  int find(char* T) {
    int n = strlen(T);
    int j = 0;
    for(int i = 0; i < n; i++) {
      int c = idx(T[i]);
      while(j && !ch[j][c]) j = fail[j]; 
      j = ch[j][c];
      if(val[j]) print(j);
      else if(last[j]) print(last[j]); 
    }
  }

  void getFail() {
    queue<int> q;
    fail[0] = 0;
    for(int c = 0; c < SIGMA_SIZE; c++) {
      int u = ch[0][c];
      if(u) { fail[u] = 0; q.push(u); last[u] = 0; }
    }
    while(!q.empty()) {
      int r = q.front(); q.pop();
      for(int c = 0; c < SIGMA_SIZE; c++) {
        int u = ch[r][c];
        if(!u) continue;
        q.push(u);
        int v = fail[r];
        while(v && !ch[v][c]) v = fail[v];
        fail[u] = ch[v][c];
        last[u] = val[fail[u]] ? fail[u] : last[fail[u]];
      }
    }
  }
};

