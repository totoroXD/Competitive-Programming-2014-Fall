struct pa{
    int val, i, od, sum;
    pa(int _val=0, int _i=0, int _od=0){
        val = _val;
        i = _i;
        od = _od;
    }
};
bool cmp1(const pa &a, const pa &b){
    if(a.val!=b.val)return a.val<b.val;
    return a.i<b.i;
}
bool cmp2(const pa &a, const pa &b){
    return a.i<b.i;
}
struct segment_tree{
    struct segment{
        vector<pa> num;
    }sg[401111];
    int n;
    void init(int _n){
        n=_n;
    }
    void init(vector<int> &ar){
        n = ar.size();
        vector<pa> a(n);
        for(int i=0; i<401111; i++)sg[i].num.clear();
        for(int i=0; i<n; i++){
            a[i].val=ar[i];
            a[i].i=i;
        }
        sort(a.begin(), a.end(), cmp1);
        for(int i=0; i<n; i++)a[i].od=i;
        sort(a.begin(), a.end(), cmp2);
        sg[0].num = a;
        build(0,n,0);
    }
    void build(int l, int r, int k){
        static bool clr[401111];
        vector<pa> &num = sg[k].num;
        if(num.size()<=1)return;
        for(int i=0; i<num.size(); i++){
            if(num[i].od < (l+r)/2){
                sg[k*2+1].num.push_back(num[i]);
                clr[num[i].i]=1;
            }else{
                sg[k*2+2].num.push_back(num[i]);
                clr[num[i].i]=0;
            }
        }
        for(int i=0; i<num.size(); i++){
            num[i].sum = clr[num[i].i];
            if(i)num[i].sum+= num[i-1].sum;
        }
        build(l, (l+r)/2, k*2+1);
        build((l+r)/2, r, k*2+2);
    }
    int kth_element(int l, int r, int k){
        return kth_element(l,r,k,0);
    }
    int kth_element(int l, int r, int k, int sid){
        vector<pa> &num = sg[sid].num;
        if(num.size()<=1)return num[0].val;
        int n_left = num[r-1].sum - (l?num[l-1].sum:0);
        int mid=num.size()/2, n=num.size();
        if(k<n_left){
            int nl = (l?num[l-1].sum:0) ;
            int nr = num[r-1].sum ;
            return kth_element(nl, nr, k, sid*2+1);
        }else{
            int nl = (l-(l?num[l-1].sum:0)) ;
            int nr = (r-num[r-1].sum) ;
            return kth_element(nl, nr, k-n_left, sid*2+2);
        }
    }
}st;
