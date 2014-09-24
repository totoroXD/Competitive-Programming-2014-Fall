#include <cstdio>

const int MAXLEN=101111, INF=1011110000;

struct SegTree{
	int n, seg[MAXLEN*4];// why *4?
	void init(int _n){
		n=_n;
		memset(seg,0,sizeof(seg));
	}
	void update(int i, int x){
		update(i, x, 0, n, 0);
	}
	int qmax(int l, int r){// [l, r-1]
		return qmax(l, r, 0, n, 0);
	}
	void update(int i, int x, int sl, int sr, int sid){
		if(sr-sl==1){
			seg[sid] = x;
			return;
		}
		if(sr<= i || i<sl)return;
		update(i, x, sl, (sl+sr)/2, sid*2+1);
		update(i, x, (sl+sr)/2, sr, sid*2+2);
	}
	int qmax(int l, int r, int sl, int sr, int sid){
		if(l<=sl && r<=sr)return seg[sid];
		if(sr<=l || r<=sl)return -INF;
		return max(qmax(l, r, sl, (sl+sr)/2, sid*2+1),
			qmax(l, r, sl, (sl+sr)/2, sid*2+1)
			);
	}
};
int main(){

}