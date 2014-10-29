#include <iostream>     // std::cout
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>
#include <queue>

using namespace std;
typedef long long LL;
int n, m, k;
LL dp[111][11111];

struct Pair{
	LL x, y;// (welcome value, time cost),
			// (prefix ..    , prefix..)
	Pair(LL _x=0, LL _y=0):x(_x), y(_y){}
}road[111][11111];

bool input(){
	scanf("%d%d%d",&n,&m,&k);
	if(n+m+k==0)return 0;
	for(int i=0; i<=n; i++)
		for(int j=0; j<m; j++)
			scanf("%lld",&road[i][j].x);
	for(int i=0; i<=n; i++)
		for(int j=0; j<m; j++)
			scanf("%lld",&road[i][j].y);
	return 1;
}

void solve(){
	deque<Pair> q;
	memset(dp,0,sizeof(dp));
	for(int i=0; i<=n; i++){
		q.clear();
		for(LL j=0,xsum=0,ysum=0; j<=m; j++){// from left-hand side
			xsum+= (j==0? 0: road[i][j-1].x);
			ysum+= (j==0? 0: road[i][j-1].y);
			LL xpre = (i==0?0:dp[i-1][j])-xsum;
			while(!q.empty() && q.back().x<=xpre) q.pop_back();
			q.push_back(Pair(xpre,ysum));
			while(ysum-q.front().y>k)q.pop_front();
			dp[i][j] = max(dp[i][j], xsum+q.front().x);
		}/*
			// from right-hand side
		*/

	}
	LL res = 0;
	for(int i=0; i<=m; i++)
		res = max(res, dp[n][i]);
	printf("%lld\n",res);
}
int main () {
	while(input())solve();
	return 0;
}
/*

2 3 2
7 8 1
4 5 6
1 2 3
1 1 1
1 1 1
1 1 1
0 0 0

*/
