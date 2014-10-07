/*
 * totoroXD
 *
 */
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <queue>
#include <stack>
using namespace std;
typedef long long LL;
const int INF = 1011110000, MOD=1000000000;
const int dx[]={1,0,-1,0}, dy[]={0,1,0,-1};
const double EPS = 1e-6;
typedef long long LL;
int n, m, k, t;

char str[111];
int bb[44];

bool input(){
	memset(bb,0,sizeof(bb));
	if(scanf("%d",&n)==EOF)return 0;
	for(int i=0; i<n; i++){
		scanf("%s",str);
		for(int j=0; j<strlen(str); j++)
			bb[i] |= 1<<(str[j]-'A');
	}
	return 1;
}
void solve(){
	map<int, pair<int,int> > nBone;// nBone[ alphaSet ]: (max # of Bones, boneSet)
	for(int s=0; s<(1<<(n/2)); s++){// s: all subset of the first n/2 bones
		unsigned int k=0,cnt=0;
		for(int j=0; j<(n/2); j++){
			if(s&(1<<j)) k^= bb[j], cnt++;
		}
		pair<int,int> t(cnt,s);// (# of Bones, boneSet)
		nBone[k] = max(nBone[k], t);
	}
	for(map<int, pair<int,int> >::iterator it=nBone.begin(); it!=nBone.end(); it++){// DFS in BST
		int a = it->first;
		int nb = it->second.first;
		int s = it->second.second;
		printf("[%d] = (%d, %d)\n",a,nb,s);
	}
}
int main(){
	while(input()){
		solve();
	}
    return 0;
}

/*
6
ABD
EG
GE
ABE
AC
BCD

10
A
B
C
D
E
F
G
H
I
J

*/

