int n;
int G[MAXN][MAXN];
int dp[MAXN], best;

bool dfs(int id[],int top,int cnt)
{
    if(!top){
        if(best<cnt){
            best = cnt;
            return true;
        }
        return false;
    }
    int a[MAXN];
    for(int i=0;i<top;i++){
        if(cnt+top-i<=best) return false;
        if(cnt+dp[id[i]]<=best) return false;
        int k = 0;
        for(int j=i+1;j<top;j++)
            if(G[id[i]][id[j]])
                a[k++] = id[j];
        if(dfs(a,k,cnt+1)) return true;
    }
    return false;
}

int solve(void)
{
    int id[MAXN];
    best = 0;
    for(int i=n-1;i>=0;i--){
        int top = 0;
        for(int j=i+1;j<n;j++)
            if(G[i][j])
                id[top++] = j;
        dfs(id,top,1);
        dp[i] = best;
    }
    return best;
}