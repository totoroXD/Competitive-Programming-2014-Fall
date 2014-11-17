int mp(int tpl[], int n, int pat[], int m){
	if(n<m)return 0;
	int res=0;
	fail[0]=-1;
	for(int i=1,j=-1; i<m; i++){
		while(0<=j && pat[j+1]!=pat[i]) j=fail[j];
		if(pat[i]==pat[j+1])j++;
		fail[i] = j;
	}
	for(int i=0,j=-1; i<n; i++){
		while(0<=j && tpl[i]!=pat[j+1]) j=fail[j];
		if(pat[j+1]==tpl[i])j++;
		if(m-1<=j){
			res++;
			j=fail[j];
		}
	}
	return res;
}