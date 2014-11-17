#define MAX_L 1111
int v[MAX_L], prev[MAX_L], sa[MAX_L], rsa[MAX_L], lcpa[MAX_L];
bool cmp(int a, int b){ return v[a]<v[b]; }
void cptSAandLCPA(char s[], int sa[], int lcpa[]){
    int len=strlen(s);
    //SA
    for(int i=0; i<len; i++)sa[i]=i;
    for(int p=0,sh=1; sh<len; p++,sh*=2){
        if(p==0)
            for(int i=0; i<len; i++)v[i]=s[i];
        else{
            for(int i=0; i<len; i++){
                v[i]= prev[i]*MAX_L;
                if(i+sh/2<len)v[i]+= prev[i+sh/2];
            }
        }
        sort(sa, sa+len, cmp);
        for(int i=0; i<len; i++){
            prev[sa[i]] = i;
            if(0<i) if(v[sa[i-1]]==v[sa[i]]) prev[sa[i]] = prev[sa[i-1]];
        }
    }
    for(int i=0; i<len; i++)rsa[sa[i]]=i;
    // LCPA
    for(int lcp=0,i=0; i<len;i++){
        if(rsa[i]==0)
            lcpa[rsa[i]]=0;
        else{
            if(lcp>0)lcp--;
            int j=sa[rsa[i]-1];// behind i in SA order
            while(s[i+lcp] == s[j+lcp])lcp++;
            lcpa[rsa[i]]=lcp;
        }
    }
}
