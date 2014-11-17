LL powMod(LL a, LL p, LL m){
	if(p==0)return 1%m;
	LL res=powMod(a, p/2, m);
	res = res*res%m;
	if(p%2==1)
		res = res*a%m;
	return res;
}
LL gcd(LL a, LL b){
	if(b==0)return a;
	return gcd(b, a%b);
}
// |x|+|y| min
LL extGcd(LL a, LL b, LL &x, LL &y){
	LL res;
	if(b==0){
		x=1;
		y=0;
		res = a;
	}else{
		LL res = extGcd(b, a%b, y, x);
		y-= x*(a/b);
	}
	return res;
}
LL invMod(LL a, LL m){
	LL x, y;
	extGcd(a, m, x, y);
	return x;
}
int eulerPhi(int n){
	int m = sqrt(n+0.5);
	int res=n;
	for(int i=2; i<=m; i++){
		if(n%i==0){
			res = res*(i-1)/i;
			while(res%i==0)res/=i;
		}
	}
	if(n>1) res = res*(n-1)/n;
}
vector<int> phiTable(int n){
	vector<int> phi(n, 0);
	phi[1] = 1;
	for(int i=2; i<=n; i++){
		for(int j=i; j<=n; j++){
			if(!phi[j])phi[j] = j;
			phi[j] = phi[j]*(i-1)/i;
		}
	}
	return phi;
}
// a[i]*x = b[i](mod m[i])
// return pair(x, M)
typedef vector<int> vI;
pair<int,int> linearCongruence(const vI &a, const vI &b, const vI &m){
	// initially x=0(mod 1)
	int x=0, M=1;
	for(int i=0; i<a.size(); i++){
		int a_ = a[i]*M, b_=b[i]-a[i]*x, d=gcd(m[i],a_);
		if(b_%d!=0)return make_pair(0,-1);
		int t = b_ / d * invMod(a_/d, m[i]/d) %(m[i]/d);
		x = x+M*t;
		M*= m[i]/d;
	}
	return make_pair(x%M, M);
}
// x = a[i] (mod m[i]), a coprime each other
LL china(const vector<int> &a, const vector<int> &m){
	LL M=1, d, res=0, x, y;
	for(int i=0; i<m.size(); i++)M*= m[i];
	for(int i=0; i<a.size(); i++){
		LL w = M/m[i];
		extGcd(m[i], w, x, y);
		res = (res + w*y*a[i])%M;
	}
	return res;
}
// a^x = b (mod p), n is prime
int logMod(int a, int b, int p){
	int m=sqrt(p+0.5), v, e=1;
	v = invMod(powMod(a, m, p), p);
	map<int,int> mp;
	mp[1] = 0;
	for(int i=1; i<m ;i++){
		e = (e*a)%p;
		if(!mp.count(e))mp[e] = i;
	}
	for(int i=0; i<m; i++){
		if(mp.count(b))return i*m + mp[b];
		b = (b*v)%p;
	}
	return -1;
}
// n! (mod p) n! = a*p^e (p is prime)
// based on Wilson's theorem
const MAX_P = 101;
int fact[MAX_P];
int factMod(int n, int p, int &e){
	e = 0;
	if(n==0)return 1;
	int res = factMod(n/p, p, e);
	e+= n/p;
	if(n / p  %2 != 0)return res*(p-fact[n%p])%p;
	return res*fact[n%p]%p;
}
int combMod(int n, int m, int p){
	if(n<0 || k<0 || n<k)return 0;
	int e1, e2, e3;
	int a1 = factMod(n, p, e1), a2 = factMod(k, p, e2), a3 = factMod(n-k, p, e3);
	if(e1 > e2+e2)return 0;
	return a1* invMod(a2*a3%p, p) %p;
}
// Miller-Rabin
// fake prime test always correct when n<2^32
// if n < 4,759,123,141, it is enough to test a = 2, 7, and 61;
// if n < 3,825,123,056,546,413,051, it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, and 23.

bool MRTest(LL x, LL n){
	LL y = n-1;
	while(y%2==0)y/=2;
	LL z = powMod(x, y, n);
	if(z==1)return 1;
	while(y<n-1 && z!=1 && z!=n-1)
		z=(z*z)%n, y*=2;
	return z==n-1;
}
bool isPrime(LL n){
	if(n==2 || n==7 || n==61)return 1;
	if(n==1 || n%2==0) return 0;
	return MRTest(2, n) && MRTest(7, n) && MRTest(61, n);
}
//integration
double simpson(double (*f)(double), double a, double b){
	return (f(a)+4*f((a+b)/2)+f(b))*(b-a)/6;
}
double adaptiveSimpson(double (*f)(double), double a, double b, double eps){
	double c=(a+b)/2;
	double A=simpson(f,a,b), L=simpson(f,a,c), R=simpson(f,c,b);
	if(fabs(L+R-A)<=15*eps) return L+R+(L+R-A)/15;
	return adaptiveSimpson(f,a,c,eps/2) + adaptiveSimpson(f,c,b,eps/2);
}
//moebius function
map<int,int> moebius(int n){
	map<int,int> res;
	vector<int> primes;
	for(int i=2; i*i<=n; i++){
		if(n%i==0){
			primes.push_back(i);
			while(n%i==0)n/=i;
		}
	}
	if(n!=1) primes.push_back(n);
	int m = primes.size();
	for(int i=0; i<(1<<m); i++){
		int mu=1, d=1;
		for(int j=0; j<m; j++){
			if(i>>j &1){
				mu*= -1;
				d*= primes[j];
			}
		}
		res[d] = mu;
	}
	return res;
}
