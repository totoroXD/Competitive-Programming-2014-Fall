void makeZ( int z[], char s[], int l ){
	z[0] = l;
	int L=0, R=0, i, x;
	for( i = 1 ; i < l ; i++ ){
		if( R < i || z[i-L] >= R-i+1 ){
			R < i ? x = i : x = R+1;
			while( x < l && s[x] == s[x-i] ) x++;
			z[i] = x-i;
			if( i < x ){ L = i; R = x-1; }
		}
		else z[i] = z[i-L];
	}
}
