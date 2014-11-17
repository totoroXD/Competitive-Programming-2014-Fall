cD bta[1<<19], gar[1<<19], al[1<<19];
cD ba[1<<19], bb[1<<19], bg[1<<19];
double nu[1<<19];

void fft(cD xx[], int kh, int n, cD yy[], bool inv){
  if(n==1){
    yy[0] = xx[0];
    return ;
  }
  int h = (n+1)/2;
  cD *xeven = new cD[h], *xodd = new cD[h];
  cD *yeven = new cD[h], *yodd = new cD[h];
  cD w(1, 0), wn;
  for(int i=0;i<h;i++){
    xeven[i] = xx[i+i];
    if(i+i+1<n)xodd[i] = xx[i+i+1];
  }
  fft(xeven,kh/2, h, yeven, inv);
  fft(xodd, kh/2, n-h, yodd, inv);
  if(inv) wn = cD(cos(-2*pi/kh), sin(-2*pi/kh));
  else wn = cD(cos(2*pi/kh), sin(2*pi/kh));
  for(int i=0;i<h;i++){
    yy[i] = yeven[i] + w*yodd[i];
    if(i+h < n)yy[i + h] = yeven[i] - w*yodd[i];
    w *= wn;
  }
  delete xeven; delete xodd;
  delete yeven; delete yodd;
}

int main()
{
  int n, nn, nn2;
  int i;
  double t;

  while(scanf("%d",&n)!=EOF){
    for(nn2=1;nn2<n+n;nn2*=2);
    nn = n+n;
    assert(nn <= (1<<18));

    for(i=0;i<nn;i++){
      al[i] = bta[i] = gar[i] = cD(0,0);
      ba[i] = bb[i] = bg[i] = cD(0,0);
    }

    for(i=0;i<n;i++){
      scanf("%lf", &t);
      bta[i] = cD(t,0);
    }

    for(i=0;i<n;i++){
      scanf("%lf", &t);
      gar[i] = cD(t, 0);
    }
    fft(bta, nn2, nn, bb, false);
    fft(gar, nn2, nn, bg, false);
    for(i=0;i<nn;i++)
      ba[i] = bg[i]/bb[i];

    fft(ba, nn2, nn, al, true);
    printf("al:");
    for(i=0;i<nn;i++)
      printf("%.4lf ", double((al[i]).real())/double(nn));
    puts("");
    for(i=0;i<n;i++)
      printf("%.4lf\n", double((al[i]+al[i+nn/2]).real())/double(nn));

  }
  return 0;
}
