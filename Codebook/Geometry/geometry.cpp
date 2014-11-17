const double EPS = 1e-8, PI = 3.14159265358979323846;
int dcmp(double x){
        if(x<-EPS)return -1;
        else if(EPS)return 1;
        else return 0;
}
bool between(double x, double a, double b){// closed
        if(b < a)swap(a,b);
        return 0<=dcmp(x-a) && 0<=dcmp(b-x);
}
struct Vector;
struct Line;
struct Circle;
typedef Vector Point;
typedef Line Seg;
struct Vector{
    double x, y;
    Vector(double _x=0, double _y=0){
        this->x = _x;
        this->y = _y;
    }
    Vector operator + (const Vector &that)const{
        return Vector(x+that.x, y+that.y);
    }
    Vector operator - (const Vector &that)const{
        return Vector(x-that.x, y-that.y);
    }
    Vector operator * (double that)const{
        return Vector(x*that, y*that);
    }
    Vector operator / (double that)const{
        return Vector(x/that, y/that);
    }
    double operator ^ (const Vector &that)const{// dot
        return x*that.x + y*that.y;
    }
    double operator * (const Vector &that)const{// cross
        return x*that.y - y*that.x;
    }
    bool operator <(const Vector &that)const{
        return dcmp(x-that.x)<0 || (dcmp(x-that.x)==0 && dcmp(y-that.y)<0);
    }
    bool operator == (const Vector &that)const{
        return dcmp(x-that.x)==0 && dcmp(y-that.y)==0;
    }
    double length()const{
        return sqrt(x*x+y*y);
    }
    double angle()const{//
        return atan2(y,x);
    }
    double interAngle(const Vector &that)const{
        return acos(((*this)^that)/this->length()/that.length());
    }
    double disToPoint(const Point &that)const{
        return (that-*this).length();
    }
    double disToLine(const Line &that)const;
    double disToSeg(const Seg &that)const;
    bool onLine(const Line &that)const;
    bool onSeg(const Seg &that)const;
    Vector projectionToLine(const Line & that)const;
    Vector rotate(double rad)const{
        return Vector(x*cos(rad)-y*sin(rad), x*sin(rad)+y*cos(rad));
    }
    Vector unit()const{
        return (*this)/this->length();
    }
    Vector normal()const{
        return Vector(-y,x).unit();
    }
};
struct Line{
    Point o;
    Vector v;
    Line(){}
    Line(Point _o,Vector _v){
        this->o = _o;
        this->v = _v;
    }
    Vector e()const{
        return o+v;
    }
    bool operator < (const Line &that)const{
        if( dcmp(v*that.v) == 0 )return dcmp(v*(that.o-o)) > 0;
        return dcmp(v*that.v) > 0;
    }
    bool operator == (const Line &that)const{
        return !(*this<that) && !(that<*this);
    }
    Point intersectionWithLine(const Line &that)const;
    vector<Point> intersectionWithCircle(const Circle &that)const;
};
struct Circle{
    Point c;
    double r;
    Circle(Point _c, double _r):c(_c), r(_r){}
    Point point(double ang)const{
        return Point(c.x + r*cos(ang), c.y+r*sin(ang));
    }
    bool operator == (const Circle &that)const{
        return c==that.c && dcmp(r-that.r)==0;
    }
    vector<Point> intersectionWithLine(const Line &that)const;
    vector<Point> intersectionWithCircle(const Circle &that)const;
};
double Point::disToLine(const Line &that)const{
    return fabs((that.v * (*this-that.o)) / that.v.length());
}
double Point::disToSeg(const Seg &that)const{
    double crossa = that.v ^ (*this-that.o);
    double crossb = that.v ^ (*this-that.e());
    if(dcmp(crossa*crossb) < 0)
        return this->disToLine(that);
    else
        return min(this->disToPoint(that.o), this->disToPoint(that.e()));
}
bool Point::onLine(const Line &that)const{
    return dcmp(this->disToLine(that))==0;
}
bool Point::onSeg(const Seg &that)const{
    return dcmp(this->disToSeg(that))==0;
}
Point Point::projectionToLine(const Line &that)const{
    const Point &o = that.o;
    const Vector &v = that.v;
    return o + v*( (v^(*this-o))/(v^v));
}
Point Line::intersectionWithLine(const Line &that)const{
    Vector u = o-that.o;
    double t = (that.v*u) / (v*that.v);
    return o+v*t;
}
vector<Point> Line::intersectionWithCircle(const Circle &that)const{
    vector<Point> res;
    Point c = that.c, p;
    double d = c.disToLine(*this), r = that.r;
    if(dcmp(r-d)>=0){
        p = c.projectionToLine(*this);
        double l = sqrt(r*r-d*d);
        if(dcmp(r-d)==0){
            res.push_back(p);
        }else{
            res.push_back(p+v.unit()*l);
            res.push_back(p-v.unit()*l);
        }
    }
    return res;
}
vector<Point> Circle::intersectionWithLine(const Line &that)const{
    return that.intersectionWithCircle(*this);
}
vector<Point> Circle::intersectionWithCircle(const Circle &that)const{
    vector<Point> res;
    double d = (c-that.c).length();
    if(dcmp(d)==0);
    else if(dcmp(r+that.r-d)<0);
    else if(dcmp(fabs(r-that.r)-d)>0);
    else{
        double ang = (that.c-c).angle();
        double da = acos((r*r+d*d-that.r*that.r)/(2*r*d));
        Point p1 = point(ang-da), p2 = point(ang+da);
        res.push_back(p1);
        if(!(p1==p2))
            res.push_back(p2);
    }
    return res;
}
double triangleArea(const Point &a, const Point &b, const Point &c){
    Vector va = a-c, vb = b-c;
    return abs(va*vb)/2;
}
vector<Point> halfplane(vector<Line> l)
{
    sort(l.begin(),l.end());
    l.resize(int(unique(l.begin(),l.end())-l.begin()));
    vector<Line> q;
    q.resize(int(l.size()));
    int L=0,R=0;
    for( int i = 0; i < int(l.size()); i++ ){
        while( R-L >= 2 ){
            Point p = q[R-1].intersectionWithLine(q[R-2]);
            if( dcmp((l[i].v)*(p-l[i].o)) <= 0 )R--;
        }
        while( R-L >= 2 ){
            Point p = q[L].intersectionWithLine(q[L+1]);
            if( dcmp((l[i].v)*(p-l[i].o)) <= 0 )L++;
        }
        q[R++] = l[i];
    }
    while( R-L >= 2 ){
        Point p = q[R-1].intersectionWithLine(q[R-2]);
        if( dcmp((q[L].v)*(p-q[L].o)) <= 0 )R--;
    }
    vector<Point> res;
    if( R-L < 3 )return res;
    for( int i = L; i < R-1; i++ )
        res.push_back(q[i].intersectionWithLine(q[i+1]));
    res.push_back(q[L].intersectionWithLine(q[R-1]));
    return res;
}
