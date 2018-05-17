#include<bits/stdc++.h>
using namespace std;
#define pi (2*acos(0.0))


class Point{
public:
    double x,y,z;

    void normalize(){
        double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        x/=r;
        y/=r;
        z/=r;
    }
};

class Matrix{
public:
    vector<vector<double> > mat;

    Matrix(int row,int col){
        for(int i=0;i<row;i++){
            vector<double> v;
            for(int j=0;j<col;j++)v.push_back(0);
            mat.push_back(v);
        }
    }

    void IdentityMatrix(){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                if(i==j) mat[i][j]=1;
            }
        }
    }

    void ColoumnMatrix(Point p){
        mat.clear();
        for(int i=0;i<4;i++){
            vector<double> v;
            v.push_back(1);
            mat.push_back(v);
        }
        mat[0][0] = p.x;
        mat[1][0] = p.y;
        mat[2][0] = p.z;
    }

    void showMatrix(int row,int col){
        for(int i=0;i<row;i++){
            for(int j=0;j<col;j++) cout << mat[i][j] << " ";
        cout << endl;
        }
    }

    void TranslationMatrix(int x,int y,int z){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                if(i==j)mat[i][j] = 1;
            }
        }
        mat[0][3] = x;
        mat[1][3] = y;
        mat[2][3] = z;
        //showMatrix(4,4);
    }

    void ScaleMatrix(int x,int y,int z){
        mat[0][0] = x;
        mat[1][1] = y;
        mat[2][2] = z;
        mat[3][3] = 1;
        //showMatrix(4,4);
    }

    void RotateMatrix(Point c1,Point c2,Point c3){
        //showMatrix(4,4);
        mat[0][0] = c1.x;
        mat[1][0] = c1.y;
        mat[2][0] = c1.z;
        mat[0][1] = c2.x;
        mat[1][1] = c2.y;
        mat[2][1] = c2.z;
        mat[0][2] = c3.x;
        mat[1][2] = c3.y;
        mat[2][2] = c3.z;
        mat[3][3] = 1;
        //showMatrix(4,4);
        //cout << endl;
    }

    void RMatrix(Point r,Point u,Point l){
        mat[0][0] = r.x;
        mat[0][1] = r.y;
        mat[0][2] = r.z;

        mat[1][0] = u.x;
        mat[1][1] = u.y;
        mat[1][2] = u.z;

        mat[2][0] = -l.x;
        mat[2][1] = -l.y;
        mat[2][2] = -l.z;

        mat[3][3] = 1;

    }

    void ProjectionMatrix(double far,double near,double r,double t){
        mat[0][0] = near/r;
        mat[1][1] = near/t;
        mat[2][2] = -(far+near)/(far-near);
        mat[2][3] = -(2*far*near)/(far-near);
        mat[3][2] = -1;
    }
};

Matrix multiplyMatrices(Matrix a, Matrix b) {

    int row1 = a.mat.size();
    int col1 = a.mat[0].size();
    int row2 = b.mat.size();
    int col2 = b.mat[0].size();

    Matrix product(row1, col2);

    for(int i=0; i<row1; i++) {
        for(int j=0; j<col2; j++) {
            for(int k=0; k<col1; k++) {
                product.mat[i][j] += a.mat[i][k] * b.mat[k][j];
            }
        }
    }

    return product;

}

Point transformPoint(Matrix m,Point p){
    Matrix P(4,1);
    P.ColoumnMatrix(p);

    Matrix multiple = multiplyMatrices(m,P);

    Point p1;

    p1.x = multiple.mat[0][0]/multiple.mat[3][0];
    p1.y = multiple.mat[1][0]/multiple.mat[3][0];
    p1.z = multiple.mat[2][0]/multiple.mat[3][0];

    return p1;
}

Point R(Point u, Point a, double angle) {

    double cosine = cos(angle * pi/180);
    double sine = sin(angle * pi/180);
    double m = (1 - cosine) * (a.x * u.x + a.y * u.y + a.z * u.z);

    Point p;

    /*p.x = cosine * u.x + m * a.x; + sine*(a.y * u.z - a.z * u.y);
    p.y = cosine * u.y + m * a.y + sine * (a.z * u.x - a.x * u.z);
    p.z = cosine * u.z + m * a.z + sine * (a.x * u.y - a.y * u.x);*/

    p.x = cosine*u.x;
    p.y = cosine*u.y;
    p.z = cosine*u.z;

    p.x += m * a.x;
    p.y += m * a.y;
    p.z += m * a.z;

    p.x += sine * (a.y * u.z - a.z * u.y);
    p.y += sine * (a.z * u.x - a.x * u.z);
    p.z += sine * (a.x * u.y - a.y * u.x);

    //cout << "Rotate Col:" << endl;
    //cout << p.x << " " << p.y << p.z << endl;
    //cout << endl;
    return p;

}


int main()
{
    ifstream input;
    ofstream output,output2,output3;
    output.open("stage1.txt", std::ios_base::out);
    output2.open("stage2.txt", std::ios_base::out);
    output3.open("stage3.txt", std::ios_base::out);
    output.precision(7);
    output2.precision(7);
    output3.precision(7);

    Point look,eye,up;
    Point l,r,u;
    double fovY,aspectRatio, near,far ;

    stack<Matrix> S;
    stack<stack <Matrix> > buffer;

    input.open("scene3.txt");

    input >> eye.x >> eye.y >> eye.z;
    input >> look.x >> look.y >> look.z;
    input >> up.x >> up.y >> up.z;
    input >> fovY >> aspectRatio >> near >> far;

    l.x = look.x - eye.x;
    l.y = look.y - eye.y;
    l.z = look.z - eye.z;
    l.normalize();

    r.x = l.y*up.z - l.z*up.y;
    r.y = l.z*up.x - l.x*up.z;
    r.z = l.x*up.y - l.y*up.x;
    r.normalize();

    u.x = r.y * l.z - r.z * l.y;
    u.y = r.z * l.x - r.x * l.z;
    u.z = r.x * l.y - r.y * l.x;

    Matrix T(4,4);
    T.TranslationMatrix(-eye.x,-eye.y,-eye.z);

    Matrix Rt(4,4);
    Rt.RMatrix(r,u,l);

    Matrix VT(4,4);
    VT = multiplyMatrices(Rt,T);

    double fovX = fovY*aspectRatio;
    double t = near*tan((fovY * pi) / 360);
    double ro = near*tan((fovX * pi) / 360);

    Matrix ProM(4,4);
    ProM.ProjectionMatrix(far,near,ro,t);
    //Pr.showMatrix(4,4);

    Matrix m = Matrix(4,4);
    m.IdentityMatrix();
    S.push(m);

    while(1){
        string command;

        input >> command;
        //setprecision(8);

        if(command=="triangle"){
            Point p[3];

            for(int i=0;i<3;i++){
                input >> p[i].x >> p[i].y >> p[i].z;
                Point p1 = transformPoint(S.top(),p[i]);
                Point p2 = transformPoint(VT,p1);
                Point p3 = transformPoint(ProM,p2);
                output << fixed << p1.x << " " << p1.y << " " << p1.z << endl;
                output2 << fixed << p2.x << " " << p2.y << " " << p2.z << endl;
                output3 << fixed << p3.x << " " << p3.y << " " << p3.z << endl;
                //cout << fixed << p1.x << " " << p1.y << " " << p1.z << endl;
                //cout << fixed << p2.x << " " << p2.y << " " << p2.z << endl;
                //cout << fixed << p3.x << " " << p3.y << " " << p3.z << endl;
            }
            output << endl;
            output2 << endl;
            output3 << endl;
            //cout << endl;
        }
        else if(command=="translate"){
            double dx,dy,dz;
            input >> dx >> dy >> dz;
            Matrix m(4,4);
            m.TranslationMatrix(dx,dy,dz);
            S.push(multiplyMatrices(S.top(),m));
        }
        else if(command=="scale"){
            double sx,sy,sz;
            input >> sx >> sy >> sz;
            Matrix m(4,4);
            m.ScaleMatrix(sx,sy,sz);
            S.push(multiplyMatrices(S.top(),m));
        }
        else if(command=="rotate"){
            double angle;
            Point a;

            input >> angle >> a.x >> a.y >> a.z;
            double r = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

            /*a.x /= r;
            a.y /= r;
            a.z /= r;*/
            a.normalize();

            Point i,j,k;
            i.x = 1;
            i.y = 0;
            i.z = 0;
            j.x = 0;
            j.y = 1;
            j.z = 0;
            k.x = 0;
            k.y = 0;
            k.z = 1;

            Point c1 = R(i,a,angle);
            Point c2 = R(j,a,angle);
            Point c3 = R(k,a,angle);

            /*cout << "Rotation : \n";
            cout << fixed << c1.x << " " << c1.y << " " << c1.z << endl;
            cout << fixed << c2.x << " " << c2.y << " " << c2.z << endl;
            cout << fixed << c3.x << " " << c3.y << " " << c3.z << endl;
            cout << endl;*/

            Matrix R(4,4);
            R.RotateMatrix(c1,c2,c3);
            S.push(multiplyMatrices(S.top(),R));
        }

        else if(command=="push"){
            buffer.push(S);
        }
        else if(command=="pop"){
            S = buffer.top();
            buffer.pop();
        }

        else if(command=="end") return 0;
    }
}
