#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include </home/anders/programs/cppLibrary/lib.cpp>
#define h 1.0/((double)n+1)

using namespace std;
using namespace arma;

class Matrix {
    fstream file;
    clock_t start, finish;
    int n, v0, v1;
    double *a,*b,*c,*f,*v, **matrix;
public:
    Matrix(int n, int v0, int v1){
        this->n = n;
        this->v0 = v0;
        this->v1 = v1;
        a = new double[n];
        b = new double[n];
        c = new double[n];
        f = new double[n];
        v = new double[n];
        matrix = new double*[n];
        for(int i = 0; i<n;i++){
            matrix[i] = new double[n];
        }
        init_vec(a,-1.0);
        init_vec(b,2.0);
        init_vec(c,-1.0);
        init_vec(f);
        init_matrix(a, b, c);
        //a[0]=0;
        //c[n-1]=0;
        /*
        print(a);
        print(b);
        print(c);
        print(f);
        print(v);
*/
        //print();
    }
public:
    int ludec(){
        //void ludcmp(double **a, int n, int *indx, double *d)
        start = clock();
        int *t = new int[n];
        double *s = new double[n];
        ludcmp(matrix, n, t, s);
        finish = clock();
        cout<<"Time (LU-dec): "<<((finish - start)/CLOCKS_PER_SEC)<<endl;
        return 0;
    }

private:
    int init_matrix(double *a, double *b, double *c) {
        int i;
        matrix[0][0] = b[0];
        matrix[0][1] = c[0];
        for(i = 1; i<n-1; i++) {
            matrix[i][i-1] = a[i];
            matrix[i][i] = b[i];
            matrix[i][i+1] = c[i];
        }
        matrix[i][i-1] = a[i];
        matrix[i][i] = b[i];
        return 0;
    }
public:
    int printmatrix() {
        cout << "Matrix: "<<endl;
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++){
                cout <<matrix[i][j];
            }
            cout <<endl;
        }
        return 0;
    }

public:
    int print(double *temp) {
        for(int i = 0; i<n; i++)
            cout<<temp[i]<<endl;
        return 0;
    }
public:
    int print_f(){
        for(int i = 0; i<n; i++){
            cout<<f[i]<<endl;
        }
        return 0;
    }

public:
    int print(){
        for(int i = 0; i<n; i++)
            cout<<v[i]<<endl;
        return 0;
    }
private:
    /*Metode for å initialisere f(x_i)*/
    int init_vec(double *temp){
        for(int i = 0; i<n; i++){
            temp[i] = h*h*100*exp(-10*i*h);
        }
        return 0;
    }

private:
    void init_vec(double *temp,double value) {
        for(int i = 0; i<n; i++) {
            temp[i] = value;
        }
    }

public:
    int solver(){
        start = clock();
        int i;
        double *temp;
        double btemp;
        temp = new double[n];
        //Forward
        btemp = b[0];
        cout<<"F(0) = "<<f[0]<<endl;
        v[0] = f[0]/btemp;
        cout<<"v(0) = "<<v[0]<<endl;

        for(i = 1; i<n; i++){
            temp[i] = c[i-1]/btemp;
            btemp = b[i]-a[i]*temp[i];
            v[i] = (f[i]-a[i]*v[i-1])/btemp;
        }
        //Backward
        for(i = n-2; i>0; i--) {
            v[i] -=temp[i+1]*v[i+1];
        }
        finish = clock();
        double diff = ((double)finish-(double)start)/1000000;
        double seconds = (double)diff/CLOCKS_PER_SEC;
        cout<<"time: "<<seconds<<endl;
        cout<<"Time: "<<((finish - start)/CLOCKS_PER_SEC)<<endl;

        return 0;
    }
public:
    int write_tofile(char *filename){
        file.open(filename, fstream::in |fstream::out|fstream::trunc);
        if(!file.is_open()){
            cout <<"error"<<endl;
        }
        double epsilon;
        double temp,max;
        max = 0;
        for(int i = 0; i<n; i++){
            file<<setiosflags(ios::showpoint | ios::uppercase);
            file<<setw(15)<<setprecision(8)<<i*h;
            file<<setw(15)<<setprecision(8)<<v[i];
            temp = 1-(1-exp(-10))*i*h-exp(-10*i*h);
            file<<setw(15)<<setprecision(8)<<temp<<endl;
            if(temp == 0)
               epsilon = -10;
            else
                epsilon = log10(abs((v[i]-temp)/temp));
            //Finner maks error
            if(max<epsilon)
                max=epsilon;
        }
        file.close();
        cout<<"Max: "<<max<<endl;
        return 0;
    }
};
int main()
{
    int n, v0,v1;
    //Dirichlet boundary conditions
    v0 = v1 = 0;
    cout <<"Insert dimension n: "<<endl;
    cin >>n;
    Matrix *m = new Matrix(n, v0, v1);
    m->solver();
    //m->print();
    m->write_tofile("pdata1.dat");
    //m->printmatrix();
    m->ludec();
    return 0;
}
