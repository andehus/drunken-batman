#include <iostream>
#include <math.h>
#include <stdlib.h>
//#include </uio/hume/student-u32/andehus/fys3150/compphys/programs/cppLibrary/lib.cpp>
#include </home/anders/Project1/Project2/lib.cpp>
#include <armadillo>
#include <time.h>
#include <fstream>
using namespace std;

class Matrix {
    ofstream f;
public:
    Matrix(){
    }

    /*Allocate memory for an nxn matrix
     *Return: The allocated matrix
     */
    double** alloc_matrix(double** temp, int n) {
        int i;
        temp = new double*[n];
        for(i = 0; i<n; i++) {
            temp[i] = new double[n];
        }
        return &(*temp);
    }

    /*Put all the elements to zero
     *Return: The initialized matrix
     */
    double** zero_matrix(double** temp, int n) {
        int i,j;
        for(i = 0; i <n; i++) {
            for(j = 0; j<n; j++) {
                temp[i][j] = 0.0;
            }
        }
        return &(*temp);
    }

    /*Creates an identity matrix*/
    double** identity_matrix(double** temp, int n) {
        int i,j;
        for(i = 0; i<n; i++) {
            for(j = 0; j<n; j++) {
                if(i == j) {
                    temp[i][j] = 1.0;
                }
                else temp[i][j] = 0.0;
            }
        }
        return &(*temp);
    }

    /*Sorting diagonal elements in a matrix from low to high
     *Parameter: nxn matrix, dimension n
     *Return: array of n double elements
     */
    double* sort(double* temp,int n) {
        double t;
        int i;
        int swapped = 1;
        while(swapped) {
            swapped = 0;
            for(i = 0;i<n-1; i++){
                if(temp[i]>temp[i+1]) {
                    t = temp[i];
                    temp[i] = temp[i+1];
                    temp[i+1] = t;
                    swapped = 1;
                }
            }
        }
        return temp;
    }

    /*Collect the values on the diagonal of an nxn matrix
     *Parameter: nxn double array, length n
     *Return: list of n doubles
     */
    double* get_diagonal(double** temp, int n){
        double* list = new double[n];
        for(int i = 0; i<n; i++)
            list[i] = temp[i][i];
        return list;
    }

    /*Finding the index in *array which contains the value t
     *Return: index which contains the value t
     */
    int get_index(double t,double*array,int n) {
        int i = 0;
        int index = 0;
        while(i<n && t!=array[i])
            i++;
        cout<<"Index: "<<i<<endl;
        return i;
    }

    /*Puts the row elements of A in column "index" into vector v
     *Returns the vector v
     */
    double* get_vec(int index, double** M, double* v,int n){
        for(int i = 0; i<n; i++) {
            v[i] = M[i][index];
        }
        return v;
    }

    /*Prints array to file
     *Return: void
     */
    void atf(double* vec,double* vec2, int n) {
        f.open("/uio/hume/student-u32/andehus/fys3150/project2_2/vector.dat");
        for(int i = 0; i<n; i++) {
            f<<vec[i]<<" "<<vec2[i]<<endl;
        }
        f.close();
    }

    /*Prints the elements in the array
     *Parameter: double array, length
     *Return: void
     */
    void print_array(double* temp,int n) {
        for(int i = 0; i<3; i++)
            cout<<temp[i]<<endl;
    }

    /*Prints the matrix
     *Parameter: nxn matrix, dim n
     */
    void print(double** temp, int n) {
        int i,j;
        for(i = 0; i<n; i++) {
            for(j = 0; j<n; j++) {
                cout<<temp[i][j]<<" ";
            }
            cout<<endl;
        }
    }

};

class Project2 {
    double *p, *v, *d, e,*e_array, *u, **A,**B, **I, h, epsilon,*eigenvalue,*diagonal;
    double p_min,p_max;
    double *vector1, *vector2, *vector3;
    int n,k,l,index;
    clock_t start,finish;
    Matrix *matrix;
public:
    Project2() {
        matrix = new Matrix();
        p_min = 0;
        p_max = 5;
        n = 196;
        h = (double)((p_max - p_min)/(n+1));
        epsilon =1.0e-9;
        p = new double[n];
        v = new double[n];
        d = new double[n];
        e = -1/(h*h);
        e_array = new double[n];
        u = new double[n];

        init();
        A = matrix->alloc_matrix(&(*A),n);
        construct_tridmatrix();
        I = matrix->alloc_matrix(&(*I),n);
        I = matrix->identity_matrix(&(*I),n);
        start = clock();
        jakobi(A,&(*I),n);
        //tqli(d,e_array,n,&(*A));
        diagonal = matrix->get_diagonal(&(*A),n);
        eigenvalue = matrix->get_diagonal(&(*A),n);
        //cout<<"nye egenverdier"<<endl;
        //for(int k = 0; k<n; k++)
        //    cout<<d[k]<<endl;
        //eigenvalue = d;

        eigenvalue = matrix->sort(eigenvalue,n);
        vector1 = new double[n];
        index = matrix->get_index(eigenvalue[0],diagonal,n);
        vector1 = matrix->get_vec(index,I,vector1,n);
        matrix->atf(vector1,p,n);
        //matrix->print(I,n);
        cout<<"Eigenvalues: "<<endl;
        matrix->print_array(eigenvalue,n);
        finish = clock();
        cout<<"Time: "<<(double)((finish-start)/CLOCKS_PER_SEC)<<"sec"<<endl;
        //print();
    }
    int print() {
        int i;
        cout<<"p";
        cout<<"v";
        cout<<"d";
        cout<<"u"<<endl;
        for(i = 0; i <n; i++) {
            cout<<p[i]<<" ";
            cout<<v[i]<<" ";
            cout<<d[i]<<" ";
            cout<<u[i]<<endl;
        }
        return 0;
    }
    int rotate(double** A,double **I,int k, int l, int n){
        int i;
        double s,c,t,tau,a_kk,a_ll,a_il,a_ik,i_ik,i_il;
        if(A[k][l] != 0.0) {
            tau = (A[l][l]-A[k][k])/(2*A[k][l]);
            if(tau>0) {
                t = 1.0/(tau+sqrt(1.0+tau*tau));
            }
            else {
                t = -1.0/(-tau+sqrt(1.0+tau*tau));
            }
            c = 1/sqrt(1+t*t);
            s = c*t;
        }
        else {
            c = 1.0;
            s = 0.0;
        }
        a_kk = A[k][k];
        a_ll = A[l][l];
        A[k][k] = c*c*a_kk-2.0*c*s*A[k][l]+s*s*a_ll;
        A[l][l] = s*s*a_kk+2.0*c*s*A[k][l]+c*c*a_ll;
        A[k][l] = 0.0;
        A[l][k] = 0.0;
        for(i = 0; i<n; i++) {
            if(i != k && i != l) {
                a_ik = A[i][k];
                a_il = A[i][l];
                A[i][k] = c*a_ik - s*a_il;
                A[k][i] = A[i][k];
                A[i][l] = c*a_il + s*a_ik;
                A[l][i] = A[i][l];
            }
            i_ik = I[i][k];
            i_il = I[i][l];
            I[i][k] = c*i_ik - s*i_il;
            I[i][l] = c*i_il + s*i_ik;
        }
        return 0;
    }

    double max_offdiagonal(double **A, int*k, int*l, int n) {
        int i,j;
        double max = 0.0;
        for(i = 0; i<n; i++) {
            for(j = i+1; j<n; j++) {
                if(fabs(A[i][j])>max) {
                    max = fabs(A[i][j]);
                    *l = i;
                    *k = j;
                }
            }
        }
        return max;
    }

    int jakobi(double**A, double**I, int n) {
        int k, l;
        int max_iter = n*n*n;
        int iter = 0;
        double max_off = max_offdiagonal(A,&k,&l,n);
        while(fabs(max_off)>epsilon && iter < max_iter) {
            max_off = max_offdiagonal(A,&k,&l,n);
            rotate(A,I,k,l,n);
            iter++;
        }
        cout<<"Number of iter"<<iter<<endl;
        return 0;
    }

    int construct_tridmatrix() {
        int i;
        A[0][0] = d[0];//d1
        A[0][1] = e;//e1
        for(i  = 1; i<n-1; i++) {
            A[i][i-1] = e;//elm(1,0)
            A[i][i] = d[i];//elm(1,1)
            A[i][i+1] = e;//elm(1,2)
        }
        A[i][i-1] = e;
        A[i][i] = d[i];
        return 0;
    }

    int init() {
        int i;
        double omega = 1.0;//0.01,0.5,1.0,5.0
        double beta_e2 = 1.44;

        for(i = 0; i<n; i++) {
            p[i] = p_min +(i+1)*h;//p(1)
            v[i] = p[i]*p[i];//v(1)
            //v[i] = omega*omega*p[i]*p[i]+1/p[i]+beta_e2/p[i];
            d[i] = (double)2/(h*h)+(double)v[i];//+v(1)
            u[i] = 0;
            e_array[i] = -1/(h*h);
        }
        return 0;
    }
};

int main()
{
    Project2 *p2 = new Project2();
    return 0;
}

