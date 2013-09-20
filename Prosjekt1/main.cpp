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
/*Klasse med konstruktør og metoder for oppgave e)
 *Kommenter vekk opprettelsen av objektet i main,
 *hvis ikke oppgave e) skal regnes ut. Multiplikasjonen
 *blir kalt fra konstruktøren.
 */
class Multiplication {
    double **mata, **matb, **matc;
    int n;
    clock_t start, finish;
public:
    Multiplication(int n) {
        this->n = n;
        mat_alloc(mata);
        mat_alloc(matb);
        mat_alloc(matc);
        mat_init(matb);
        mat_init(matc);
        rm_mult();
        //DEBUG
        //printmat(mata);
        //printmat(matb);
        //printmat(matc);
    }
public:
    /*Metode om multipliserer A = BC.
     *Beregner også tiden og skriver ut på skjerm
     */
    int rm_mult(){
        start = clock();
        for(int j = 0; j<n; j++){
            for(int i = 0; i<n; i++){
                for(int k = 0; k<n; k++){
                    mata[i][j] += matb[i][k]*matc[k][j];
                }
            }
        }
        finish = clock();
        cout<<"Time (Multiplication): "<<((finish-start)/(CLOCKS_PER_SEC/1000))<<endl;
        return 0;
    }
private:
    /*Setter random matriseverdier ved hjelp av
     *ran0() metoden fra lib.cpp
     *Usikker på om jeg brukte denne metoden riktig.
     *Den ser ut til å generere to like matriser, men
     *det har forsåvidt ikke noe å si.
     */
    int mat_init(double**&temp) {
        long *t,s;
        s = 0;
        t = &s;
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++) {
                temp[i][j] = ran0(t);
            }
        }
        return 0;
    }
    /*Dynamisk minneallokering for matrisene
     *mata, matb, matc.
     */
    int mat_alloc(double**&temp) {
        temp = new double*[n];
        for(int i = 0; i<n; i++) {
            temp[i] = new double[n];
        }
        return 0;
    }
    int printmat(double**&temp) {
        cout<<"Matrix: "<<endl;
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++) {
                cout<<temp[i][j];
            }
            cout<<endl;
        }
        return 0;
    }
};
/*Denne klassen brukes i oppgave a)-d)
 *Kunne slettet en del av print metodene,
 *men jeg brukte de som debugging.
 *Oppretter vektorer og matrisen i konstruktøren.
 *Initialiserer disse også i konstruktøren.
 *Kaller metodene solver() fra main() for å
 *utføre min algoritme.
 *Kaller metoden ludec() fra main() for å kjøre metodene
 *fra lib.cpp
 *Noen av metodene er merket updated osv.. Det er bare
 *fordi jeg først oppdaterte koden ved å bruke vektorer
 *med lengde n+2 i stedet for n.
 */
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
        //Vektorene går fra 0 -> n+1
        a = new double[n+2];
        b = new double[n+2];
        c = new double[n+2];
        f = new double[n+2];
        v = new double[n+2];
        matrix = new double*[n];
        for(int i = 0; i<n;i++){
            matrix[i] = new double[n];
        }
        init_vec(a,-1.0);
        init_vec(b,2.0);
        init_vec(c,-1.0);
        init_vec(f);
        a[0]=a[1]=a[n+1] = 0;
        c[0]=c[n+1] = 0;
        b[0]=b[n+1] = 0;
        f[0]=f[n+1] = 0;
        v[0]=v[n+1] = 0;
        init_matrix(a, b, c);
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
        lubksb(matrix, n, t, s);
        finish = clock();
        cout<<"Time (LU-dec): "<<((finish - start)/(CLOCKS_PER_SEC/1000))<<endl;
        return 0;
    }
//Updated
private:
    int init_matrix(double *a, double *b, double *c) {
        int i;
        matrix[0][0] = b[1];
        matrix[0][1] = c[1];
        for(i = 1; i<n-1; i++) {
            matrix[i][i-1] = a[i+1];
            matrix[i][i] = b[i+1];
            matrix[i][i+1] = c[i+1];
        }
        matrix[i][i-1] = a[i+1];
        matrix[i][i] = b[i+1];
        return 0;
    }
//Not updated
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
//Updated
public:
    int print(double *temp) {
        for(int i = 0; i<n+2; i++)
            cout<<temp[i]<<endl;
        return 0;
    }
//Updated
public:
    int print_f(){
        for(int i = 0; i<n+2; i++){
            cout<<f[i]<<endl;
        }
        return 0;
    }
//Updated
public:
    int print(){
        for(int i = 0; i<n+2; i++)
            cout<<v[i]<<endl;
        return 0;
    }
//Updated
private:
    /*Metode for å initialisere f(x_i)*/
    int init_vec(double *temp){
        for(int i = 1; i<n+2; i++){
            temp[i] = h*h*100*exp(-10*i*h);
        }
        return 0;
    }
//Updated
private:
    void init_vec(double *temp,double value) {
        for(int i = 0; i<n+2; i++) {
            temp[i] = value;
        }
    }
//
public:
    int solver(){
        start = clock();
        int i;
        double *temp;
        double btemp;
        temp = new double[n];
        //Forward
        btemp = b[1];
        //cout<<"F(0) = "<<f[0]<<endl;
        //v[1] = f[1]/btemp;
        //cout<<"v(0) = "<<v[0]<<endl;

        for(i = 2; i<n+1; i++){
            temp[i] = c[i-1]/btemp;
            btemp = b[i]-a[i]*temp[i];
            v[i] = (f[i]-a[i]*v[i-1])/btemp;
        }
        //Backward
        v[n] = f[n]/b[n];
        for(i = n-1; i>0; i--) {
            v[i] -=temp[i+1]*v[i+1];
        }
        finish = clock();
        //double diff = ((double)finish-(double)start)/1000000;
        //double seconds = (double)diff/CLOCKS_PER_SEC;
        //cout<<"time: "<<seconds<<endl;
        cout<<"Time (Algorithm): "<<((finish - start)/((double)CLOCKS_PER_SEC/1000))<<endl;

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
        max = -11;
        for(int i = 0; i<n+2; i++){
            file<<setiosflags(ios::showpoint | ios::uppercase);
            file<<setw(15)<<setprecision(8)<<i*h;
            file<<setw(15)<<setprecision(8)<<v[i];
            temp = 1-(1-exp(-10))*i*h-exp(-10*i*h);
            file<<setw(15)<<setprecision(8)<<temp<<endl;
            if(i == 0 || i == n+1) {
                continue;
            }
            if(temp == 0)
                epsilon = -10;
            else
                epsilon = log10(abs((v[i]-temp)/temp));
            //cout<<"Underveis: "<<epsilon<<endl;
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

    /*Oppretter et objekt av klassen Matrix
    *og en peker til dette objektet.
    *Kan kalle metodene som er deklarert
    *som public
    */
    Matrix *m = new Matrix(n, v0, v1);

    /*Oppretter objekt til oppgave e)
    *og kjører kall på matrisemultiplikasjon fra
    *konstruktøren. Kommenter vekk for å ikke gjøre
    *oppgave e)
    */
    //Multiplication *mp = new Multiplication(n);
    //Kjører min algoritme
    m->solver();
    //Debugging (se bort i fra denne)
    //m->print();
    //Skriver ut til fil
    //m->write_tofile("pdata1.dat");
    //Debugging (se bort i fra denne)
    //m->printmatrix();
    //Kjører LU-dekomposisjon fra lib.cpp
    m->ludec();
    return 0;
}
