#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#define pi atan(1)*4

using namespace std;

class Planet {
    double x0, y0, vx0, vy0; //Initial coordinates and velocities
    double mass;//kg
    double distance; //Distance from origo(sun) in astrominical units
    char* name;//Name of planet

public:
        Planet(double mass, double distance, char* name, double velocity){
        this->mass = mass;
        this->name = name;
        this->x0 = distance;
        this->y0 = 0.0;
        this->vx0 = 0.0;
        this->vy0 = kmt_to_au(velocity);//average orbital velocity in au/yr
    }
    double get_x0() {
        return x0;
    }
    double get_y0() {
        return y0;
    }
    double get_vx0() {
        return vx0;
    }
    double get_vy0() {
        return vy0;
    }
    double get_mass() {
        return mass;
    }

    /*Parameter: velocity in km/h
     *Return: velocity in au/y*
     */
    double kmt_to_au(double v){
        double au = 149.6*pow(10,6);
        double y = 3600*24*365;
        return (double) v*y/au;
    }
};

class Solver {
    int n;
    double h, t_i, t_f, x_0, y_0,E_0;
    double *q,*y,G;//q={x0,y0,vx0,vy0,x1,y1,vx1,vy1,....}
    ofstream file;
public:
    Solver(int n,Planet **p) {
        this->n = n;
        t_i = 0;
        t_f = 2*pi;
        h = (t_f-t_i)/100;
        this->G = 6.67384*pow(10,-11)/pow(1.49*pow(10,11),3)*pow(3600*24*365,2); //G in au^3 kg^-1 yr^-2
        cout<<"G: "<<G<<endl;
        q = new double[4*n];
        y = new double[4*n];
        init_q(p,n); //Putting the initial values into vector q
        while(t_i<t_f){
        rk4(p);//Runge-Kutta
        t_i +=h;
        }
        file.close();
       }
public:
    int init_q(Planet **p,int n) {
        int j = 0;
        for(int i = 0; i<4*n; i=i+4) {
            q[i] = p[j]->get_x0();
            q[i+1] = p[j]->get_y0();
            q[i+2] = p[j]->get_vx0();
            q[i+3] = p[j]->get_vy0();
            j++;
        }
        return 0;
    }
public:
    int rk4(Planet **p) {
        double *k1, *k2, *k3, *k4;

        k1 = f(q,p);
        for(int i = 0; i<4*n; i++){
            k1[i] = q[i] + k1[i]/ 2.0;
        }
        k2 = f(k1,p);
        for(int i = 0; i<4*n; i++){
            k2[i] = q[i] + k2[i]/ 2.0;
        }
        k3 = f(k2,p);
        for(int i = 0; i<4*n; i++) {
            k3[i] = q[i] + k3[i];
        }
        k4 = f(k3,p);

        y[0] = q[0];
        for(int i = 0; i<4*n; i++) {
           // y[i+1] =y[i]+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
            q[i] =q[i]+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }

        for(int i = 0; i<4*n; i+=4) {
           //cout<<y[i+1]<<" "<<y[i+2]<<" "<<y[i+3]<<" "<<y[i+4]<<endl;
           cout<<q[i+0]<<" "<<q[i+1]<<" "<<q[i+2]<<" "<<q[i+3]<<endl;
        }
        //Writing to file
        atf(q,n);
        return 0;
    }

    double* f(double* q,Planet **p){
        double* k = new double[4*n];
        int r = 0;

        for (int i = 0; i < 4*n; i+=4){
            k[i+0] = q[i+2];
            k[i+1] = q[i+3];
            if(i == 0) {
                k[i+2] = 0.0;
                k[i+3] = 0.0;
            }
            else {
                //r = sqrt(q[i+0]*q[i+0]+q[i+1]*q[i+1]);
                r= sqrt(p[1]->get_x0()*p[1]->get_x0()+p[1]->get_y0()*p[1]->get_y0());
                k[i+2] = -G*(p[0]->get_mass())/pow(r,3)*q[i+0];
                k[i+3] = -G*(p[0]->get_mass())/pow(r,3)*q[i+1];
                //k[i+2] = 4*pi*pi/pow(r,2);
                //k[i+3] = 4*pi*pi/pow(r,2);
            }
        }
        return k;
    }

    int print_array(double *temp, int n){
        cout<<"Array: "<<endl;
        for(int i = 0; i<4*n; i++){
            cout<<temp[i]<<endl;
        }
        return 0;
    }

  //Array to file
    int atf(double *temp,int n) {
        if(!file.is_open())
            file.open("/uio/hume/student-u32/andehus/Project3/data.dat",ios::trunc);
        for(int i = 4; i<4*n; i+=8) {
            file<<temp[i]<<" "<<temp[i+1]<<" "<<temp[i+2]<<" "<<temp[i+3]<<endl;
        }
//        file.close();
        return 0;
    }

};

int main()
{
    int n = 2;
    Planet **p = new Planet*[n];
    //Planet(mass,distance,name,v0)
    p[0] = new Planet(2.0*pow(10,30),0,"sun",0);
    p[1] = new Planet(6.0*pow(10,24),1,"earth",29.78);
    //p[1] = new Planet(2.4*pow(10,23),0.39, "mercury",47.87);
    //p[2] = new Planet(4.9*pow(10,24), 0.72, "venus",35.02);
    //p[3] = new Planet(6.0*pow(10,24),1,"earth",29.78);
    //p[4] = new Planet(6.6*pow(10,23),1.52, "mars",24.13);
    //p[5] = new Planet(1.9*pow(10,27), 5.20, "jupiter",13.07);
    //p[6] = new Planet(5.5*pow(10,26), 9.54,"saturn",9.69);
    //p[7] = new Planet(8.8*pow(10,25), 19.19,"uranus",6.81);
    //p[8] = new Planet(1.03*pow(10,26),30.06,"neptune",5.43);
    //p[9] = new Planet(1.31*pow(10,22),39.53, "pluto",4.7);
    Solver* s = new Solver(n,p);
    return 0;
}

