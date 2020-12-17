#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
using namespace std;


int iterationslimit;
double error,u0,u1;
double u;
double funcion1;
double funcion2;
double t0,h,w,y,epsilon,pi,eccentricity,hmax,hmin,l1,l2,l3,l4,l5,l6,x1,x2, it,w0;
double timeinterval,R,w1,w2,k1,k2,k3,k4,k5,k6,raiz1,raiz2,emax,q, pendiente, m;
int i, option,j,k,l;
bool recta;
//Auxiliar functions

double raiz(double t, double estimate, double e)//Solution of Kepler equation
{
	iterationslimit=0;
	error=1.0;
	u0=t;
	while(error>=estimate&&iterationslimit<=100)
	{
		u1=u0-((t-u0+(e*sin(u0)))/((e*cos(u0))-1.0));
		error=abs(u1-u0);
		u0=u1;
		iterationslimit=iterationslimit+1;
	};
	return u0;
};

double ratio(double e, double t, double estimate)//Radio of primaries
{ 
	if (e==0.0){
		return (1.0/2.0);
	}else{
		u=raiz(t,estimate,e);
		return (1.0/2.0)*(1.0-(e*cos(u)));
	}
};

double f(double a)//Velocity of the system (first equation)
{
	funcion2=a;
	return funcion2;
};

double g(double z, double e, double estimate, double t)//Aceleration of the system (second equation)
{
	funcion1=-1.0*(z/pow(pow(z,2.0)+pow(ratio(e,t,estimate),2.0),(3.0/2.0)));
	return funcion1;
};



//Boddy
int main()
{

//Variables declaration

pi=4.0*atan(1.0);
emax=1.0e-13;
epsilon=1.0e-21;
eccentricity=0.8;
q=0.8;

cout<<"1. Poincaré Map"<<"\n";
cout<<"2. Numerical Solution"<<"\n";
cout<<"3. Exit"<<"\n";

cin>>option;

if(option==3)
{
	cout<<"End"<<"\n";
}else{
	if(option==2) //Numerical solution
	{
		//Data file
		ofstream archivo;
		archivo.open("exc(0,9)init(1,1).txt",ios::out);
		if(archivo.fail()){
			cout<<"No se pudo abrir el archivo."<<"\n";
			exit(1);
		}

		//Variables declaration
		t0=0.0;
		h=0.006;
		hmax=0.006;
		hmin=0.001;
		i=0;l=0;
		w=1.0;//initial value of position
		y=1.0;//initial value of velocity

		//Method
		do{
			k1 = h*f(y);
			l1 = h*g(w,eccentricity,epsilon,t0);
			k2 = h*f(y+((1.0/4.0)*l1)); 
			l2 = h*g(w+(k1/4.0),eccentricity,epsilon,t0+((1.0/4.0)*h)); 
			k3 = h*f(y+((3.0/32.0)*l1)+((9.0/32.0)*l2));
			l3 = h*g(w+((3.0/32.0)*k1)+((9.0/32.0)*k2),eccentricity,epsilon,t0+((3.0/8.0)*h));
			k4 = h*f(y+((1932.0/2197.0)*l1)-((7200.0/2197.0)*l2)+((7296.0/2197.0)*l3));
			l4 = h*g(w+((1932.0/2197.0)*k1)-((7200.0/2197.0)*k2)+((7296.0/2197.0)*k3),eccentricity,epsilon,t0+((12.0/13.0)*h));
			k5 = h*f(y+((439.0/216.0)*l1)-(8.0*l2)+((3680.0/513.0)*l3)-((845.0/4104.0)*l4));
			l5 = h*g(w+((439.0/216.0)*k1)-(8.0*k2)+((3680.0/513.0)*k3)-((845.0/4104.0)*k4), eccentricity, epsilon,t0+h);
			k6 = h*f(y-((8.0/27.0)*l1)+(2.0*l2)-((3544.0/2565.0)*l3)+((1859.0/4104.0)*l4)-((11.0/40.0)*l5));
			l6 = h*g(w-((8.0/27.0)*k1)+(2.0*k2)-((3544.0/2565.0)*k3)+((1859.0/4104.0)*k4)-((11.0/40.0)*k5), eccentricity,epsilon,t0+((1.0/2.0)*h));

			//solution z
			w1 = w + (25.0*(k1/216.0))+(1408.0*(k3/2565.0))+(2197.0*(k4/4104.0))-(k5/5.0);
			w2 = w + (16.0*(k1/135.0))+(6656.0*(k3/12825.0))+(28561.0*(k4/56430.0))-(9.0*(k5/50.0))+(2.0*(k6/55.0));
			//solution \dot{z}
			x1 = y + (25.0*(l1/216.0))+(1408.0*(l3/2565.0))+(2197.0*(l4/4104.0))-(l5/5.0);
			x2 = y + (16.0*(l1/135.0))+(6656.0*(l3/12825.0))+(28561.0*(l4/56430.0))-(9.0*(l5/50.0))+(2.0*(l6/55.0));
			
			if((abs(w1-w2))>(abs(x1-x2))){
				R = abs(w1-w2);
				//cout<<"R1: "<<R<<"\n";
			}else{
				R = abs(x1-x2);
				//cout<<"R2: "<<R<<"\n";
			}

			if(R>emax)
			{
				h=q*pow(emax/R,1.0/6.0)*h;
				l=l+1;
			}else{
				if (R==0.0){
					//cout<<"i: "<<i<<"\n";
					//cout<<"Time: "<<t0<<"\n";
					archivo<<t0<<";"<<w<<";"<<y<<"\n";			
					t0 = t0+h;
					w = w1;
					y = x1;
					i = i+1;
				}else{
					h=q*pow(emax/R,1.0/6.0)*h;
					//cout<<"i: "<<i<<"\n";
					//cout<<"Time: "<<t0<<"\n";
					archivo<<t0<<";"<<w<<";"<<y<<"\n";
					t0 = t0+h;
					w = w1;
					y = x1;
					i = i+1;
				}
			}
		}while(t0<=(64*pi)&&i<1000000);
		cout<<"Eccentricity: "<<eccentricity<<"\n";
		cout<<"Solution iterations: "<<i<<"\n";cout<<"Times that recalculate h: "<<l<<"\n";cout<<"Time: "<<t0<<"\n";
		archivo.close();


	}else{
		if(option==1)//Poincare Map
		{
			//Data file
			ofstream archivo;
			archivo.open("PMexc(0,8)slope(1).txt",ios::out);
			if(archivo.fail()){
				cout<<"No se pudo abrir el archivo."<<"\n";
				exit(1);
			}
			//Variables declaration
			
			h=0.006;
			hmax=0.006;
			hmin=0.001;
			pendiente=1.0;
			m=1000;
			w0=-1.0;
			cout<<"Eccentricity: "<<eccentricity<<"\n";
			do
			{
				for(int j=0; j<=m; j=j+1)
				{
					//cout<<"j: "<<j<<"\n";
					if(j==0)
					{
						w=w0;//initial value of position
						y=pendiente*w;//initial value of velocity	
					}
					//cout<<"Time: "<<t0<<", z: "<<w<<", dz: "<<y<<", h: "<<h<<"\n";
					archivo<<t0<<";"<<w<<";"<<y<<"\n";
					t0=0.0;
					i=0;
					//Method
					do{
						k1 = h*f(y);
						l1 = h*g(w,eccentricity,epsilon,t0);
						k2 = h*f(y+((1.0/4.0)*l1)); 
						l2 = h*g(w+(k1/4.0),eccentricity,epsilon,t0+((1.0/4.0)*h)); 
						k3 = h*f(y+((3.0/32.0)*l1)+((9.0/32.0)*l2));
						l3 = h*g(w+((3.0/32.0)*k1)+((9.0/32.0)*k2),eccentricity,epsilon,t0+((3.0/8.0)*h));
						k4 = h*f(y+((1932.0/2197.0)*l1)-((7200.0/2197.0)*l2)+((7296.0/2197.0)*l3));
						l4 = h*g(w+((1932.0/2197.0)*k1)-((7200.0/2197.0)*k2)+((7296.0/2197.0)*k3),eccentricity,epsilon,t0+((12.0/13.0)*h));
						k5 = h*f(y+((439.0/216.0)*l1)-(8.0*l2)+((3680.0/513.0)*l3)-((845.0/4104.0)*l4));
						l5 = h*g(w+((439.0/216.0)*k1)-(8.0*k2)+((3680.0/513.0)*k3)-((845.0/4104.0)*k4), eccentricity, epsilon,t0+h);
						k6 = h*f(y-((8.0/27.0)*l1)+(2.0*l2)-((3544.0/2565.0)*l3)+((1859.0/4104.0)*l4)-((11.0/40.0)*l5));
						l6 = h*g(w-((8.0/27.0)*k1)+(2.0*k2)-((3544.0/2565.0)*k3)+((1859.0/4104.0)*k4)-((11.0/40.0)*k5), eccentricity,epsilon,t0+((1.0/2.0)*h));

						//solution z
						w1 = w + (25.0*(k1/216.0))+(1408.0*(k3/2565.0))+(2197.0*(k4/4104.0))-(k5/5.0);
						w2 = w + (16.0*(k1/135.0))+(6656.0*(k3/12825.0))+(28561.0*(k4/56430.0))-(9.0*(k5/50.0))+(2.0*(k6/55.0));
						//solution \dot{z}
						x1 = y + (25.0*(l1/216.0))+(1408.0*(l3/2565.0))+(2197.0*(l4/4104.0))-(l5/5.0);
						x2 = y + (16.0*(l1/135.0))+(6656.0*(l3/12825.0))+(28561.0*(l4/56430.0))-(9.0*(l5/50.0))+(2.0*(l6/55.0));
						
						if((abs(w1-w2))>(abs(x1-x2))){
							R = abs(w1-w2);
							//cout<<"R1: "<<R<<"\n";
						}else{
							R = abs(x1-x2);
							//cout<<"R2: "<<R<<"\n";
						}

						if(R>emax)
						{
							h=q*pow(emax/R,1.0/6.0)*h;//cout<<"h: "<<h<<"\n";
						}else{
							if (R==0.0){
								t0 = t0+h;
								w = w1;
								y = x1;
								i = i+1;
							}else{
								h=q*pow(emax/R,1.0/6.0)*h;
								t0 = t0+h;
								w = w1;
								y = x1;
								i = i+1;
							}
						}
						if(h<=0.00001)
						{t0=10.0;cout<<"Composition: "<<j<<"\n";j=m+1.0;}
					}while(t0<=(2.0*pi));
				}
				cout<<"w0: "<<w0<<"\n";
				w0=w0+0.01;
			}while(w0<=1.0);
			archivo.close();
		}	
	}
}
}

