// NeedleFreeImuniasation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <thread>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>

// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements

private:
	std::vector<double> v;
};

//Operator overloads

//a function to compare MVector sizes
void MVector_Size_Comparison(MVector rhs, MVector lhs)
{
	if (lhs.size() != rhs.size())
	{
		std::cout << "error: MVectors are of different sizes :" << rhs.size() << " & " << " " << lhs.size() << std::endl;
		exit(0);
	}
}

// Operator overload for "scalar*vector" 
inline MVector operator*(const double& lhs, const MVector& rhs)
{
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs*rhs[i];
	return returnVal;
}

// Operator overload for "vector+vector" 
inline MVector operator+(const MVector& lhs, const MVector& rhs)
{
	MVector_Size_Comparison(rhs, lhs);
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs[i] + rhs[i];
	return returnVal;
}

// Operator overload for "vector-vector" 
inline MVector operator-(const MVector& lhs, const MVector& rhs)
{
	MVector_Size_Comparison(rhs, lhs);
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs[i] - rhs[i];
	return returnVal;
}

// Operator overload for "vector*vector" 
inline MVector operator*(const MVector& lhs, const MVector& rhs)
{
	MVector_Size_Comparison(rhs, lhs);
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs[i] * rhs[i];
	return returnVal;
}

// Operator overload for "vector/scalar" 
inline MVector operator/(const MVector& lhs, const double& rhs)
{
	MVector returnVal(lhs);
	for (int i = 0; i < lhs.size(); i++) returnVal[i] = lhs[i] / rhs;
	return returnVal;
}

//Operator overload for throwing MVector terms
std::ostream &operator<<(std::ostream &fss, const MVector &m)
{
	fss << m[0];
	for (int i = 1; i < m.size(); i++)
	{
		fss << "," << m[i];
	}
	return fss;
}

struct MFunction
{
	virtual MVector operator()(const double& x, const MVector& y) = 0;
};

double power(double x, int p) //a function to compute powers
{
	double temp = 1;
	if (p >= 0) { for (p; p > 0; p--) { temp = temp*x; } return temp; }
	else { p = -p; for (p; p > 0; p--) { temp = temp*x; } return 1 / temp; }
}

// A Runge-Kutta ODE solver, specify a file name and it will create a .csv
int RKSolve(int steps, double a, double b, MVector &y, MFunction &f, std::string fname = "")
{
	double h = ((b - a) / steps);
	if (h < 0) return 1; //[a,b] is not a valid domain
	MVector k1, k2, k3, k4;

	if (fname != "") //the version that creates a .csv
	{
		std::ofstream file;
		file.open(fname + ".csv");
		if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }
		file.precision(17);
		file << a << "," << y << std::endl;

		for (double x = a + h; x <= b + h / 2; x += h) //we add h/2 to our upper bound to compensate for floating point error
		{
			k1 = f(x, y);
			k2 = f(x + h / 2, y + h*k1 / 2);
			k3 = f(x + h / 2, y + h*k2 / 2);
			k4 = f(x + h, y + h*k3);
			y = y + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
			file.precision(17);
			file << x << "," << y << std::endl;
		}
		file.close();
	}
	else //the regular version of the code
	{
		for (double x = a + h; x <= b + h / 2; x += h) //we add h/2 to our upper bound to compensate for floating point error
		{
			k1 = f(x, y);
			k2 = f(x + h / 2, y + h*k1 / 2);
			k3 = f(x + h / 2, y + h*k2 / 2);
			k4 = f(x + h, y + h*k3);
			y = y + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
			//std::cout << y << std::endl << std::endl;
			//std::cin.get();
		}
	}
	return 0;
}

double NDist();

// A Runge-Kutta ODE solver, specify a file name and it will create a .csv
double RKOvershoot(int interval, double a, double b, MVector &y, MFunction &f, int colisions = 0)
{
	double h = ((b-a) / std::abs(interval));
	MVector k1, k2, k3, k4;
	for (double t = a + h; t <= b + h / 2; t += h) //we add h/2 to our upper bound to compensate for floating point error
	{
		k1 = f(t, y);
		k2 = f(t + h / 2, y + h*k1 / 2);
		k3 = f(t + h / 2, y + h*k2 / 2);
		k4 = f(t + h, y + h*k3);
		y = y + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
		if (colisions == 1) { if ((rand() % 1000) > 970 && std::abs(y[0]) < 0.03*NDist()) { return b; } } //simulating the random breaks in the upper layer
		if (colisions == 1) { if ((rand() % 1000) > 980) { y[1] = y[1] / (1.0 + NDist()); y[2] = y[2] / (1.0 + NDist()); } } //simulating random colisions in all layers
		if (std::abs(y[1]) < 0.0000001) return t; //breaking if at any point in the solve the velocity is less than 1e-8
	}
	return b;
}

class NeedleFreeS : public MFunction
{
public:
	NeedleFreeS(double rho_ = 0.1, double b_ = 1.0) { rho = rho_; b = b_; }
	virtual MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(6);
		temp[0] = y[1];//x
		temp[1] = y[2];//v, x'
		temp[2] = rho*power(y[1], 2) + b; //a, v', x''
		temp[3] = y[1];//z
		temp[4] = y[2];//z'
		temp[5] = rho*2.0*y[4] * y[1]; //z''
		return temp;
	}
	void SetRho(double rho_) { rho = rho_; }
private:
	double rho, b;
};


class NeedleFree : public MFunction
{
public:
	NeedleFree(double rho_ = 1, double c_ = 0.01, double b_ = 1) { rho = rho_; c = c_; b = b_; }
	virtual MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(6);
		temp[0] = y[1];//x
		temp[1] = y[2];//v, x'
		temp[2] = rho*power(y[1], 2) * (std::exp(-c*y[0])) + b; //a, v', x''
		temp[3] = y[4];//z
		temp[4] = y[5];//z'
		temp[5] = rho*(2 * std::exp(-c*y[0])*y[4] * y[1] - c*y[3] * std::exp(-c*y[0])*power(y[1], 2));//z''
		return temp;
	}
	void SetRho(double rho_) { rho = rho_; }
private:
	double rho, c, b;
};

class NeedleFreeA : public MFunction
{
public:
	NeedleFreeA(double rho_ = 1, double c_ = 0.01, double b_ = 1, double d_ = 1) { rho = rho_; c = c_; b = b_; d = d_; }
	virtual MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(6);
		temp[0] = y[1];//x
		temp[1] = y[2];//v, x'
		temp[2] = rho*power(y[1], 2) * (std::exp(-c*y[0])) + y[1] + b; //a, v', x''
		temp[3] = y[4];//z
		temp[4] = y[5];//z'
		temp[5] = rho*(2 * std::exp(-c*y[0])*y[4] * y[1] - c*y[3] * std::exp(-c*y[0])*power(y[1], 2))+d*y[4];//z''
		return temp;
	}
	void SetRho(double rho_) { rho = rho_; }
private:
	double rho, c, b, d;
};

//A Poisson distribution generator
double PDist() 
{
	std::random_device rd;
	std::mt19937 gen(rd());

	// if an event occurs 500 times a minute on average
	// how often is it that it occurs n times in one minute?
	std::poisson_distribution<> d(100);
	double out = 0.0;
	for (int i = 0; i < 10; i++)
	{
		out += d(gen);
	}
	return out/1000.0;
}

//A normal distribution generator
double NDist()
{
	std::random_device rd;
	std::mt19937 gen(rd());

	// values near the mean are the most likely
	// standard deviation affects the dispersion of generated values from the mean
	std::normal_distribution<> d(1.25, 0.3);
	return d(gen);
}

int BVPSolve(MFunction &f, int maxNewtonSteps, double a, double b, double y0, double y1, double farBC, double &guess, double tol = 1e-8)
{
	MVector y(6);
	double phi, phidash;
	for (int i = 0; i < maxNewtonSteps; i++)
	{
		// y[0] = y, y[1] = y', y[2] = y'', y[3] = Z_1, y[4] = Z_2, y[5]= Z_3
		y[0] = y0; y[1] = y1; y[2] = guess; y[3] = 0.0; y[4] = 0.0; y[5] = 1.0;
		RKOvershoot(100, a, b, y, f); // solve IVP
		phi = y[1] - farBC; // calculate residual
		//std::cout << "Final velocity, " << phi << std::endl; //TESTING
		phidash = y[4]; // 'Jacobian' phidash = Z_1(t - > infty)
		//std::cout << "phi, " << phi << std::endl; //TESTING
		//std::cout << "phi', " << phidash << std::endl << std::endl; //TESTING
		if (std::abs(phi) < tol) return 0; // exit if converged
		guess -= phi / std::abs(phidash); // apply newton step
		//std::cout << "Current guess, "<< guess << std::endl; //TESTING
	}
	//what happens if the guess doenst converge within the max steps
	std::cout << std::endl << "Guess did not converge to within tolerance, ";
	y[0] = y0; y[1] = y1; y[2] = guess; y[3] = 0.0; y[4] = 0.0; y[5] = 1.0;
	RKSolve(100, a, b, y, f);
	if (std::abs(y[1] - farBC) <= std::abs(phi)) //cheking for convergence, making do if we have it
	{
		std::cout << "guess is convering. choosing current guess : " << guess << std::endl;
		std::cout << "Current guess has absolute error: " << std::abs(y[1] - farBC) << std::endl << std::endl;
		return 0;
	}
	else { std::cout << "guess is not converging. Current Guess is: " << guess << std::endl; return 1; } //return an error if we have divergence.
}


int PenetrationTest(MFunction &f, int maxNewtonSteps, double a, double b, double &y0, double &y1, double farBC, double &guess, double tol = 1e-8)
{
	MVector y(6);
	double phi, phidash;
	double alpha = y0, beta = y1;
	for (int i = 0; i < maxNewtonSteps; i++)
	{
		// y[0] = y, y[1] = y', y[2] = y'', y[3] = Z_1, y[4] = Z_2, y[5]= Z_3
		y[0] = y0; y[1] = y1; y[2] = guess; y[3] = 0.0; y[4] = 0.0; y[5] = 1.0;
		double t = RKOvershoot(100, a, b, y, f);
		if(t!=b) return 2; // solve IVP
		phi = y[1] - farBC; // calculate residual
		phidash = y[4]; // 'Jacobian' phidash = Z_1(t - > infty)
		if ((std::abs(phi) < tol) & (y[0] > 0.3)) return 0; // exit if converged
		guess -= phi / std::abs(phidash); // apply newton step
	}
	//what happens if the guess doenst converge within the max steps
	std::cout << std::endl << "Guess did not converge to within tolerance, ";
	y[0] = y0; y[1] = y1; y[2] = guess; y[3] = 0.0; y[4] = 0.0; y[5] = 1.0;
	RKOvershoot(100, a, b, y, f);
	if (std::abs(y[1] - farBC) <= std::abs(phi)) //cheking for convergence, making do if we have it
	{
		std::cout << "guess is convering. choosing current guess : " << guess << std::endl;
		std::cout << "Current guess has absolute error: " << std::abs(y[1] - farBC) << std::endl << std::endl;
		return 0;
	}
	else 
	{ 
		std::cout << "guess is not converging. Current Guess is: " << guess << std::endl; 

		return 1; 
	} //return an error if we have divergence.
}
#endif

int main()
{
	/*/
	for (int i = 0; i < 100; i++)
	{
		std::cout << NDist() << std::endl;
	}
	/*/
	std::string fname;
	std::cout << "input filename: ";
	std::cin >> fname;
	std::cout << std::endl;
	std::ofstream file;
	file.open(fname + ".csv");
	if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }
	int datapoints;
	std::cout << "input number of data points: ";
	std::cin >> datapoints;
	std::cout << std::endl;
	for (int i = 0; i < datapoints; i++)
	{
		double r = NDist()/1000000, V = power((NDist()-1.25),2)+NDist(); //Area, Velocity
		double A = power(r, 2);
		double a = 0, b = 0.50;//time domain (1.0 ~ infty)
		NeedleFreeA f(A,0.1,A,r);
		MVector y(6);
		y[0] = 0.0;
		y[1] = std::abs(V);
		y[2] = -1.0; // guess
		BVPSolve(f, 100, a, b, y[0], y[1], 0.0, y[2], 0.0000001);
		RKOvershoot(100, a, b, y, f, 1);
		file << std::abs(r*V) << "," << y[0] << std::endl;
		std::cout << ((i + 1.0) / datapoints) * 100.0 << "%" << std::endl;
	}
	std::cout << "Finished" << std::endl;
	file.close();
	/**/
	return 0;
}