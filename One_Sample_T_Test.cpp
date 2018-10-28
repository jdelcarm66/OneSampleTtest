
#include "pch.h"
#include <iostream>
#include <math.h>
#include <statistics.h>

// CPP Program to implement one sample t-test. 
// based on Numerical Recipes in C
// F:\D\ctraderx\SCONX2\NRC206
using namespace std;

double gammln(double xx)
{
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005*ser / x);
}

double beta(double z, double w)
{
	double gammln(double xx);

	return exp(gammln(z) + gammln(w) - gammln(z + w));
}

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

// using namespace 

double betacf(double a, double b, double x)
{
	// void nrerror(char error_text[]);
	int m, m2;
	double aa, c, d, del, h, qab, qam, qap;

	qab = a + b;
	qap = a + 1.0;
	qam = a - 1.0;
	c = 1.0;
	d = 1.0 - qab * x / qap;
	if (fabs(d) < FPMIN) d = FPMIN;
	d = 1.0 / d;
	h = d;
	for (m = 1; m <= MAXIT; m++) {
		m2 = 2 * m;
		aa = m * (b - m)*x / ((qam + m2)*(a + m2));
		d = 1.0 + aa * d;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		h *= d * c;
		aa = -(a + m)*(qab + m)*x / ((a + m2)*(qap + m2));
		d = 1.0 + aa * d;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS) break;
	}
	if (m > MAXIT) // nrerror("a or b too big, or MAXIT too small in betacf");
		cout << "a or b too big, or MAXIT too small in betacf" << endl;
	return h;
}
double betai(double a, double b, double x)
{
	double betacf(double a, double b, double x);
	double gammln(double xx);
	// void nrerror(char error_text[]);
	double bt;

	if (x < 0.0 || x > 1.0) // nrerror("Bad x in routine betai");
		cout << "Bad x in routine betai" << endl;
	if (x == 0.0 || x == 1.0) bt = 0.0;
	else
		bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
	if (x < (a + 1.0) / (a + b + 2.0))
		return bt * betacf(a, b, x) / a;
	else
		return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

// Function to find mean. 
double Mean(double arr[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum = sum + arr[i];
	return sum / n;
}

// Function to find standard 
// deviation of given array. 
double standardDeviation(double arr[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum = sum + (arr[i] - Mean(arr, n)) *
		(arr[i] - Mean(arr, n));

	return sqrt(sum / (n - 1));
}

// Function to find t-test of 
// two set of statistical data. 
double tStatistic(double arr1[], int n)
{
	double mean1 = Mean(arr1, n);
	// float mean2 = Mean(arr2, m);
	double sd1 = standardDeviation(arr1, n);
	// float sd2 = standardDeviation(arr2, m);

	// Formula to find t-test 
	// of two set of data. 
	double t_test = (mean1) / (sd1 / sqrt(n));
	return t_test;
}

// Driver function.
int main()
{
	double arr1[] = { 0.00018, -0.46247,  -0.07640, 1.71032, 2.12357, 1.43444, 1.52432, 0.54683, 0.16138,	 0.32739,
					 3.20060,	0.09702,	1.11096, 2.81644, 1.48612, 0.63236, 0.57432, 0.61730, 0.92102, -0.42061,
					 1.45099, -0.40313,  -0.97548, 1.14233, 2.09003, 1.89080, 0.33051, 0.23936, 1.36818,   0.92570 };

	// Calculate size of first array. 
	int n = sizeof(arr1) / sizeof(double);
	double mean1 = Mean(arr1, n);
	double sd1 = standardDeviation(arr1, n);
	double bothtails, lefttail, righttail, p_value, cdf, pdf, np, fac, t, degrees_freedom, significance { 0 };
	
	// Using alglib library and their data structure
	alglib::real_1d_array sample;
	sample.setcontent(n, arr1);
	for (size_t i = 0; i < n; i++)
	{
		cout << " Algblib sample array [" << i << "] = " << sample[i] << "\t" << endl;
	}
	
	alglib::studentttest1(sample, n, 0, bothtails, lefttail, righttail);

	t = tStatistic(arr1, n);
	const double pi { 3.14159265358979324 };

	// Based on Numerical recipes the cumulative distribution function is as follow:

	cdf = 0.5*betai(0.5*n, 0.5, n / (n + sqrt((t - mean1) / sd1)));
	
	// Based on Numerical recipes the probability density function is as follow:
	np = 0.5 *(n + 1);
	fac = gammln(np) - gammln(0.5*n);
	pdf = exp(-np * log(1 + sqrt((t - mean1) / sd1) / n) + fac) / (sqrt(pi*n)*sd1);

	// Based on Numerical recipes the righ tail significance for one sample t test is as follow:

	degrees_freedom = n - 1;
	significance = 0.5 * betai(0.5*degrees_freedom, 0.5, degrees_freedom / (degrees_freedom + t * t));

	// Function call. 
	cout << "\n  The sample size is                 : \t" << n;

	cout << "\n  The sample average is              : \t" << mean1;

	cout << "\n  The sample standard desviation is  : \t" << sd1;

	cout << "\n  The t statistics is                : \t" << t ;

	cout << "\n  The both tails p Value (Alglib) asociated is   : \t" << bothtails << endl;

	cout << "\n  The left tail p Value (Alglib) asociated is    : \t" << lefttail << endl;

	cout << "\n  The right tail p Value (Alglib) asociated is   : \t" << righttail << endl;

	cout << "\n  The right tail p Value (Numerical Recipes) asociated is   : \t" << significance << endl;

	std::cout << endl << "Anything else!\n";
}
