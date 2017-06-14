// Black-Scholes European Option Pricing Code
// Adapted from Prof. Odegaard's Notes
// Written by Prof. Sreenivas for IE523: Financial Computing
// Black-Scholes European Option Pricing Code
// Adapted from Prof. Odegaard's Notes
// Written by Prof. Sreenivas for IE523: Financial Computing


#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include "normdist.h"          // this defines the normal distribution from Odegaard's files
using namespace std;

typedef std::chrono::high_resolution_clock Clock;


double up_factor, uptick_prob, risk_free_rate, strike_price, downtick_prob, notick_prob;
double initial_stock_price, expiration_time, volatility, R;

int no_of_divisions;
double **memo;

double max(double a, double b) {
	return (b < a) ? a : b;
}


void make_memo() {
	memo = new double *[no_of_divisions]; // makeing memo function. using dynamic array.
	for (int i = 0; i <= no_of_divisions; i++) // 
		memo[i] = new double[(2*i+1)]; // make triangular 2D array to prevent from wasting the array
}


void clean_memo(double **memo_array) {			//  cleaning memo
	for (int i = 0; i <= no_of_divisions; i++)
		for (int j = 0; j <= 2 * i + 1; j++)  // put zeros to all.
			memo_array[i][j] = -1;
}

void print_memo(double **memo_array) {
	for (int i = 0; i < no_of_divisions; i++){
		for (int j = 0; j < 2 * i + 1; j++) {  // put zeros to all.
			cout << memo_array[i][j] << " ";
		}
		cout << endl;				
	}
}


double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt; // equations
	double d2 = d1 - (sigma*time_sqrt);
	return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S*N(d1) - K*exp(-r*time)*N(d2);
};

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double european_call_option(int k, int i) {

	if (memo[k][i] != -1) {  // if memo is already written
		return memo[k][i]; //recall memo
	}
	else if (k >= no_of_divisions)
	{
		memo[k][i] = max(0.0, (initial_stock_price*pow(up_factor, ((double)(i - k)))) - strike_price);//store the value 
		return memo[k][i];
	}
	else
	{
		memo[k][i] = ((uptick_prob*european_call_option(k + 1, i + 2) + notick_prob*european_call_option(k + 1, i + 1) + downtick_prob*european_call_option(k + 1, i )) / R);
		return memo[k][i]; // reculsive algorithm. 
	}
}


double european_put_option(int k, int i) {
	if (memo[k][i] != -1) {
		return memo[k][i];
	}
	else if (k == no_of_divisions) {
		memo[k][i] = max(0.0, strike_price - (initial_stock_price*pow(up_factor, ((double)(i - k)))));
		return memo[k][i];
	}
	else {
		memo[k][i] = ((uptick_prob*european_put_option(k + 1, i + 2) + notick_prob*european_put_option(k + 1, i + 1) + downtick_prob*european_put_option(k + 1, i)) / R);
		return memo[k][i];
	}
}


int main(int argc, char* argv[])
{

	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%d", &no_of_divisions);
	sscanf(argv[3], "%lf", &risk_free_rate);
	sscanf(argv[4], "%lf", &volatility);
	sscanf(argv[5], "%lf", &initial_stock_price);
	sscanf(argv[6], "%lf", &strike_price);


	make_memo();
	clean_memo(memo);


	up_factor = exp(volatility*sqrt(2 * expiration_time / ((double)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((double)no_of_divisions));
	uptick_prob = pow(((sqrt(R) - (1 / (sqrt(up_factor)))) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	downtick_prob = pow(((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	notick_prob = 1 - (downtick_prob + uptick_prob);


	cout << "Recursive Binomial European Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "--------------------------------------" << endl;

	auto t_1 = Clock::now();
	double call_price = european_call_option(0, 0);
	auto t_2 = Clock::now();

	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "Trinomial Price of an European Call Option = " << call_price << "			"<<time_diff_call<<" sec"<<endl;
	cout << "Call Price according to Black-Scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;

	clean_memo(memo); // clean the memo so that use the array again

	cout << "--------------------------------------" << endl;
	t_1 = Clock::now();
	double put_price = european_put_option(0, 0);
	 t_2 = Clock::now();

	double time_diff_put = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "Trinomial Price of an European Put Option = " << put_price <<"			" << time_diff_put << " sec" << endl;
	cout << "Put Price according to Black-Scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	cout << endl;
	cout << "Total Time spent is " << (time_diff_call + time_diff_put) <<" sec"<< endl;
	cout << "--------------------------------------" << endl;
}
