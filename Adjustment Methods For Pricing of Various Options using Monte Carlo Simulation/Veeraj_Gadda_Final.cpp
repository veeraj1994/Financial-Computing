
// Pricing an European Option using Simulation
// Written by Prof. Sreenivas

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include "normdist.h" 
#include <random>
#include <time.h>
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;
typedef std::chrono::high_resolution_clock Clock; // to use clock!
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price;

double put_option_monte_price = 0.0;
double call_option_monte_price = 0.0;
double adj_call_price = 0.0;
double adj_put_price = 0.0;

int no_of_trials, no_of_divisions;


double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
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

double option_price_delta_call_black_scholes(const double& S,     // spot price
	const double& K,     // Strike (exercise) price,
	const double& r,     // interest rate
	const double& sigma, // volatility
	const double& time) {  // time to maturity
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double delta = N(d1);
	return delta;
};

double option_price_delta_put_black_scholes(const double& S, // spot price
	const double& K, // Strike (exercise) price,
	const double& r,  // interest rate
	const double& sigma,
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double delta = -N(-d1);
	return delta;
}

double closed_form_down_and_out_european_call_option()
{
	// I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
	double K = (2 * risk_free_rate) / (volatility*volatility);
	double A = option_price_call_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	double B = (barrier_price*barrier_price) / initial_stock_price;
	double C = pow(initial_stock_price / barrier_price, -(K - 1));
	double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
	return (A - D*C);
}


double closed_form_down_and_out_european_put_option()
{
	// just making it easier by renaming the global variables locally
	double S = initial_stock_price;
	double r = risk_free_rate;
	double T = expiration_time;
	double sigma = volatility;
	double H = barrier_price;
	double X = strike_price;

	// Took these formulae from some online reference
	double lambda = (r + ((sigma*sigma) / 2)) / (sigma*sigma);
	double temp = 2 * lambda - 2.0;
	double x1 = (log(S / H) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double y = (log(H*H / (S*X)) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double y1 = (log(H / S) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double put_down_in = (-S*N(-x1) + X*exp(-r*T)*N(-x1 + sigma*sqrt(T)) +
		S*pow(H / S, 2 * lambda)*(N(y) - N(y1)) -
		X*exp(-r*T)*pow(H / S, temp)*(N(y - sigma*sqrt(T)) - N(y1 - sigma*sqrt(T))));

	return option_price_put_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time) - put_down_in;
}
double max(double a, double b) {
	return (b < a) ? a : b;
}

double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

double barrier(double price) {
	if (price < barrier_price)
		return 0;
	else
		return price;
}

double adjust_probability(double S) {
	if (S < barrier_price)
		return 0;
	else
		return (double)(1-(exp(-((2 * log(initial_stock_price / barrier_price)*log(S / barrier_price))) / (expiration_time*pow(volatility, 2)))));
}


void adjust_price() {
	double delta_T = expiration_time / ((double)no_of_divisions);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);

	for (int j = 0; j < no_of_trials; j++) {

		double S1 = initial_stock_price;
		double S2 = initial_stock_price;
		double S3 = initial_stock_price;
		double S4 = initial_stock_price;


		for (int i = 0; i < no_of_divisions; i++)
		{
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			S1 = S1 * exp(delta_R + delta_SD*a);
			S2 = S2 * exp(delta_R - delta_SD*a);
			S3 = S3 * exp(delta_R + delta_SD*b);
			S4 = S4 * exp(delta_R - delta_SD*b);
		}
		S1 = barrier(S1);
		S2 = barrier(S2);
		S3 = barrier(S3);
		S4 = barrier(S4);
	
		double S1_adjusted_call = adjust_probability(S1)* max(0.0, S1 - strike_price);
		double S2_adjusted_call = adjust_probability(S2)* max(0.0, S2 - strike_price);
		double S3_adjusted_call = adjust_probability(S3)* max(0.0, S3 - strike_price);
		double S4_adjusted_call = adjust_probability(S4)* max(0.0, S4 - strike_price);
	//	cout << adjust_probability(S1)*barrier(S1) << endl;
		adj_call_price += (S1_adjusted_call + S2_adjusted_call + S3_adjusted_call + S4_adjusted_call) / 4.0;
	
		double S1_adjusted_put = adjust_probability(S1)* max(0.0, strike_price - S1);
		double S2_adjusted_put = adjust_probability(S2)* max(0.0, strike_price - S2);
		double S3_adjusted_put = adjust_probability(S3)* max(0.0, strike_price - S3);
		double S4_adjusted_put = adjust_probability(S4)* max(0.0, strike_price - S4);

		adj_put_price += (S1_adjusted_put + S2_adjusted_put + S3_adjusted_put + S4_adjusted_put) / 4.0;
	}
	adj_call_price = exp(-risk_free_rate*expiration_time)*(adj_call_price / ((double)no_of_trials));
	adj_put_price = exp(-risk_free_rate*expiration_time)*(adj_put_price / ((double)no_of_trials));

}

void montecarlo() {
	double delta_T = expiration_time / ((double)no_of_divisions);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);

	for (int j = 0; j < no_of_trials; j++) {

		double S1 = initial_stock_price;
		double S2 = initial_stock_price;
		double S3 = initial_stock_price;
		double S4 = initial_stock_price;


		for (int i = 0; i < no_of_divisions; i++)
		{
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			if (S1 != 0) 	S1 = barrier(S1) * exp(delta_R + delta_SD*a);
			if (S2 != 0) 	S2 = barrier(S2) * exp(delta_R - delta_SD*a);
			if (S3 != 0) 	S3 = barrier(S3) * exp(delta_R + delta_SD*b);
			if (S4 != 0) 	S4 = barrier(S4) * exp(delta_R - delta_SD*b);

			if ((S1 == 0) && (S2 == 0) && (S3 == 0) && (S4 == 0)) break;

			}
		call_option_monte_price += (max(0.0, S1 - strike_price) +
			max(0.0, S2 - strike_price) +
			max(0.0, S3 - strike_price) +
			max(0.0, S4 - strike_price)) / 4.0;
		put_option_monte_price += (max(0.0, strike_price - S1) +
			max(0.0, strike_price - S2) +
			max(0.0, strike_price - S3) +
			max(0.0, strike_price - S4)) / 4.0;
		}
	call_option_monte_price = exp(-risk_free_rate*expiration_time)*(call_option_monte_price / ((double)no_of_trials));
	put_option_monte_price = exp(-risk_free_rate*expiration_time)*(put_option_monte_price / ((double)no_of_trials));
	}
	



int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_divisions);
	sscanf_s(argv[8], "%lf", &barrier_price);

	
	

	/*
	cout << "--------------------------------" << endl;
	cout << "European Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;//
	cout << "Number of Trials = " << no_of_trials << endl;//
	cout << "Number of Divisions = " << no_of_divisions << endl;//
	cout << "--------------------------------" << endl;
	*/

	auto t_1 = Clock::now();




	//cout << adj_call_price << "   " << adj_put_price << endl;
	adjust_price();
	montecarlo();




	auto t_2 = Clock::now();

	/*
	ad_prob = adjust_probability(val);
	adj_call_price += (1 - ad_prob)*call_temp;
	adj_put_price += (1 - ad_prob)*put_temp;
	*/
	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);





	cout << "The average Call Price by explicit simulation " << call_option_monte_price << endl;
	cout << "The call price using the (1-p)-adjustment term = " << adj_call_price << endl;
	cout << "Theoretical Call Price		(Joonha)		= " << closed_form_down_and_out_european_call_option() << endl;
	cout << "--------------------------------" << endl;
	cout << "The average Put Price by explicit simulation " << put_option_monte_price << endl;
	cout << "The put price using the (1-p)-adjustment term = " << adj_put_price << endl;
	cout << "Theoretical Put Price		(Joonha)		= " << closed_form_down_and_out_european_put_option() << endl;
	cout << "--------------------------------" << endl;
	cout << "Total  " << time_diff_call << " sec" << endl;

/*
	for (int i = 0; i < 10; i++) {
		cout<<get_uniform()<<endl;
	}
	*/
}
