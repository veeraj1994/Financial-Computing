

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>          // this defines the normal distribution from Odegaard's files
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
	memo = new double *[no_of_divisions];
	for (int i = 0; i <= no_of_divisions; i++)
		memo[i] = new double[(2 * i + 1)];
}


void clean_memo(double **memo_array) {			//왜 Array 크기를 0으로만들어도 Array가 생기는지 질문!
												// 왜 Array를 만들면, -1, -2, -3 에도 자리가 생기는지 질문.
	for (int i = 0; i <= no_of_divisions; i++)
		for (int j = 0; j <= 2 * i + 1; j++)  // put zeros to all.
			memo_array[i][j] = -1;
}

void print_memo(double **memo_array) {
	for (int i = 0; i < no_of_divisions; i++) {
		for (int j = 0; j < 2 * i + 1; j++) {  // put zeros to all.
			cout << memo_array[i][j] << " ";
		}
		cout << endl;
	}
}


double american_call_option(int k, int i) {

	if (memo[k][i] != -1) {
		return memo[k][i];
	}
	else if (k >= no_of_divisions)
	{
		memo[k][i] = max(0.0, (initial_stock_price*pow(up_factor, ((double)(i - k)))) - strike_price);
		return memo[k][i];
	}
	else
	{
		memo[k][i] = max(((initial_stock_price*pow(up_factor, ((double)(i-k))))- strike_price), ((uptick_prob*american_call_option(k + 1, i + 2) + notick_prob*american_call_option(k + 1, i + 1) + downtick_prob*american_call_option(k + 1, i)) / R));
		return memo[k][i];
	}
}


double american_put_option(int k, int i) {
	if (memo[k][i] != -1) {
		return memo[k][i];
	}
	else if (k == no_of_divisions) {
		memo[k][i] = max(0.0, strike_price - (initial_stock_price*pow(up_factor, ((double)(i - k)))));
		return memo[k][i];
	}
	else {
		memo[k][i] = max(((strike_price - initial_stock_price*pow(up_factor, ((double)(i - k))))), ((uptick_prob*american_put_option(k + 1, i + 2) + notick_prob*american_put_option(k + 1, i + 1) + downtick_prob*american_put_option(k + 1, i)) / R));
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
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "--------------------------------------" << endl;

	auto t_1 = Clock::now();
	double call_price = american_call_option(0, 0);
	auto t_2 = Clock::now();

	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "Trinomial Price of an American Call Option = " << call_price << "			" << time_diff_call << " sec" << endl;

	clean_memo(memo);

	t_1 = Clock::now();
	double put_price = american_put_option(0, 0);
	t_2 = Clock::now();

	double time_diff_put = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "Trinomial Price of an American Put Option = " << put_price << "			" << time_diff_put << " sec" << endl;
	
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	cout << " It seems Put-Call parity does not exists" << endl;
	cout << endl;
	cout << "Total Time spent is " << (time_diff_call + time_diff_put) << " sec" << endl;
	cout << "--------------------------------------" << endl;
}
