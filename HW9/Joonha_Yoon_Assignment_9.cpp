
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cstdlib>

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

float up_factor, probability_up, probability_down, probability_stay, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R, time_diff;

int no_of_divisions;


float max(float a, float b) {
	return (b < a) ? a : b;
}

float ACO(int k, int i) { // American Call Option
	if (k == no_of_divisions) {
		return max(0.0, (initial_stock_price*pow(up_factor, ((float)i))) - strike_price);
	}
	else
	{
		return max(((initial_stock_price*pow(up_factor, ((float)i))) - strike_price), ((probability_up*ACO(k + 1, i + 1) + probability_stay*ACO(k + 1, i) + probability_down*ACO(k + 1, i - 1)) / R));
	}
}


float APO(int k, int i) {
	if (k == no_of_divisions) {
		return max(0.0, (strike_price - initial_stock_price*pow(up_factor, ((float)i))));
	}
	else
	{
		return max(((strike_price - initial_stock_price*pow(up_factor, ((float)i)))), ((probability_up*APO(k + 1, i + 1) + probability_stay*APO(k + 1, i) + probability_down*APO(k + 1, i - 1)) / R));
	}
}



int main(int argc, char* argv[])
{

	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%d", &no_of_divisions);
	sscanf(argv[3], "%f", &risk_free_rate);
	sscanf(argv[4], "%f", &volatility); // sigma
	sscanf(argv[5], "%f", &initial_stock_price);
	sscanf(argv[6], "%f", &strike_price);

	
	up_factor = exp(volatility*sqrt(2*expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	probability_up =  pow(((sqrt(R)-(1/(sqrt(up_factor))))/(sqrt(up_factor)-(1/ sqrt(up_factor)))), 2);
	probability_down = pow(((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	probability_stay = 1 - (probability_down + probability_up);

	cout << "European Down and Out Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << probability_up << endl;
	cout << "Downtick Probability = " << probability_down << endl;
	cout << "Notick Probability = " << probability_stay << endl;
	cout << "--------------------------------------" << endl;

	cout << "[Recurrsive Algorithm]" << endl;
	auto t_1 = Clock::now();
	cout << "Trinomial Price of an American Call Option = " << ACO(0, 0);
	auto t_2 = Clock::now();
	time_diff = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "		" << time_diff << " Seconds " << endl;

	t_1 = Clock::now();
	cout << "Trinomial Price of an American Put Option = " << APO(0, 0);
	t_2 = Clock::now();
	time_diff = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "		" << time_diff << " Seconds " << endl;
}

