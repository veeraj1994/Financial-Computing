
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include "newmat.h"
#include "newmatap.h" 

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

float up_factor, probability_up, probability_down, probability_stay, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R , time_diff;


int no_of_divisions;


float max(float a, float b) {
	return (b < a) ? a : b;
}

Matrix generate_probability_matrix(Matrix &probabilities) {

	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			probabilities(i, j) = 0.0; // initializing the matrix.

	// boundary values of the probabilities need to be entered
	probabilities(1, 1) = probability_stay;
	probabilities(1, 2) = 1.0 - probability_stay; // because state 1 has only two choices, stay or down

	probabilities(2 * no_of_divisions + 1, 2 * no_of_divisions + 1) = probability_stay;
	probabilities(2 * no_of_divisions + 1, 2 * no_of_divisions) = 1.0 - probability_stay;


	for (int i = 2; i <= 2 * no_of_divisions; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++) { // do not consider the first and the last rows.
			if (probabilities(i, j) == 0) {
				if (j == i) // i = j means staying
					probabilities(i, j) = probability_stay;
				if (j == (i - 1))
					probabilities(i, j) = probability_down;
				if (j == (i + 1))
					probabilities(i, j) = probability_up;
			}
		}
	return probabilities; // return the matrix
}

float ACO(int k, int i) { // American Call Option
	if (k == no_of_divisions) { // the boundary condition. 
		return max(0.0, (initial_stock_price*pow(up_factor, ((float)i))) - strike_price);
		}else{ // since this is american option, we have to compare the values for every steps.
		return max(((initial_stock_price*pow(up_factor, ((float)i))) - strike_price), ((probability_up*ACO(k + 1, i + 1) + probability_stay*ACO(k + 1, i) + probability_down*ACO(k + 1, i - 1)) / R));
	}
}


float APO(int k, int i) { // american Put option
	if (k == no_of_divisions) { // the boundary condition. 
		return max(0.0, (strike_price - initial_stock_price*pow(up_factor, ((float)i))));
		}else{
	return max(((strike_price - initial_stock_price*pow(up_factor, ((float)i)))), ((probability_up*APO(k + 1, i + 1) + probability_stay*APO(k + 1, i) + probability_down*APO(k + 1, i - 1)) / R));
	}
}



float ACO_dynamic(Matrix transition_probability)
{
	Matrix V_t(2 * no_of_divisions + 1, 1);
	// value at expiration
	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		V_t(i, 1) = max(0.0, (initial_stock_price*pow(up_factor, i - 1 - no_of_divisions)) - strike_price);

	// value at intermediate stages
	Matrix discounted_one_step_forward_value(2 * no_of_divisions + 1, 1);
	Matrix value_if_option_is_exercised_now(2 * no_of_divisions + 1, 1);
	for (int i = no_of_divisions; i > 0; i--) {
		// going backwards from expiration to zero-time
		discounted_one_step_forward_value = (transition_probability * V_t) / R;
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			value_if_option_is_exercised_now(j, 1) = max(0, (initial_stock_price*pow(up_factor, j - 1 - no_of_divisions) - strike_price));
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			V_t(j, 1) = max(value_if_option_is_exercised_now(j, 1), discounted_one_step_forward_value(j, 1));
	}
	return (V_t(no_of_divisions + 1, 1));
}


double APO_dynamic(Matrix transition_probability)
{
	Matrix V_t(2 * no_of_divisions + 1, 1);
	// value at expiration
	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		V_t(i, 1) = max(0.0, strike_price - (initial_stock_price*pow(up_factor, i - 1 - no_of_divisions)));
	// value at intermediate stages
	Matrix discounted_one_step_forward_value(2 * no_of_divisions + 1, 1);
	Matrix value_if_option_is_exercised_now(2 * no_of_divisions + 1, 1);
	for (int i = no_of_divisions; i > 0; i--) {
		// going backwards from expiration to zero-time
		discounted_one_step_forward_value = (transition_probability * V_t) / R;
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			value_if_option_is_exercised_now(j, 1) = max(0, strike_price - (initial_stock_price*pow(up_factor, j - 1 - no_of_divisions)));
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			V_t(j, 1) = max(value_if_option_is_exercised_now(j, 1), discounted_one_step_forward_value(j, 1));
	}

	return (V_t(no_of_divisions + 1, 1));
}


int main(int argc, char* argv[])
{

	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%d", &no_of_divisions);
	sscanf(argv[3], "%f", &risk_free_rate);
	sscanf(argv[4], "%f", &volatility); // sigma
	sscanf(argv[5], "%f", &initial_stock_price);
	sscanf(argv[6], "%f", &strike_price);

	Matrix transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1);
	up_factor = exp(volatility*sqrt(2 * expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	probability_up = pow(((sqrt(R) - (1 / (sqrt(up_factor)))) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	probability_down = pow(((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	probability_stay = 1 - (probability_down + probability_up);

	generate_probability_matrix(transition_probability);

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
	cout << "--------------------------------------" << endl << endl;
	/*
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
	*/
	
	cout <<endl<< "[Dynamic Programming]" << endl;
	auto t_1 = Clock::now();
	cout << "Dynamic Pricing of an American Call Option = " << ACO_dynamic(transition_probability);
	auto t_2 = Clock::now();
	time_diff = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "		" << time_diff << " Seconds " << endl;

	 t_1 = Clock::now();
	cout << "Dynamic Pricing of an American Put Option = " << APO_dynamic(transition_probability);
	 t_2 = Clock::now();
	time_diff = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);
	cout << "		" << time_diff << " Seconds " << endl;

}
