#include "../../src/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>

/*
@summary generates a random number.
@param rangeStart: The lower bound of the range to be generated in.
@param rangeEnd: The upper bound of the range to be generated in.
@return A random double.
*/
double generate_random(double rangeStart, double rangeEnd) {
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_real_distribution<double>  distr(rangeStart, rangeEnd);

	return distr(generator);
}

/*
@summary Converts an array to a vector.
@param arr: The array that has to be converted.
@param length: Length of the array that has to be converted.
@return The converted vector.
*/
std::vector<double> arr2vec(double arr[], int length) {
	std::vector<double> v;
	for (int i = 1; i <= length; ++i) {
		v.push_back(arr[i]);
	}
	return v;
}

/*
@summary Adds the values of two arrays into a third.
@param array1: The first array to be added.
@param array2: The second array to be added.
@param result: The result array, automatic call by reference.
@param size: The size of arrays array1, array2 and result.
*/
void add_array(double array1[], double array2[], double result[], int size) {
	for (int i = 0; i <= size; ++i) {
		result[i] = array1[i] + array2[i];
	}
}

/*
@summary Copies array input into output.
@param input: The array to be copied.
@param output: The result array, automatic call by reference.
@param size: The size of arrays input and output.
*/
void copy_array(double input[], double output[], int size) {
	for (int i = 0; i <= size; ++i) {
		output[i] = input[i];
	}
}

/*
@summary Performes a partial derivative of the objective function to x in the m-th dimension.
@param shark: The positions of the individuals.
@param m: The results of the optimization will be logged here.
@param M: The dimension to which the objective function should be derived.
@param problem: The problem to be optimized.
@param logger: The results of the optimization will be logged here.
@return The partial derivative that is computed here.
*/
double diff_objective(double shark[], int m, int M, std::shared_ptr<IOHprofiler_problem<double>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	double delta;
	double xLow, xHigh, yLow, yHigh;
	double sharkCopy[M];
	double range = problem->IOHprofiler_get_upperbound()[m-1] - problem->IOHprofiler_get_lowerbound()[m-1];
	std::vector<double> temp;
	delta = range / 1000000;
	xLow = shark[m] - delta;
	xHigh = shark[m] + delta;
	copy_array(shark, sharkCopy, M);
	sharkCopy[m] = xLow;
	for (int i = 0; i < M; ++i) {
		temp.push_back(sharkCopy[i]);
	}
	yLow = problem->evaluate(temp);
	logger->write_line(problem->loggerCOCOInfo());
	sharkCopy[m] = xHigh;
	temp.clear();
	for (int i = 0; i < M; ++i) {
		temp.push_back(sharkCopy[i]);
	}
	yHigh = problem->evaluate(temp);
	logger->write_line(problem->loggerCOCOInfo());
	return (yHigh-yLow)/delta;
}


/*
@summary Main algorithm function.
@param problem: The problem to be optimized.
@param logger: The results of the optimization will be logged here.
*/
void shark_smell_search(std::shared_ptr<IOHprofiler_problem<double>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	std::vector<double> x(problem->IOHprofiler_get_number_of_variables());
	double y;

	//Variable assignments
	const int N = 50; //Populdiff_objectiveation size
	const int M = problem->IOHprofiler_get_number_of_variables(); //dimensions
	const int O = 12; //amount of random searches
	const int kMax = 100; //Maximum amount of stages


	//TODO make vector
	double eta = 0.9; //Velocity multiplier
	double beta = 4.0; //Velocity ratio limiter
	double alpha = 0.1; //Inertia coefficient

	double positions[kMax+2][N+1][M+1]; //position of individual
	double positionsTemp[kMax+2][N+1][M+1]; //provisional position of individual
	double velocities[kMax+2][N+1][M+1]; //velocity values of individual
	double rotational[kMax+2][N+1][O+1][M+1]; //Random positions of individual, to be evaluated 
	double tempResult[M+1]; //temporary copy array
	int k = 1; //Current stage

	double R1;
	double R2;
	double R3; //Used as if it is a vector
	double temp1;
	double temp2;

	/*Get the lower and upper bound of each dimension's search domain*/
	std::vector<double> lowerBounds = problem->IOHprofiler_get_lowerbound();
	std::vector<double> upperBounds = problem->IOHprofiler_get_upperbound();

	/*Generate random starting positions for individuals*/
	for (int n = 1; n <= N; ++n) {
		for (int m = 1; m <= M; ++m) {
			positions[1][n][m] = generate_random(lowerBounds[m], upperBounds[m]);
			velocities[1][n][m] = 0;
		}
	}

	/*algorithm main loop*/
	for (; k <= kMax; ++k) { //iteration over all stages
		for (int n = 1; n <= N; ++n) { //Calculation of movement per individual, dimension
			for (int m = 1; m <= M; ++m) {
				R1 = generate_random(0, 1);
				R2 = generate_random(0, 1);
				velocities[k][n][m] = eta * R1 * (diff_objective(positions[k][n], m, M, problem, logger));
				velocities[k][n][m] += alpha * R2 * velocities[k-1][n][m]; //function 1
				if (velocities[k][n][m] > beta * velocities[k-1][n][m]) { //function 2
					velocities[k][n][m] = beta * velocities[k-1][n][m];
				}
			}
		}
		for (int n = 1; n <= N; ++n) { //Provisional assignment of new position per individual 
			add_array(positions[k][n], velocities[k][n], tempResult, M); //function 3
			copy_array(tempResult, positionsTemp[k+1][n], M); 
		}
		for (int n = 1; n <= N; ++n) {//Random local search from provisional position
			for (int o = 1; o <= O; ++o) {						//per individual, dimension
				for (int m = 1; m <= M; ++m) {
					R3 = generate_random(-1, 1);
					rotational[k+1][n][o][m] = positionsTemp[k+1][n][m] + R3 * velocities[k][n][m]; //Deviation from original algorithm / function 4
				}
			}
		}
		for (int n = 1; n <= N; ++n) {//Definitive assignment of new position per individual
			temp1 = 0;
			for (int o = 1; o <= O; ++o) { //iteration over the random local search points
				temp2 = problem->evaluate(arr2vec(rotational[k+1][n][o], M));
				logger->write_line(problem->loggerCOCOInfo());
				if (temp1 < temp2) {
					temp1 = temp2;
					copy_array(rotational[k+1][n][o], positions[k+1][n], M);
				}
			}
			temp1 = problem->evaluate(arr2vec(positionsTemp[k+1][n], M));
			logger->write_line(problem->loggerCOCOInfo());
			temp2 = problem->evaluate(arr2vec(positions[k+1][n], M));
			logger->write_line(problem->loggerCOCOInfo());
			if (temp1 > temp2) {
				copy_array(positionsTemp[k+1][n], positions[k+1][n], M);
			}
		} //this for-loop is equal to function 5
	}
} //main function
