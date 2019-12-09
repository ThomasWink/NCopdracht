#ifndef GROUP08_CPP
#define GROUP08_CPP

#include "../../src/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>

/*Defenitions for multiple-order vectors*/
typedef std::vector<double> vec1;
typedef std::vector<vec1> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<vec3> vec4;

/*
@summary generates a random number.
@param rangeStart: The lower bound of the range to be generated in.
@param rangeEnd: The upper bound of the range to be generated in.
@return A random double.
*/
double generate_random(double rangeStart, double rangeEnd) {
	std::random_device rand_dev; //Random device
	std::mt19937 generator(rand_dev()); //Random generator
	//Generates a uniform distribution of real values
	std::uniform_real_distribution<double>  distr(rangeStart, rangeEnd);
	return distr(generator);
}

/*
@summary Adds the values of two arrays into a third.
@param array1: The first array to be added.
@param array2: The second array to be added.
@param result: The result array, automatic call by reference.
@param size: The size of arrays array1, array2 and result.
*/
std::vector<double> add_vector(std::vector<double> a, std::vector<double> b) {
	std::vector<double> result;
	for (int i = 0; i < a.size(); ++i) {
		result.push_back(a[i] + b[i]);
	}
	return result;
}

/*
@summary Performes a partial derivative of the objective function to x
         in the m-th dimension.
@param shark: The positions of the individuals.
@param m: The results of the optimization will be logged here.
@param M: The dimension to which the objective function should be derived.
@param problem: The problem to be optimized.
@param logger: The results of the optimization will be logged here.
@return The partial derivative that is computed here.
*/
double diff_objective(std::vector<double>& shark, int m, std::shared_ptr<IOHprofiler_problem<double>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	double delta;
	double xLow, xHigh, yLow, yHigh;
	double range = problem->IOHprofiler_get_upperbound()[m] - problem->IOHprofiler_get_lowerbound()[m];
	std::vector<double> sharkCopy (shark);
	delta = range / 1000000;
	xLow = shark[m] - delta;
	xHigh = shark[m] + delta;
	//sharkCopy = shark;
	sharkCopy[m] = xLow;
	yLow = problem->evaluate(sharkCopy);
	logger->write_line(problem->loggerCOCOInfo());
	sharkCopy[m] = xHigh;
	yHigh = problem->evaluate(sharkCopy);
	logger->write_line(problem->loggerCOCOInfo());


	//if (isnan((yHigh-yLow)/delta)){


		std::cout << "=========================" << std::endl;
		std::cout << (yHigh-yLow)/delta << std::endl;
		std::cout << yHigh << std::endl;
		std::cout << yLow << std::endl;
		std::cout << delta << std::endl;
		std::cout << xHigh << std::endl;
		std::cout << xLow << std::endl;
		std::cout << shark[m] << std::endl;

	//}

	if (isnan((yHigh-yLow)/delta)) exit(0);

	return (yHigh-yLow)/delta;
}


/*
@summary Main algorithm function.
@param problem: The problem to be optimized.
@param logger: The results of the optimization will be logged here.
*/
void shark_smell_search(std::shared_ptr<IOHprofiler_problem<double>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	std::vector<double> x(problem->IOHprofiler_get_number_of_variables());

	/*Assign constant parameters*/
	const int N = 50; //Population size
	const int M = problem->IOHprofiler_get_number_of_variables(); //dimensions
	const int O = 12; //amount of random searches
	const int kMax = 100; //Maximum amount of stages
	const double eta = 0.9; //Velocity multiplier
	const double beta = 4.0; //Velocity ratio limiter
	const double alpha = 0.1; //Inertia coefficient

	/*Assign temporary variables*/
	double r1; //Random continuous value between 0 and 1
	double r2; //Random continuous value between 0 and 1
	double R3; //Used as if it is a vector, random value between -1 and 1
	double temp1; //Arbitrary random variable
	double temp2; //Arbitrary random variable
	bool maximization;
	double a;

	/*Get the lower and upper bound of each dimension's search domain*/
	//std::vector<double> lowerBounds = problem->IOHprofiler_get_lowerbound();
	std::vector<double> lowerBounds = problem->IOHprofiler_get_lowerbound();
	std::vector<double> upperBounds = problem->IOHprofiler_get_upperbound();

	/*Vectors containing position and velocity vectors per stage per individual*/
	vec3 positions; //position of individual
	vec3 positionsTemp; //provisional position of individual
	vec3 velocities; //velocity values of individual
	vec4 rotational; //Random positions of individual, to be evaluated

	/*Allocate space on vectors*/
	//Allocate first-order vectors
	positions.resize(kMax+2);
	positionsTemp.resize(kMax+2);
	velocities.resize(kMax+2);
	rotational.resize(kMax+2);
	//Allocate second-order vectors
	for (int k = 0; k <= kMax+1; ++k) {
		positions[k].resize(N);
		positionsTemp[k].resize(N);
		velocities[k].resize(N);
		rotational[k].resize(N);
		//Allocate third-order vectors
		for (int n = 0; n < N; ++n) {
			positions[k][n].resize(M);
			positionsTemp[k][n].resize(M);
			velocities[k][n].resize(M);
			rotational[k][n].resize(O);
			//Allocate fourth-order vectors
			for (int o = 0; o < O; ++o)  {
				rotational[k][n][o].resize(M);
			}
		}
	}

	/*Generate random starting positions for individuals*/
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < M; ++m) {
			//std::cout << lowerBounds[m] << std::endl;
			positions[1][n][m] = generate_random(lowerBounds[m], upperBounds[m]);
			velocities[0][n][m] = (upperBounds[m] - lowerBounds[m])/100;
		}
	}
	/*Determine if algorithm is maximization or minimization*/
	if (problem->IOHprofiler_get_optimization_type() == 1) {
		maximization = 1;
	} else {
		maximization = 0;
	}

	/*algorithm main loop*/
	//iteration over all stages
	for (int k = 1; k <= kMax; ++k) {
		//std::cout << "=======================================" << std::endl
		//					<< "=======================================" << std::endl;
		//Calculation of movement per individual, dimension
		for (int n = 0; n < N; ++n) {
			for (int m = 0; m < M; ++m) {
				r1 = generate_random(0, 1);
				r2 = generate_random(0, 1);

				//function 1
				if (maximization) {
					a = diff_objective(positions[k][n], m, problem, logger);
				} else {
					a = -diff_objective(positions[k][n], m, problem, logger);
				}
				//if (isnan(a)){
					std::cout << a << std::endl;
					std::cout << "k-1: " <<positions[k-1][n][m] << std::endl;
					std::cout << "k-: " << positions[k][n][m] << std::endl;
					std::cout << "k: " << k << std::endl;
				//}

				velocities[k][n][m] = eta * r1 * a;
				velocities[k][n][m] += alpha * r2 * velocities[k-1][n][m];

				//TODO delete prints when finished
				if (n == 0 && m == 0){
					std::cout << "afgeleide: " << a << std::endl;
					std::cout << "V[k-1]: " << velocities[k-1][n][m] << std::endl;
					std::cout << "V[k-0]: " << velocities[k][n][m] << std::endl;
					std::cout << "X[k-1]: " << positions[k-1][n][m] << std::endl;
					std::cout << "X[k-0]: " << positions[k][n][m] << std::endl;
					std::cout << "k: " << k << std::endl;
					std::cout << "n: " << n << std::endl;
					std::cout << "m: " << m << std::endl;
				}
				//function 2
				if (abs(velocities[k][n][m]) > abs(beta * velocities[k-1][n][m])) {
					velocities[k][n][m] = beta * velocities[k-1][n][m] * velocities[k][n][m]/abs(velocities[k][n][m]);
				}
			}
		}
		//Provisional assignment of new position per individual
		for (int n = 0; n < N; ++n) {
			//function 3
			positionsTemp[k+1][n] = add_vector(positions[k][n], velocities[k][n]);
		}
		//Random local search from provisional position
		for (int n = 0; n < N; ++n) {
			//per individual, dimension
			for (int o = 0; o < O; ++o) {
				for (int m = 0; m < M; ++m) {
					R3 = generate_random(-1, 1);
					//Derivation from original algorithm / function 4
					rotational[k+1][n][o][m] = positionsTemp[k+1][n][m] + R3 * velocities[k][n][m];
				}
			}
		}
		//Definitive assignment of new position per individual (funcition 5)
		for (int n = 0; n < N; ++n) {
			//iteration over the random local search points
			for (int o = 0; o < O; ++o) {
				if (o == 0){
					temp1 = problem->evaluate(rotational[k+1][n][o]);
					logger->write_line(problem->loggerCOCOInfo());
					positions[k+1][n] = rotational[k+1][n][o];
				}
				else {
					temp2 = problem->evaluate(rotational[k+1][n][o]);
					logger->write_line(problem->loggerCOCOInfo());
					if (maximization) {
						if (temp1 > temp2) {
							temp1 = temp2;
							positions[k+1][n] = rotational[k+1][n][o];
						}
					} else {
						if (temp1 < temp2) {
							temp1 = temp2;
							positions[k+1][n] = rotational[k+1][n][o];
						}
					}
				}
			}
			temp1 = problem->evaluate(positionsTemp[k+1][n]);
			logger->write_line(problem->loggerCOCOInfo());
			temp2 = problem->evaluate(positions[k+1][n]);
			logger->write_line(problem->loggerCOCOInfo());
			if (maximization) {
				if (temp1 > temp2) {
					positions[k+1][n] = positionsTemp[k+1][n];
				}
			} else {
				if (temp1 < temp2) {
					positions[k+1][n] = positionsTemp[k+1][n];
				}
			}
		}
	}
} //main function

#endif //GROUp08_CPP
