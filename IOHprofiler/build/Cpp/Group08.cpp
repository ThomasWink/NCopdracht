#include "../../src/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>

double generate_random(double rangeStart, double rangeEnd) {
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_real_distribution<double>  distr(rangeStart, rangeEnd);

	return distr(generator);
}
std::vector<double> arr2vec(double arr[], int length) {
	std::vector<double> v;
	for (int i = 0; i <= length; ++i) {
		v.push_back(arr[i]);
	}
	return v;
}

void add_array(double array1[], double array2[], double result[], int size) {
	for (int i = 0; i <= size; ++i) {
		result[i] = array1[i] + array2[i];
	}
}

void copy_array(double input[], double output[], int size) {
	for (int i = 0; i <= size; ++i) {
		output[i] = input[i];
	}
}

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



void shark_smell_search(std::shared_ptr<IOHprofiler_problem<double>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	std::vector<double> x(problem->IOHprofiler_get_number_of_variables());
	double y;

	const int N = 50; //Population size
	const int M = problem->IOHprofiler_get_number_of_variables(); //dimensions
	const int O = 12;
	const int kMax = 100; //Maximum amount of stages

	//TODO make vector
	double eta = 0.9; //Velocity multiplier
	double beta = 4.0; //Velocity ratio limiter
	double alpha = 0.1; //Inertia coefficient

	double positions[kMax+2][N+1][M+1];
	double positionsTemp[kMax+2][N+1][M+1];
	double velocities[kMax+2][N+1][M+1];
	double rotational[kMax+2][N+1][O+1][M+1];
	double tempResult[M+1];
	int k = 1; //Current stage

	double R1;
	double R2;
	double R3;
	double temp1;
	double temp2;

	/*Get the lower and upper bound of each dimension's search domain*/
	std::vector<double> lowerBounds = problem->IOHprofiler_get_lowerbound();
	std::vector<double> upperBounds = problem->IOHprofiler_get_upperbound();

	/*Generate random starting positions*/
	for (int n = 1; n <= N; ++n) {
		for (int m = 1; m <= M; ++m) {
			positions[1][n][m] = generate_random(lowerBounds[m], upperBounds[m]);
			velocities[1][n][m] = 0;
		}
	}

	/*algorithm main loop*/
	for (; k <= kMax; ++k) {
		std::cout << "Starting round " << k << std::endl;
		for (int n = 1; n <= N; ++n) {
			for (int m = 1; m <= M; ++m) {
				R1 = generate_random(0, 1);
				R2 = generate_random(0, 1);
				velocities[k][n][m] = eta * R1 * (diff_objective(positions[k][n], m, M, problem, logger));
				velocities[k][n][m] += alpha * R2 * velocities[k-1][n][m];
				if (velocities[k][n][m] > beta * velocities[k-1][n][m]) {
					velocities[k][n][m] = beta * velocities[k-1][n][m];
				}
			}
		}
		for (int n = 1; n <= N; ++n) {
			add_array(positions[k][n], velocities[k][n], tempResult, M);
			copy_array(tempResult, positionsTemp[k+1][n], M);
		}
		for (int n = 1; n <= N; ++n) {
			for (int o = 1; o <= O; ++o) {
				for (int m = 1; m <= M; ++m) {
					R3 = generate_random(-1, 1);
					rotational[k+1][n][o][m] = positionsTemp[k+1][n][m] + R3 * velocities[k][n][m];
				}
			}
		}
		for (int n = 1; n <= N; ++n) {
			temp1 = 0;
			for (int o = 1; o <= O; ++o) {
				temp2 = problem->evaluate(arr2vec(rotational[k+1][n][o], M));
				logger->write_line(problem->loggerCOCOInfo());
				if (temp1 < temp2) {
					temp1 = temp2;
					for (int m = 1; m <= M; ++m) {
						positions[k+1][n][m] = rotational[k+1][n][o][m];
					}
				}
			}
			temp1 = problem->evaluate(arr2vec(positionsTemp[k+1][n], M));
			logger->write_line(problem->loggerCOCOInfo());
			temp2 = problem->evaluate(arr2vec(positions[k+1][n], M));
			logger->write_line(problem->loggerCOCOInfo());
			if (temp1 > temp2) {
				for (int m = 1; m <= M; ++m) {
					positions[k+1][n][m] = positionsTemp[k+1][n][m];
				}
			}
		}
	}
} //niffo is gehutst
