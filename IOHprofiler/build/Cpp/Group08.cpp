#include "../../src/Template/Experiments/IOHprofiler_experimenter.hpp"

double generate_random(rangeStart, rangeEnd) {
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_int_distribution<double>  distr(range_from, range_to);

	std::cout << distr(generator) << '\n';
}

double[] add_array(double array1[], double array2[], int size);
	double result[size];
	for (int i = 0; i < size; ++i) {
		result[i] = array1[i] + array2[i];
	}
	return result;
}

void copy_array(double input[], double output[], int size) {
	for (int i = 0; i < size; ++i) {
		output[i] = input[i];
	}
}

double diff_objective(double shark[], int M; IOHprofiler_problem<double>> problem) {
	double delta;
	double xLow, xHigh, yLow, yHigh;
	const int size = sizeof(shark)/sizeof(shark[0]));
	double sharkCopy[size];
	double range = problem->IOHprofiler_get_upperbound()[m-1] - problem->IOHprofiler_get_lowerbound()[m-1];
	delta = range / 1000000;
	xLow = shark[m] - delta;
	xHigh = shark[m] + delta;
	copy_array(shark, sharkCopy, size);
	yLow = problem->evaluate(std::vector<double>(sharkCopy, size);

	yHigh = problem->evaluate(std::vector<double>(sharkCopy, size);
}



void shark_smell_search(std::shared_ptr<IOHprofiler_problem<double>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	std::vector<double> x(problem->IOHprofiler_get_number_of_variables());
	double y;

	const int N = 5; //Population size
	const int M = problem->IOHprofiler_get_number_of_variables(); //dimensions
	const int O = 12;
	const int kMax = 7; //Maximum amount of stages

	//TODO make vector
	double eta = 0.9; //Velocity multiplier
	double beta = 4.0; //Velocity ratio limiter
	double alpha = 0.1; //Inertia coefficient

	double positions[kMax+1][N][M];
	double positionsTemp[kMax+1][N][M];
	double velocities[kMax][N][M];
	double rotational[kMax+1][N][O][M];
	int k = 1; //Current stage

	double R1;
	double R2;
	double R3[M];
	double temp;

	/*Get the lower and upper bound of each dimension's search domain*/
	std::vector<double> lowerBounds = problem->IOHprofiler_get_lowerbound();
	std::vector<double> upperBound = problem->IOHprofiler_get_upperbound();

	/*Generate random starting positions*/
	for (int n = 1; n <= N; ++n) {
		for (int m = 1; m < M; ++M) {
			positions[1][n][m] = generate_random(lowerBounds[n], upperBounds[n]);
			velocities[1][n][m] = 0;
		}
	}

	/*algorithm main loop*/
	for (k; k < kMax; ++k) {
		for (int n = 1; n < N; ++n) {
			for (int m = 1; m < M; ++m) {
				R1 = generate_random(0, 1);
				R2 = generate_random(0, 1);
				velocities[k][n][m] = eta * R1 * (diff_objective(positions[k][n]));
				velocities[k][n][m] += alpha * R2 * velocities[k-1][n][m];
				if (velocities[k][n][m] > beta * velocities[k-1][n][m]) {
					velocities[k][n][m] = beta * velocities[k-1][n][m];
				}
			}
		}
		for (int n = 1; n <= N; ++n) {
			copy_array(positionsTemp[k+1][n], add_array(positions[k][n], velocities[k][n], M), M);
		}
		for (int n = 1; n <= N; ++n) {
			for (int o = 1; o <= O; ++o) {
				for (int m = 1; m <= M; ++m) {
					R3 = generate_random(-1, 1);
					rotational[k+1][n][o][m] = positionTemp[k+1][n][m] + R3 * velocities[k][n][m];
				}
			}
		}
		for (int n = 1; n <= N; ++n) {
			temp = 0;
			for (int o = 1; o <= O; ++o) {
				if (temp < position->evaluate(std::vector<double>(rotational[k+1][n][o], M)) {
					temp = position->evaluate(std::vector<double>(rotational[k+1][n][o], M));
					for (int m = 1; m <= M; ++m) {
						position[k+1][n][o][m] = rotational[k+1][n][o][m];
					}
				}
			}
			if (position->evaluate(std::vector<double>(positionTemp[k+1][n][o], M)) > position->evaluate(std::vector<double>(position[k+1][n][o], M))) {
				for (int m = 1; m <= M; ++m) {
					position[k+1][n][o][m] = positionTemp[k+1][n][o][m];
				}
			}
		}
	}
	logger->write_line(problem->loggerCOCOInfo());
} //niffo is gehutst
