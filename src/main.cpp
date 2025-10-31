#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <random>
#include <time.h>


#include "optimizer/hexalyoptimizer.h". // Now Hexaly --> needs to be updated 


using namespace localsolver;
using namespace std;

int main() {
	int instance = 15;

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ READ DATA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int i, j, k, b, t, p, l, s;
	//WEIGHT FACTORS--------------------------------------------------------------------------
	float c1 = 0.33f;
	float c2 = 0.33f;
	float c3 = 1 - c1 - c2;

	int MTime = 3600;

	//Define parameters-----------------------------------------------------------------------
	const int B = 4; // number of buses available
	const int N = 6; // number of mandatory stations
	const int M = 5; // number of stations in cluster
	const int S = (N - 1) * M + N; // number of Stations
	const int R = 20; // number of passenger requests
	const int R1 = int(R / 2), R2 = R - R1;

	const int C = 20; // Bus capcity
	const int xt = 60 * 20; //Frequency max value (time between two buses)
	const float pspeed = 1.0f; // passengers speed in meter per scond
	const float bspeed = 40.0f / 3.6f; //bus speed in m/s
	const int dw = 20 * 60; // threshold of individual walking time in sec
	const int d_dl = 15 * 60;// Departing LATE
	const int d_de = 5 * 60; // Departing EARLY
	const int d_ae = 15 * 60; // Arriving EARLY
	const int d_al = 5 * 60; // Arriving LATE
	const int TS = 3600 * 4.5; //Planning horizon 


	// Read in locations 
	double** passengers = new double* [R];
	for (i = 0; i < R; i++) {
		passengers[i] = new double[2];
	}
	ifstream filep("data/input/passengers" + to_string(R) + ".txt");
	i = 0;
	while (i < R) {
		filep >> passengers[i][0] >> passengers[i][1];
		i++;
	}

	double** mandatory = new double* [N];//mandatory Stops
	for (i = 0; i < N; i++) {
		mandatory[i] = new double[2];
	}
	ifstream filem("data/input/mandatory.txt");
	i = 0;
	while (i < N) {
		filem >> mandatory[i][0] >> mandatory[i][1];
		i++;
	}

	double** optional = new double* [(N - 1) * M]; // optinal stops
	for (i = 0; i < (N - 1) * M; i++) {
		optional[i] = new double[2];
	}
	ifstream fileo("data/input/optional" + to_string(M) + ".txt");
	i = 0;
	while (i < (N - 1) * M) {
		fileo >> optional[i][0] >> optional[i][1];
		i++;
	}

	// Arrival times of the passengers 
	double* arrivals = new double[R1];
	ifstream filea("data/input/arrivals" + to_string(R) + ".txt");
	i = 0;
	while (i < R1) {
		filea >> arrivals[i];
		i++;
	}
	// Departure times of the passengers 
	double* departures = new double[R2];
	ifstream filed("data/input/departures" + to_string(R) + ".txt");
	i = 0;
	while (i < R2) {
		filed >> departures[i];
		i++;
	}

	//calculate travel time using manhattan distance
	double** traveltimep = new double* [R];
	for (i = 0; i < R; i++) {
		traveltimep[i] = new double[S];
	}
	double** traveltimes = new double* [S];
	for (i = 0; i < S; i++) {
		traveltimes[i] = new double[S];
	}

	// travel times of passengers to bus stops
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < N; j++) {
			traveltimep[i][j] = (abs(passengers[i][0] - mandatory[j][0]) + abs(passengers[i][1] - mandatory[j][1])) * 1000 / pspeed;
		}

		for (int j = N; j < S; j++) {
			traveltimep[i][j] = (abs(passengers[i][0] - optional[j - N][0]) + abs(passengers[i][1] - optional[j - N][1])) * 1000 / pspeed;
		}
	}

	//filter unncessesary data (stops farther away than dw or optional stops farther away than a mandatory stop)
	for (int i = 0; i < R; i++) {
		double minp = 10000000000;
		for (int k = 0; k < N; k++) {
			if (minp > traveltimep[i][k])minp = traveltimep[i][k];
		}
		for (int j = 0; j < S; j++) {
			if (traveltimep[i][j] > dw || (j > N && traveltimep[i][j] > minp))  traveltimep[i][j] = INT32_MAX;
		}
	}

	// travel times of buses between stops 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			traveltimes[i][j] = (abs(mandatory[i][0] - mandatory[j][0]) + abs(mandatory[i][1] - mandatory[j][1])) * 1000 / bspeed;
		}
		for (int j = N; j < S; j++) {
			traveltimes[i][j] = (abs(mandatory[i][0] - optional[j - N][0]) + abs(mandatory[i][1] - optional[j - N][1])) * 1000 / bspeed;
		}
	}

	for (int i = N; i < S; i++) {
		for (int j = 0; j < N; j++) {
			traveltimes[i][j] = (abs(optional[i - N][0] - mandatory[j][0]) + abs(optional[i - N][1] - mandatory[j][1])) * 1000 / bspeed;
		}
		for (int j = N; j < S; j++) {
			traveltimes[i][j] = (abs(optional[i - N][0] - optional[j - N][0]) + abs(optional[i - N][1] - optional[j - N][1])) * 1000 / bspeed;
		}
	}


	//remove memory
	for (i = 0; i < N; i++) {
		delete mandatory[i];
	}
	delete mandatory;

	for (i = 0; i < (N - 1) * M; i++) {
		delete optional[i];
	}
	delete optional;

	for (i = 0; i < R; i++) {
		delete passengers[i];
	}
	delete passengers;
	//exit(0);

	// ---------------- Shortest route
	double short_route = 0;
	for (i = 0; i < N - 1; i++) {
		short_route += traveltimes[i][i + 1];
	}
	cout << "short route: " << round(short_route / 60 * 100) / 100 << " min" << endl;

	double maxap = -1, maxdp = -1;
	double minp = 10000000000000000;
	for (p = 0; p < R1; p++) {
		if (maxap < arrivals[p])maxap = arrivals[p];
	}
	for (p = 0; p < R2; p++) {
		if (maxdp < departures[p])maxdp = departures[p];
		if (minp > departures[p])minp = departures[p];
	}

	int T = int(TS / (2 * short_route)) ;
	T = 3;
	cout << "Max number of trips: " << T << endl;
	double maxp = max(maxap + d_al + 2 * xt * (T), maxdp + d_dl + (M) * short_route * 2 + xt * T);
	minp -= (d_de + short_route * (M - 1));
	cout << int(maxp/60) << "  " << int(minp/60) << endl;
	//minp = 70 * 60;
	//maxp = 280 * 60;

	try {
		// +++++++++++++++++++++++++++++++++++++  Declares the optimization model. ++++++++++++++++++++++++++++++++++++++++++++++++++ 
		LocalSolver localsolver;
		LSModel model = localsolver.getModel();

		// ++++++++++++++++++++++++++++++++++++++++++++ VARIABLES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		 // routing decisions for each bus b on each trip t
		vector<vector<LSExpression>> x(B);
		for (b = 0; b < B; b++) {
			vector<LSExpression> x1(T);
			for (t = 0; t < T; t++) {
				x1[t] = model.listVar(S); // number of possible stops (can be different for each bus on eeach trip) 
				x1[t].setName("x_" + to_string(b) + "_" + to_string(t));
			}
			x[b] = x1;
		}
		//*/

		// 0-1 assignment decisions: assigns a passenger to one bus on one trip and one stop
		vector<vector<vector<vector<LSExpression>>>> y(R);
		for (p = 0; p < R; p++) {
			vector<vector<vector<LSExpression>>> y1(B);
			for (b = 0; b < B; b++) {
				vector<vector<LSExpression>> y2(T);
				for (t = 0; t < T; t++) {
					vector<LSExpression> y3(S);
					for (i = 0; i < S; i++) {
						y3[i] = model.boolVar();
						y3[i].setName("y_" + to_string(p) + "_" + to_string(b) + "_" + to_string(t) + "_" + to_string(i));
					}
					y2[t] = y3;
				}
				y1[b] = y2;
			}
			y[p] = y1;
			//y[p][b][t].setName("y" + to_string(p) + "_" + to_string(b) + "_" + to_string(t));
		}

		// departure times at mandatory stops
		vector<vector<vector<LSExpression>>> d(B);
		for (b = 0; b < B; b++) {
			vector<vector<LSExpression>> d1(T);
			for (t = 0; t < T; t++) {
				vector<LSExpression> d2(N);
				for (i = 0; i < 1; i++) {
					d2[i] = model.floatVar(minp, maxp);
					d2[i].setName("d_" + to_string(b) + "_" + to_string(t) + "_" + to_string(i));
				}
				d1[t] = d2;
			}
			d[b] = d1;
		}

		// Diff in departure times
		vector<vector<vector<LSExpression>>> Dd(B);
		for (b = 0; b < B; b++) {
			vector<vector<LSExpression>> Dd1(T);
			for (t = 0; t < T; t++) {
				vector<LSExpression> Dd2(N);
				for (i = 0; i < N; i++) {
					//Dd2[i] = model.floatVar(0, xt);
					//Dd2[i].setName("Dd_" + to_string(b) + "_" + to_string(t) + "_" + to_string(i));
				}
				Dd1[t] = Dd2;
			}
			Dd[b] = Dd1;
		}

		// Departure time of passengers
		LSExpression dp[R];
		for (p = 0; p < R; p++) {
			//dp[p] = model.floatVar(0, INT32_MAX);
			//dp[p].setName("dp_" + to_string(p));
		}
		// Late departure time of passengers
		LSExpression d_late[R2];
		for (p = 0; p < R2; p++) {
			//d_late[p] = model.floatVar(0, d_al);
			//d_late[p].setName("d_late_" + to_string(p));
		}
		// Early departure time of passengers
		LSExpression d_early[R2];
		for (p = 0; p < R2; p++) {
			//d_early[p] = model.floatVar(0, d_de);
			//d_early[p].setName("d_early_" + to_string(p));
		}

		// Arrival time of passengers
		LSExpression ap[R];
		for (p = 0; p < R; p++) {
			//ap[p] = model.floatVar(0, INT32_MAX);
			//ap[p].setName("ap_" + to_string(p));
		}
		// Late arrival time of passengers
		LSExpression a_late[R1];
		for (p = 0; p < R1; p++) {
			//a_late[p] = model.floatVar(0, d_al);
			//a_late[p].setName("a_late_" + to_string(p));
		}
		// Early arrival time of passengers
		LSExpression a_early[R1];
		for (p = 0; p < R1; p++) {
			//a_early[p] = model.floatVar(0, d_ae);
			//a_early[p].setName("a_early_" + to_string(p));
		}

		// ++++++++++++++++++++++++++++++++++++++++++++++ CONSTRAINTS +++++++++++++++++++++++++++++++++++++++++++++++++++++++

		// Create distance as an array to be able to access it with an "at" operator
		LSExpression distanceArray = model.array();
		for (i = 0; i < S; i++) {
			LSExpression distanceArray1 = model.array();
			for (j = 0; j < S; j++) {
				distanceArray1.addOperand(traveltimes[i][j]);
			}
			distanceArray.addOperand(distanceArray1);
		}

		// All passengers need a bus a trip and a stop
		for (p = 0; p < R; p++) {
			LSExpression Ysum = model.sum();
			for (b = 0; b < B; b++) {
				for (t = 0; t < T; t++) {
					for (i = 0; i < S; i++) {
						Ysum += y[p][b][t][i];
					}
				}
			}
			model.constraint(Ysum == 1);
			//model.maximize(Ysum);
		}

		//Walking constraint: max walking time for each passenger
		for (p = 0; p < R; p++) {
			for (i = 0; i < S; i++) {
				if (traveltimep[p][i] == INT32_MAX) {
					//LSExpression sumY = model.sum();
					for (b = 0; b < B; b++) {
						for (t = 0; t < T; t++) {
							//sumY += y[p][b][t][i];
							model.constraint(y[p][b][t][i] == 0);
						}
					}
					//model.constraint(sumY == 0);
				}
			}
		}

		//capacity constraint: max number of passenger on each bus on each trip
		for (b = 0; b < B; b++) {
			for (t = 0; t < T; t++) {
				LSExpression Ysum = model.sum();
				for (p = 0; p < R; p++) {
					for (i = 0; i < S; i++) {
						Ysum += y[p][b][t][i];
					}
				}
				model.constraint(Ysum <= C);
			}
		}

		// mandatory stops
		for (b = 0; b < B; b++) {
			for (t = 0; t < T; t++) {
				// First and last stop
				model.constraint(x[b][t][0] == 0);
				LSExpression nc = model.count(x[b][t]);
				model.constraint(x[b][t][nc - 1] == N - 1);

				//visit each mandatory stop
				for (i = 1; i < N - 1; i++) {
					LSExpression cM = model.contains(x[b][t], i);
					model.constraint(cM == true);
					if (i > 1) model.constraint(model.indexOf(x[b][t], i - 1) < model.indexOf(x[b][t], i));
				}
			}
		}
		// Optional bus stop constraints: visit these stops only when passengers are assigned to the stop
		for (b = 0; b < B; b++) {
			for (t = 0; t < T; t++) {
				for (i = N; i < S; i++) {
					LSExpression sumY = model.sum();
					for (p = 0; p < R; p++) {
						sumY += y[p][b][t][i];
					}

					LSExpression cMX = model.contains(x[b][t], i);

					//LSExpression cond = model.iif(sumY != 0, cMX == true, 1); //MAKES INCONSISTENT if no "else" argument given 
					//LSExpression cond2 = model.iif(cMX == true, sumY != 0, 1);
					//model.constraint(cond);
					model.constraint(cMX == (sumY != 0));
					//model.constraint(cond2);
				}
			}
		}


		// def of bus departures at mandatory stops
		for (b = 0; b < B; b++) {
			for (t = 0; t < T; t++) {
				LSExpression distSelector = model.createLambdaFunction([&](LSExpression i) { return model.at(distanceArray, x[b][t][i - 1], x[b][t][i]); });
				for (i = 1; i < N; i++) {
					LSExpression stop1 = model.indexOf(x[b][t], i);
					//LSExpression stop0 = model.indexOf(x[b][t], i - 1);
					LSExpression tt01 = model.sum(model.range(1, stop1 + 1), distSelector);
					//LSExpression tt01 = model.at(distanceArray, x[b][t][i], x[b][t][i - 1]);
					//if(i==0) cout << "hey i:" << i<< "\n";
					//model.constraint(d[b][t][i] == d[b][t][0] + tt01);
					d[b][t][i] = (d[b][t][0] + tt01);
				}
			}
		}

		//Availiable departure time from first mandatory stop (trip 0 not included)
		for (b = 0; b < B; b++) {
			//model.constraint(d[b][0][0] >= 0);
			for (t = 1; t < T; t++) {
				model.constraint(d[b][t][0] >= d[b][t - 1][N - 1] + short_route);
			}
		}

		// def of dp
		for (p = 0; p < R; p++) {
			LSExpression sumD = model.sum();
			for (b = 0; b < B; b++) {
				for (t = 0; t < T; t++) {
					LSExpression distSelector = model.createLambdaFunction([&](LSExpression i) { return model.at(distanceArray, x[b][t][i - 1], x[b][t][i]); });

					for (i = 0; i < N; i++) {
						//LSExpression cond1 = model.iif(y[p][b][t][i] == 1, dp[p] = d[b][t][i], 1);
						sumD += y[p][b][t][i] * (d[b][t][i]);
						//model.constraint(cond1);
					}
					for (i = N; i < S; i++) {
						LSExpression stop = model.indexOf(x[b][t], i);
						LSExpression dpt = model.sum(model.range(1, stop + 1), distSelector);

						sumD += y[p][b][t][i] * (dpt + d[b][t][0]);
						//LSExpression cond1 = model.iif(y[p][b][t][i] == 1, dp[p] = dpt + d[b][t][0], 1);
						//model.constraint(cond1);
					}
				}
			}
			dp[p] = sumD;
			//model.constraint(dp[p] == sumD);
			//model.constraint(ap[p] > dp[p]);
		}


		// def of ap
		for (p = 0; p < R; p++) {
			LSExpression sumY = model.sum();
			for (b = 0; b < B; b++) {
				for (t = 0; t < T; t++) {
					//LSExpression sumY = model.sum();
					for (i = 0; i < S; i++) {
						sumY += y[p][b][t][i] * d[b][t][N - 1];
					}
					//LSExpression cond2 = model.iif(sumY == 1, ap[p] = d[b][t][N - 1], 1);
					//model.constraint(cond2);
				}
			}
			ap[p] = sumY;
			//model.constraint(ap[p] == sumY);
		}

		//for (p = 0; p < R; p++) {
			//model.constraint(ap[p] > dp[p]);
		//}

		//time window constraints + def of intermediary variables: amout of time late or early departure/ arrivals
		//LSExpression sumD = model.sum();
		for (p = 0; p < R1; p++) {
			a_early[p] = model.iif(ap[p] < arrivals[p], arrivals[p] - ap[p], 0);
			a_late[p] = model.iif(ap[p] > arrivals[p], ap[p] - arrivals[p], 0);
			//model.constraint(arrivals[p] - ap[p] + a_late[p] - a_early[p] == 0);

			model.constraint(model.iif(ap[p] < arrivals[p], arrivals[p] - ap[p] <= d_ae, 1));
			model.constraint(model.iif(ap[p] > arrivals[p], ap[p] - arrivals[p] <= d_al, 1));
			//model.constraint(a_early[p] <= d_ae);
			//model.constraint(a_late[p] <= d_al);

			//sumD += (arrivals[p] - ap[p] + a_late[p] - a_early[p]);
		}
		for (p = 0; p < R2; p++) {
			d_early[p] = model.iif(dp[p+R1] < departures[p], departures[p] - dp[p+R1], 0);
			d_late[p] = model.iif(dp[p+R1] > departures[p], dp[p+R1] - departures[p], 0);
			//model.constraint(departures[p] - dp[p + R1] + d_late[p] - d_early[p] == 0);
			
			model.constraint(model.iif(dp[p + R1] < departures[p], departures[p] - dp[p + R1] <= d_de, 1));
			model.constraint(model.iif(dp[p + R1] > departures[p], dp[p + R1] - departures[p] <= d_dl, 1));
			
			//model.constraint(d_early[p] <= d_de);
			//model.constraint(d_late[p] <= d_dl);
			//sumD += (departures[p] - dp[p + R1] + d_late[p] - d_early[p]);
		}
		//model.minimize(sumD);

		//*
		//frequency constraints: max differnce in departure time between two consecutive buses
		for (i = 0; i < N; i++) {
			for (b = 0; b < B; b++) {
				for (t = 0; t < T; t++) {
					Dd[b][t][i] = model.min();
					for (l = 0; l < B; l++) {
						for (j = 0; j < T; j++) {
							if (l != b || t != j) {
								Dd[b][t][i].addOperand(model.iif(d[l][j][i] > d[b][t][i], d[l][j][i] - d[b][t][i], INT64_MAX));
								//Dd[b][t][i].addOperand(model.iif(d[l][j][i] > d[b][t][i], d[l][j][i] - d[b][t][i], xt*10));
								//model.constraint(cond4);
								//LSExpression cond5 = model.iif(d[l][j][i] > d[b][t][i], Dd[b][t][i] >= d[l][j][i] - d[b][t][i], 1);
								//model.constraint(cond5);
								//model.constraint(model.iif(d[l][j][i] > d[b][t][i], Dd[b][t][i] <= d[l][j][i] - d[b][t][i], 1));
								//model.constraint(model.iif(d[l][j][i] > d[b][t][i], Dd[b][t][i] >= d[l][j][i] - d[b][t][i], 1));
							}
						}
					}
				}
			}
		}
		//*
		for (i = 0; i < N; i++) {
			for (b = 0; b < B; b++) {
				for (t = 0; t < T; t++) {
					model.constraint(model.iif(Dd[b][t][i] != INT64_MAX, Dd[b][t][i] <= xt, 1));
					//model.constraint(Dd[b][t][i] > 0);
				}
			}
		}
		//*/

		// ++++++++++++++++++++++++++++++++++++++++++++++ OBJECTIVE FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++
		LSExpression Travel = model.sum();
		for (p = 0; p < R; p++) {
			Travel += (ap[p] - dp[p]);
		}

		LSExpression Walk = model.sum();
		for (p = 0; p < R; p++) {
			for (b = 0; b < B; b++) {
				for (t = 0; t < T; t++) {
					for (i = 0; i < S; i++) {
						Walk += y[p][b][t][i] * traveltimep[p][i];
					}
				}
			}
		}
		LSExpression Diff = model.sum();
		for (p = 0; p < R1; p++) {
			Diff += (a_late[p] + a_early[p]);
		}
		for (p = 0; p < R2; p++) {
			Diff += (d_late[p] + d_early[p]);
		}

		// minimize the objective sum;
		//model.minimize(Diff * c3);
		model.minimize(Travel* c1 + Walk * c2 + Diff * c3);

		// ++++++++++++++++++++++++++++++++++++++++++++++ SOLVE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// close model, then solve
		model.close();
		
		// Parameterizes the solver. 
		localsolver.getParam().setTimeBetweenDisplays(60);
		localsolver.getParam().setTimeLimit(MTime);
		localsolver.solve();

		LSSolution sol = localsolver.getSolution();

		if (sol.getStatus() == 0) {
			LSInconsistency iis = localsolver.computeInconsistency();
			std::cout << iis.toString() << std::endl;
		}
		else if (sol.getStatus() == 1) {
			cout << "\n INFEASIBLE SOLUTION: expressions violated --> \n";
			int nbExpressions = model.getNbExpressions();
			for (int exprIndex = 0; exprIndex < nbExpressions; exprIndex++) {
				LSExpression expr = model.getExpression(exprIndex);
				if (expr.isViolated()) {
					cout << expr.toString() << " violated \n";
				}
			}

		}
		else {
			ofstream xsol_p("data/output/xsol_" + to_string(instance) + ".txt");
			ofstream ysol_p("data/output/ysol_" + to_string(instance) + ".txt");
			ofstream dsol_p("data/output/dsol_" + to_string(instance) + ".txt");
			cout << " Departure times, arrival times and journey times\n";
			for (p = 0; p < R; p++) {
				cout << "p_" << p << ":\t" << to_string(int(sol.getDoubleValue(ap[p]) / 60)) << " \t" << to_string(int(sol.getDoubleValue(dp[p]) / 60)) 
					<< "\tTraveltime: " << to_string(int((sol.getDoubleValue(ap[p]) - sol.getDoubleValue(dp[p])) / 60)) << endl;
			}
			cout << endl;
			cout << " Early/ late departurees and arrivals\n";
			for (p = 0; p < R1; p++) {
				if (sol.getDoubleValue(a_early[p]) > 0) cout << "p_" << p << ":\t" << to_string(int(sol.getDoubleValue(a_early[p]) / 60)) << " early\n";
				else cout << "p_" << p << ":\t" << to_string(int(sol.getDoubleValue(a_late[p]) / 60)) << " late\n";
			}
			for (p = 0; p < R2; p++) {
				if (sol.getDoubleValue(d_early[p]) > 0)cout << "p_" << p+R1 << ":\t" << to_string(int(sol.getDoubleValue(d_early[p]) / 60)) << " early\n";
				else cout << "p_" << p+R1 << ":\t" << to_string(int(sol.getDoubleValue(d_late[p]) / 60)) << " late\n";
			}
			cout << endl;
			for (b = 0; b < B; b++) {
				cout << "++++++++++++++++++++++++++++  Bus " << b << endl;
				xsol_p << "BUS " << b << endl;
				dsol_p << "BUS " << b << endl;
				for (t = 0; t < T; t++) {
					cout << "----------------------- trip " << t << endl << "Route: \n";
					LSCollection xRoute = x[b][t].getCollectionValue();
					double startb = sol.getDoubleValue(d[b][t][0]);
					for (i = 0; i < xRoute.count()-1; i++) {
						dsol_p << startb << "\t";
						cout << to_string(xRoute[i]) << "\t";
						xsol_p << xRoute[i] << "\t";
						startb += traveltimes[xRoute[i]][xRoute[i + 1]];
					}
					dsol_p << startb << "\t";
					cout << to_string(xRoute[xRoute.count() - 1]) << "\t";
					xsol_p << xRoute[xRoute.count() - 1] << "\t";

					cout << endl << "Timetable: \n";
					xsol_p << endl;
					dsol_p << endl;
					for (i = 0; i < N; i++) {
						cout << to_string(int(sol.getDoubleValue(d[b][t][i]) / 60)) << "\t";
					}
					cout << endl;
				}
			}
			cout << endl;

			ysol_p << "Bus\tTrip\tStop\n";
			for (p = 0; p < R; p++) {
				int bus = -1, stop = -1, trip = -1;
				for (b = 0; b < B; b++) {
					for (t = 0; t < T; t++) {
						for (i = 0; i < S; i++) {
							if (sol.getValue(y[p][b][t][i]) > 0) {
								bus = b;
								trip = t;
								stop = i;
								ysol_p << b << "\t" << t << "\t" << i << "\t";
								break;
							}
						}
						if (bus != -1) break;
					}
					if (bus != -1) break;
				}
				cout << "p_" << p << "\tto bus: " << bus << "\ttrip: " << trip << "\tstop: " << stop << "   walked: " << int(traveltimep[p][stop]/60) << endl;
				ysol_p << endl;
			}

			

			xsol_p.close();
			ysol_p.close();
			dsol_p.close();
		}
	}
	catch (const exception& e) {
		cerr << "Error occurred:" << e.what() << endl;
		return 1;
	}

	//remove memory
	delete arrivals;
	delete departures;

	for (i = 0; i < S; i++) {
		delete traveltimes[i];
	}
	for (i = 0; i < R; i++) {
		delete traveltimep[i];
	}

	delete traveltimes;
	delete traveltimep;

	return 0;
}