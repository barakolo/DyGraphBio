// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "parser.hpp"
#include "propgraph.hpp"
#include <stdio.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <iostream>
#include <random>

using namespace std;
using namespace std::chrono;

//#define DEBUG_MODE

void sc_print(double* scores, uint32_t size) {
#ifdef DEBUG_MODE
	std::cout << "Scores are:" << endl;
	for (int i = 0; i < size; i++) {
		if (scores[i] > 1e-4)
			std::cout << "g" << i << ": " << scores[i] << endl;
	}
#endif
}

double diff(double* scores1, double* scores2, uint32_t size) {
	double s = 0.0;
	for (uint32_t i = 0; i < size; i++)
		s += abs(scores1[i] - scores2[i]);
	return s;
}

int main()
{
	cout << "Starting..." << endl;
	node_map nm = parse_graph();
	cout << "Generating Graph..." << endl;
	DyGraph g(nm, true);
	node_map nm2 = parse_graph();
	PropGraph g2(nm2, true);
	list<Node*> prior;
	/*
	for (int j = 0; j < 100; j++)
		prior.push_back(nm[1 + j]);*/
	uint32_t HD_GENE = 2306;
	prior.push_back(nm[HD_GENE]);
	g.set_prior(prior);
	g2.set_prior(prior);
	cout << "Propagating power static:" << endl;
	double* scores2 = new double[nm.size()];
	high_resolution_clock::time_point tp1 = high_resolution_clock::now();
	g2.propagate(scores2);
	high_resolution_clock::time_point tp2 = high_resolution_clock::now();
	duration<double> ts = duration_cast<duration<double>>(tp2 - tp1);
	cout << "PropStatic: propagation time for one:" << ts.count() << endl;
	sc_print(scores2, nm.size());
	cout << "Done for power static:" << endl;
	cout << "Propagating dyn-fp scratch:" << endl;
	double* scores = new double[nm.size()];
	high_resolution_clock::time_point tfp1 = high_resolution_clock::now();
	g.propagate(scores);
	high_resolution_clock::time_point tfp2 = high_resolution_clock::now();
	ts = duration_cast<duration<double>>(tfp2 - tfp1);
	cout << "ScratchFP: propagation time for one:" << ts.count() << endl;
	sc_print(scores, nm.size());
	cout << "Done for dyn-fp scratch:" << endl;
	cout << "cur diff is:" << diff(scores, scores2, nm.size()) << endl;
	uint32_t highest_diff_gene = 0;
	uint32_t highest_diff_gene2 = 0;
	double highest_diff = 0;
	double* R = new double[nm.size()];
	double* P = new double[nm.size()];
	R = std::move(g.get_R());
	P = std::move(g.get_P());

	high_resolution_clock::time_point t1, t2; t1 = t2;
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	duration<double> time_span2 = duration_cast<duration<double>>(t2 - t1);
	high_resolution_clock::time_point tstart = high_resolution_clock::now();
	
	random_device rd;
	uniform_int_distribution<int> dist(0, 18000);
	mt19937 engine2(rd());
	for (uint32_t i = 8; i < 18000; i++) {
		uint32_t j = dist(engine2);
		if (i == HD_GENE || j == HD_GENE)
			continue;
		std::cout << "Start for dyn-fp:" << endl;
		g.set_R(R);
		g.set_P(P);
		t1 = high_resolution_clock::now();
		g.del_node(i, scores);
		//g.del_node(j, scores);
		t2 = high_resolution_clock::now();
		time_span += duration_cast<duration<double>>(t2 - t1);
		std::cout << "DynFP: Current neto-time for dels: " << time_span.count() << " seconds." << endl;
		g.enable_node(i); // restore node edges to valid and the sums accordingly.
		sc_print(scores, nm.size());
		double cur_diff = diff(scores, P, nm.size());
		if (highest_diff < cur_diff) {
			highest_diff = cur_diff;
			highest_diff_gene = i;
			highest_diff_gene2 = j;
		}
		std::cout << "Done for dyn-fp:" << endl;
		std::cout << "Cur max gene:" << highest_diff_gene << " " << highest_diff_gene2 << endl;
		
		
		PropGraph g3(nm2, true, i);
		g3.set_prior(prior);
		t1 = high_resolution_clock::now();
		g3.propagate(scores2);
		t2 = high_resolution_clock::now();
		time_span2 += duration_cast<duration<double>>(t2 - t1);
		cout << "PropStatic: Current neto-time for dels: " << time_span2.count() << " seconds." << endl;
		cout << "Curr diff:" << diff(scores, scores2, nm.size()) << endl;
		
	}

	high_resolution_clock::time_point tend = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(tend - tstart);
	std::cout << "Total time for all:" << time_span.count() << " seconds. " << endl;

	while (getchar());
    return 0;
}

