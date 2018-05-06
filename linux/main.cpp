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
typedef tuple<uint32_t, double> htup;
typedef vector<htup> hvec;
typedef hvec::iterator hiter;

void sc_print(double* scores, uint32_t size) {
#ifdef DEBUG_MODE
	cout << "Scores are:" << endl;
	for (int i = 0; i < size; i++) {
		if (scores[i] > 1e-4)
			cout << "g" << i << ": " << scores[i] << endl;
	}
#endif
}

double diff(double* scores1, double* scores2, uint32_t size) {
	double s = 0.0;
	for (uint32_t i = 0; i < size; i++)
		s += abs(scores1[i] - scores2[i]);
	return s;
}

void insert_in_place(hvec& hv, htup ht) {
	hiter hi = hv.begin();
	while (hi != hv.end() && get<1>(ht) > get<1>(*hi))
		hi++;
	hv.insert(hi, ht);
}


int main()
{
    uint32_t HD_GENE = 2306; // conform with given text file indices, it is not geneid.
    uint32_t HIGHEST_SIZE = 30;
	cout << "Starting..." << endl;

    // Parsing the graph and generating node_map.
	node_map nm = parse_graph();
	cout << "Generating Graph..." << endl;
	
    // Creating DyGraph and PropGraph objects.
    DyGraph g(nm, true);
	node_map nm2 = parse_graph();
	PropGraph g2(nm2, true);
	
    // Setting the prior nodes.
    list<Node*> prior;
    prior.push_back(nm[HD_GENE]);
	g.set_prior(prior);
	g2.set_prior(prior);

    // Starting Network Propagation for Power-Iteration and DynamicFP.
    cout << "Propagating dyn-fp scratch:" << endl;
	double* scores = new double[nm.size()];
	high_resolution_clock::time_point tfp1 = high_resolution_clock::now();
	g.propagate(scores);
	high_resolution_clock::time_point tfp2 = high_resolution_clock::now();
	duration<double> ts = duration_cast<duration<double>>(tfp2 - tfp1);
	cout << "ScratchFP: propagation time for one:" << ts.count() << endl;
	sc_print(scores, nm.size());
	cout << "Done for dyn-fp scratch:" << endl;

    cout << "Propagating power static:" << endl;
	double* scores2 = new double[nm.size()];
	high_resolution_clock::time_point tp1 = high_resolution_clock::now();
	g2.propagate(scores2);
	high_resolution_clock::time_point tp2 = high_resolution_clock::now();
    ts = duration_cast<duration<double>>(tp2 - tp1);
	cout << "PropStatic: propagation time for one:" << ts.count() << endl;
	sc_print(scores2, nm.size());
    cout << "Done for power static:" << endl;
	
	cout << "cur err-diff is:" << diff(scores, scores2, nm.size()) << endl;

    // Saving current R, P and scores.
    double* R = new double[nm.size()];
	double* P = new double[nm.size()];
	double* scores_old = new double[nm.size()];
	double* Rptr = g.get_R();
	double* Pptr = g.get_P();
    for (uint32_t i=0; i<nm.size(); i++) {
        R[i] = Rptr[i];
        P[i] = Pptr[i];
        scores_old[i] = scores[i];
    }

    // Init time variables.
	high_resolution_clock::time_point t1, t2; t1 = t2;
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	duration<double> time_span2 = duration_cast<duration<double>>(t2 - t1);
	high_resolution_clock::time_point tstart = high_resolution_clock::now();

    // Init variables for maximal tuple-ko's + highest list.
    uint32_t maxg1 = -1;
    uint32_t maxg2 = -1;
    double maxd = 0.0;
    hvec highest_list;

    // Init random devices.
    random_device rd;
	uniform_int_distribution<int> dist(0, 18000);
	mt19937 engine2(rd());

	for (uint32_t i = 0; i < nm.size(); i++) {
        // Don't remove the prior.
        if (i == HD_GENE) 
            continue;
		cout << "Start for dyn-fp:" << endl;

        // Setting R,P as the first propagation result.
		g.set_R(R);
		g.set_P(P);
		t1 = high_resolution_clock::now();

        // Knockout gene no. i.
		g.del_node(i, scores);
		t2 = high_resolution_clock::now();
		time_span += duration_cast<duration<double>>(t2 - t1);
		cout << "DynFP: Current neto-time for dels: " << time_span.count() << " seconds." << endl;
        
        // restore node edges to valid and the sums accordingly.
        g.enable_node(i); 

        // Add to the highest list.
		double curd = diff(scores_old, scores, nm.size());
        if (highest_list.size() >= HIGHEST_SIZE) {
	        double maxd = get<1>(highest_list.at(0));
	        if (curd > maxd) {
	        		highest_list.erase(highest_list.begin());
            		insert_in_place(highest_list, make_pair(i, curd));
            	}
        }
        else {
        		insert_in_place(highest_list, make_pair(i, curd));
        }

		sc_print(scores, nm.size());
		cout << "Max gene list :" << endl;
		for (vector<tuple<uint32_t, double>>::const_iterator it = highest_list.begin(); it != highest_list.end(); it++)
			cout << get<0>(*it) << ", maxd= " << get<1>(*it) << endl;
		cout << "Done for dyn-fp:" << endl;
		

        // Below in comment this is the original power-iteration starting over again on removed node graph.
        // We only measure nelow the time to propagate over the new KO graph.
		/*PropGraph g3(nm2, true, i);
		g3.set_prior(prior);
		t1 = high_resolution_clock::now();
		g3.propagate(scores2);
		t2 = high_resolution_clock::now();
		time_span2 += duration_cast<duration<double>>(t2 - t1);
		cout << "PropStatic: Current neto-time for dels: " << time_span2.count() << " seconds." << endl;
		cout << "Curr diff:" << diff(scores, scores2, nm.size()) << endl;
		*/
	}

	high_resolution_clock::time_point tend = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(tend - tstart);
	cout << "Total time for all:" << time_span.count() << " seconds. " << endl;

    return 0;
}

