/*
Here we will place all the relevant tests.
we should test here:
1. for any graph G=(V,E,w) in {h_sapiens, yeast_graph}:
	1.a. set prior p of size 200, alpha=0.4, W=w, eps=10^{-2}, ko_set of size 1,000:
	1.b. init objects for algorithms. with eps0 = sqrt(eps/nd_max) / n. calculate times/ko_touches for them.
	1.b. for any node v in ko_set:
		2.a. calculate power_iter over "knockout-v", eps'=eps/n^2. (mark it as the real_res)
		2.a. for any algorithm A in {forward_push_scratch, dynamic_fp, power_iter}:
			3.a. if A == power_iter:
				4.a. Run power_iter over eps'=eps/c,p,alpha,W, and get optimal eps' (maximal) such that (scores - real_res) < eps.
				4.b. measure (4.a) times and node touches.
			3.b. else:
				4.c. Run A over eps1 = eps/nd_max, alpha, p, W.
				4.d. Measure (4.c) times and node touches.
*/

#include "parser.hpp"
#include "propgraph.hpp"
#include <stdio.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <iostream>
#include <random>
/*
#include "boost/program_options/parsers.hpp"
#include "boost/program_options.hpp"*/
namespace 
{ 
    const size_t ERROR_IN_COMMAND_LINE = 1; 
    const size_t SUCCESS = 0; 
    const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
}

using namespace std;
using namespace std::chrono;
ofstream fres;
//list<string> g_graphs = {"yeast_anat.txt"};
list<string> g_graphs = {"human_anat.txt"};

double diff(double* scores1, double* scores2, uint32_t size) {
	if (scores1 == NULL || scores2 == NULL)
		return -INFINITY; // Identifier for error - either first/second vector is null.
	double s = 0.0;
	for (uint32_t i = 0; i < size; i++)
		s += abs(scores1[i] - scores2[i]);
	return s;
}

double test_power(PropGraph &g, double eps, double alpha, node_map& nm, list<Node*>& prior, double* scores, uint32_t ko_node=-1, double* real_scores=NULL) {
    g.clear_prior();
    g.set_prior(prior);
    g.set_alpha(alpha);
    g.set_eps(eps);
    
    cout << "propgation started" << endl;
	high_resolution_clock::time_point tfp1 = high_resolution_clock::now();
    if (ko_node == -1)
    	g.propagate(scores);
    else {
        g.del_node(ko_node, scores);
    }
	high_resolution_clock::time_point tfp2 = high_resolution_clock::now();
    cout << "propgation ended " << endl;
	duration<double> ts = duration_cast<duration<double>>(tfp2 - tfp1);
	cout << "diff for nodes, num of nodes=" << g.get_node_cnt() << endl;
    double err = diff(scores, real_scores, g.get_node_cnt());
    fres << "power_iter-ko\t" << ko_node << "\t" << g.get_alpha() << "\t" << g.get_eps() << "\t" << err << "\t" << ts.count() << "\t" << g.get_nt() << endl;
    cout << "power_iter-ko\t" << ko_node << "\t" << g.get_alpha() << "\t" << g.get_eps() << "\t" << err << "\t" << ts.count() << "\t" << g.get_nt() << endl;
    if (ko_node != -1) {
        g.restore();
    }
	return ts.count();
}


double test_fp_scratch(DyGraph& dg, list<Node*>& prior, double* scores, double eps, int32_t ko_node=-1, double* real_scores=NULL) {
	if (eps != -1) {
		dg.set_eps(eps);
	}
    uint32_t extra_nt = 0;
    if (ko_node != -1) {
        dg.backup_state();
        extra_nt += dg.disable_node(ko_node);
    }
	dg.init_prior(prior);
    cout << "starting scratch dynamics" << endl;
	high_resolution_clock::time_point tfp1 = high_resolution_clock::now();
	dg.propagate(scores);
    cout << "ending scratch dynamics" << endl;
	high_resolution_clock::time_point tfp2 = high_resolution_clock::now();
	duration<double> ts = duration_cast<duration<double>>(tfp2 - tfp1);
	cout << "diff for nodes, num of nodes=" << dg.get_node_cnt() << endl;
    double err = diff(scores, real_scores, dg.get_node_cnt());
	fres << "fp_scratch-ko\t" << ko_node << "\t" << dg.get_alpha() << "\t" << dg.get_eps() << "\t" << err << "\t" <<
    ts.count() << "\t" << (dg.get_nt() + extra_nt) << endl;
	cout << "fp_scratch-ko\t" << ko_node << "\t" << dg.get_alpha() << "\t" << dg.get_eps() << "\t" << err << "\t" <<
    ts.count() << "\t" << (dg.get_nt() + extra_nt) << endl;
	if (ko_node != -1) {
		dg.enable_node(ko_node);
		dg.restore_state();
	}
	return ts.count();
}

double test_fp_dyn(DyGraph& dg, double* scores, double eps, uint32_t ko_node=-1, double* real_scores=NULL) {
	if (eps != -1) {
		dg.set_eps(eps);
	}
    if (ko_node != -1) {
        dg.backup_state();
    }
	high_resolution_clock::time_point tfp1 = high_resolution_clock::now();
	if (ko_node == -1)
		dg.propagate(scores);
	else {
		dg.del_node(ko_node, scores);
    }
	high_resolution_clock::time_point tfp2 = high_resolution_clock::now();
	duration<double> ts = duration_cast<duration<double>>(tfp2 - tfp1);
	double err = diff(scores, real_scores, dg.get_node_cnt());
	fres << "fp_dyn-ko\t" << ko_node << "\t" << dg.get_alpha() << "\t" << dg.get_eps() << "\t" << err << "\t" << ts.count() << "\t" << dg.get_nt() << endl;
	cout << "fp_dyn-ko\t" << ko_node << "\t" << dg.get_alpha() << "\t" << dg.get_eps() << "\t" << err << "\t" << ts.count() << "\t" << dg.get_nt() << endl;
	if (ko_node != -1) {
		dg.enable_node(ko_node);
		dg.restore_state();
	}
	return ts.count();
}

list<Node*>& get_random_prior(uint32_t prior_size, node_map& nm, uint32_t node_cnt, list<Node*>& prior) {
	map<uint32_t, Node*> m;
	random_device rd;
	uniform_int_distribution<int> dist(0, node_cnt - 1);
	mt19937 engine2(rd());

    fres << "priornodes\t";

	for (uint32_t i = 0; i < prior_size;) {
		uint32_t j = dist(engine2);
		if (m.count(j) == 0) {
			m[j] = nm[j];
            fres << i << "\t";
			i++;
		}
	}

    fres << endl;

	for (map<uint32_t, Node*>::iterator it = m.begin(); it != m.end(); it++) {
		Node* n = it->second;
		prior.push_back(n);
	}
	return prior;
}

list<uint32_t>& get_ko_list(uint32_t ko_size, list<Node*>& prior, uint32_t node_cnt, list<uint32_t>& ko_list) {
	map<uint32_t, bool> invalid_nodes;
	for (Node* n : prior) {
		invalid_nodes[n->get_node()] = true;
	}
	std::random_device rd;
	std::uniform_int_distribution<> dist(0, node_cnt - 1);
	mt19937 engine2(rd());
	
    fres << "konodes\t";

    for (uint32_t i = 0; i < ko_size;) {
		uint32_t j = dist(engine2);
		if (invalid_nodes.count(j) == 0) {
			invalid_nodes[j] = true;
            fres << j << "\t";
			ko_list.push_back(j);
			i++;
		}
	}

    fres << endl;
	return ko_list;
}

void do_test1(double alpha, double eps, uint32_t prior_size, uint32_t ko_size, string res_file = "res.txt") {
    
    fres.open(res_file);
	
    for (string g : g_graphs) {
        
        list<Node*>* prior_ptr = new list<Node*>();
        list<Node*>& prior = *prior_ptr;
        list<uint32_t>* ko_list_ptr = new list<uint32_t>();
        list<uint32_t>& ko_list = *ko_list_ptr;
    	
        uint32_t d_max = 2500;
        node_map nm;
		nm = parse_graph(g);
        fres << "curgraph\t" << g << "\t" << nm.size() << endl;
    	uint32_t node_cnt = nm.size();
    	double eps0, eps1;
    	double* real_scores = new double[nm.size()];
    	double* power_scores = new double[nm.size()];
    	double* fp_dyn_scores = new double[nm.size()];
    	double* scratch_scores = new double[nm.size()];

        cout << "Total nodes count is:" << nm.size() << endl;
		cout << "Generating random prior, size=" << prior_size << endl;
		get_random_prior(prior_size, nm, nm.size(), prior);
		cout << "Generating random ko list, size=" << ko_size << endl;
		get_ko_list(ko_size, prior, nm.size(), ko_list);
		cout << "Generating results for graph " << g << endl;
        // For now, propgraph computes propgation over the original each time.
        // TODO: implement remove/add node for PropGraph also.
        PropGraph gp(nm, true, -1); 
		DyGraph dg(nm, true);
		DyGraph dg2(nm, true);
		dg.set_alpha(alpha);
		dg2.set_alpha(alpha);
		//eps0 = sqrt(eps / (node_cnt*d_max)) / node_cnt;
        eps0 = 1.0 / (node_cnt * d_max);
        eps1 = (1.0 / (d_max * d_max));
		//eps1 = eps / (node_cnt*d_max);
        //eps0 = eps1 = 1.0 / (node_cnt * d_max);
        eps = 0.05;  // minimal for the optimal solution for power-iteration.i
        //eps0 = eps1 = eps = 0.1;
		cout << "Running first propgations..." << endl;
        cout << "Parameters: alpha=" << alpha << "eps=" << eps << endl;
		test_power(gp, eps/(node_cnt), alpha, nm, prior, real_scores, -1, NULL);
		cout << "finished running first power" << endl;
        cout << "Running first fp-scratch..." << endl;
        test_fp_scratch(dg, prior, scratch_scores, eps0, -1, real_scores);
        cout << "Done with first fp-scratch" << endl;
        for (uint32_t k : ko_list) {
            cout << "Knocking out node" << k << "for all algs" << endl;
			// Getting the "real-scores" approximation.
			test_power(gp, eps / node_cnt, alpha, nm, prior, real_scores, k, NULL);
			
            
            // TODO: put here the optimal/minimal epsilon for this error to occur.
			// Getting approximation with power-iteration.
			test_power(gp, eps, alpha, nm, prior, power_scores, k, real_scores);

			// Getting approximation with forward-push from scratch.
			test_fp_scratch(dg2, prior, scratch_scores, eps0, k, real_scores);

			// Getting approximation with forward-push dynamic.
			test_fp_dyn(dg, fp_dyn_scores, eps1, k, real_scores);
		}
        cout << "finished for graph " << g << endl;
        delete prior_ptr;
        delete ko_list_ptr;
        delete[] real_scores;
        delete[] power_scores;
        delete[] fp_dyn_scores;
        delete[] scratch_scores;
	}
}

void do_test_best_tuple_rm(double alpha, double eps, list<Node*>& prior, list<Node*>& l, string g, string res_file="res.txt") {
  /* 
  Find the gene tuples/triples/.. that when we removed them the DEG's ranks were maximally changed (L1 diff of ranks vector before and after scores for example).
  that when we remove the residuals total change was the highest &&
  */
	fres.open(res_file);
	cout << "Generating results for graph " << g << endl;
	node_map nm = parse_graph(g);
	DyGraph dg(nm, true);
	dg.set_prior(prior);
	uint32_t* prior_arr = new uint32_t[nm.size()];
	for (Node* n : prior)
		prior_arr[n->get_node()] = 1;
	dg.set_alpha(alpha);
	dg.set_eps(eps);
	double* scores = new double[nm.size()];
	double max_res_added = 0.0;
	double max_res_deg = 0.0;
	for (uint32_t k = 0; k < nm.size(); k++) {
		if (prior_arr[k]) // the ko node is in the prior set.
			continue; 
		double cur_res = dg.get_res_add(k); // simulate how much residual will be added/removed upon this case.
		if (cur_res < max_res_deg ) {
			continue;
		}
		dg.del_node(k, scores);
		dg.enable_node(k);
	}
}

int main(int argc, char** argv) {
	uint32_t ko_size = 1000;
	uint32_t prior_size = 100;
	double alpha = 0.4;
	double eps = 0.01;
    /*try
    {
        /// Define and parse the program options
        namespace po = boost::program_options;
        po::options_description desc("Options");
        desc.add_options()
        ("help", "Print help messages")
        ("ksize", po::value<uint32_t>(&ko_size), "Number of nodes to knockout")
        ("psize", po::value<uint32_t>(&prior_size), "Prior size")
        ("alpha", po::value<double>(&alpha), "alpha, the reset probability")
        ("eps", po::value<double>(&eps), "epsilon, the maximal L1-error rate.");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw
    
            // --help option
         if (vm.count("help"))
         {
              std::cout << "DyGraph - Tests program" << std::endl
                       << desc << std::endl;
             return SUCCESS;
         }
        
          //po::notify(vm); // throws on error, so do after help in case
         // there are any problems
        }
        catch(po::error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }
    }
    catch(std::exception& e) 
    { 
        std::cerr << "Unhandled Exception reached the top of main: " 
                  << e.what() << ", application will now exit" << std::endl; 
        return ERROR_UNHANDLED_EXCEPTION; 
    } 
    */

	do_test1(alpha, eps, prior_size, ko_size, "res_out.txt");
    return SUCCESS;
}
