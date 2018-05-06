#include "node.hpp"
#include "dygraph.hpp"
#include <stdint.h>
#include "Eigen/Dense"

using namespace Eigen;

#ifndef PROPGRAPH_CLASS
#define PROPGRAPH_CLASS 1

typedef MatrixXd mat;
typedef VectorXd vec;

class PropGraph {
	double eps = 1e-1;
	double alpha = 0.4;
	uint32_t m_node_cnt;
	node_map m_nodes;
	double* node_sums;
	vec m_prior;
    vec m_removed; // backup for line that was removed.
    uint32_t m_removed_node;
	uint32_t m_prior_size = 0;
	bool m_sym_norm = false; // Whether it is symmetric normalization or regular stochastic.
	mat m;
	uint32_t rm_node = -1; //TODO make this object changeable instead of reconstructing the obj?
    uint32_t m_nt = 0;
    uint32_t m_edges_cnt = 0;

public:
	PropGraph(const node_map& nodes, bool sym_norm = false, uint32_t rm_node = -1) : m_nodes(nodes),
    m_sym_norm(sym_norm), m_prior(nodes.size()), m_removed(nodes.size()), m(nodes.size(), nodes.size()), rm_node(rm_node) {
        m_node_cnt = m_nodes.size();
		m_prior_size = 0;
		m_nt = 0;
		for (uint32_t i = 0; i < m_node_cnt; i++) {
			m_prior[i] = 0.0;
		}

		if (rm_node != -1) {
			//cout << "removing node" << rm_node << endl;
			m_nodes.erase(rm_node);
		}
		
        node_sums = new double[m_node_cnt]();
		for (const auto& it : m_nodes) {
			double s = 0.0;
			Node* node = it.second;
			list<Edge*>& edges = node->get_edges();
			for (const auto& it2 : edges) {
				Edge* e = it2;
				if (e->node1() == rm_node || e->node2() == rm_node) {
					continue;
				}
				s += it2->w();
			}
			node_sums[node->get_node()] = s;
		}
		//cout << "summing nodes done." << endl;

		uint32_t cnt_nodes = 0;
		for (const auto& it : m_nodes) {
			Node* node = it.second;
			//cout << "node id" << cnt_nodes << " out of " << m_node_cnt << endl;
			cnt_nodes++;
			list<Edge*>& edges = node->get_edges();
			for (const auto& it2 : edges) {
				Edge* edge = it2;
				uint32_t n1 = edge->node1();
				uint32_t n2 = edge->node2();
				if (n1 == rm_node || n2 == rm_node || (node_sums[n1] < 1e-10) || (node_sums[n2] < 1e-10))
					continue;
				if (m_sym_norm) {
					m(n1, n2) = (edge->w() / sqrt(node_sums[n1] * node_sums[n2]));
                    //m.insert_element(n1, n2, (edge->w() / sqrt(node_sums[n1] * node_sums[n2])));
				}
				else {
					m(n1, n2) = (edge->w() / node_sums[n2]);
                    //m.insert_element(n1, n2, edge->w()/node_sums[n2]);
				}
                m_edges_cnt++;
			}
		}

		//cout << "initing matrix done." << endl;
	};

    double vec_sub_sum(const vec& a, const vec& b);
	void print_vec(const vec& v);
	void set_prior(const list<Node*>& prior);
    void clear_prior();
	vec propagate(double* scores);
	static void sum_v(const vec& v);
    uint32_t get_nt();
    void set_eps(double eps_arg);
    double get_eps();
    void set_alpha(double alpha_arg);
    double get_alpha();
    uint32_t get_node_cnt();
    vec del_node(uint32_t n, double* scores); // delete one node and re-propagate.
    void restore();  // restore the graph state, before node removal.
};


#endif 
