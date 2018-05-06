#include "node.hpp"
#include <stdint.h>

#ifndef DYGRAPH_CLASS
#define DYGRAPH_CLASS 1

class DyGraph {
	double eps = 1e-7;
	double alpha = 0.4;
	uint32_t m_node_cnt = 0;
	node_map m_nodes;
	double* node_sums = NULL;
	bool* m_prior = NULL;
	uint32_t m_prior_size = 0;
	uint32_t m_nt = 0;
	double* m_R = NULL;
	double* m_P = NULL;
	Node* m_hlist = NULL; // Nodes with high residual to maintain.
	bool m_sym_norm = false; // Whether it is symmetric normalization or regular stochastic.


public:
	DyGraph(node_map nodes, bool sym_norm = false) : m_nodes(nodes), m_sym_norm(sym_norm) {
		m_node_cnt = nodes.size();
		m_prior = new bool[m_node_cnt]();
		m_prior_size = 0;
		m_R = new double[m_node_cnt]();
		m_P = new double[m_node_cnt]();
		node_sums = new double[m_node_cnt]();
		m_hlist = NULL;

		for (auto it = nodes.begin(); it != nodes.end(); it++) {
			Node* node = it->second;
			uint32_t node_id = node->get_node();
			node_sums[node_id] = sum_node(node);
			m_R[node_id] = 0;
			m_P[node_id] = 0;
			m_prior[node_id] = false;
		}
	}

	static double sum_node(Node* n);
	void set_prior(list<Node*>& prior);
	const double* propagate(double* scores);
	void _forward_push(Node* node);
	void _add_to_hlist(Node* node);
	void _del_from_hlist(Node* node);
	Node* _pop_hlist();
	const double* del_node(uint32_t node, double* scores);
	void enable_node(uint32_t node);
	void set_R(double* R);
	double* get_R();
	double* get_P();
	void set_P(double* P);
	int get_length(Node* node_l);
	void sum_r();
};

#endif