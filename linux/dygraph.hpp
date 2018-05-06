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
	double* m_R_bak = NULL;
	double* m_P = NULL;
	double* m_P_bak = NULL;
	Node* m_hlist = NULL; // Nodes with high residual to maintain.
	bool m_sym_norm = false; // Whether it is symmetric normalization or regular stochastic.

public:
	DyGraph(const node_map& nodes, bool sym_norm = false) : m_nodes(nodes), m_sym_norm(sym_norm) {
		m_node_cnt = nodes.size();
		m_prior = new bool[m_node_cnt]();
		m_prior_size = 0;
		m_R = new double[m_node_cnt]();
        m_R_bak = new double[m_node_cnt]();
		m_P = new double[m_node_cnt]();
        m_P_bak = new double[m_node_cnt]();
		node_sums = new double[m_node_cnt]();
		m_hlist = NULL;

		for (const auto& it : nodes) {
			Node* node = it.second;
			uint32_t node_id = node->get_node();
			node_sums[node_id] = sum_node(node);
			m_R[node_id] = 0;
			m_P[node_id] = 0;
			m_prior[node_id] = false;
		}
	}

    ~DyGraph() {
        delete[] m_prior;
        delete[] m_R;
        delete[] m_R_bak;
        delete[] m_P;
        delete[] m_P_bak;
        delete[] node_sums;
    }

	double sum_node(Node* n);
	void set_prior(const list<Node*>& prior);
    void init_prior(const list<Node*>& prior);
	const double* propagate(double* scores);
	void _forward_push(Node* node);
	void _add_to_hlist(Node* node);
	void _del_from_hlist(Node* node);
	Node* _pop_hlist();
	const double* del_node(uint32_t node, double* scores);
	void enable_node(uint32_t node);
    uint32_t disable_node(uint32_t node);
	void backup_state(); // backup and restore dygraph internal scores computed state.
	void restore_state();
	void set_R(const double* R);
	double* get_R();
	double* get_P();
	void set_P(const double* P);
	int get_length(Node* node_l);
	void sum_r();
    uint32_t get_nt();
    double get_res_add(uint32_t rm_node);
    void set_alpha(double alpha);
    double get_alpha();
    void set_eps(double eps);
    double get_eps();
    uint32_t get_node_cnt();
};

#endif
