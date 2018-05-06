#include "node.hpp"
#include "dygraph.hpp"
#include <stdint.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#ifndef PROPGRAPH_CLASS
#define PROPGRAPH_CLASS 1

typedef boost::numeric::ublas::compressed_matrix<double> mat;
typedef boost::numeric::ublas::vector<double> vec;

class PropGraph {
	double eps = 1e-4;
	double alpha = 0.4;
	uint32_t m_node_cnt;
	node_map m_nodes;
	double* node_sums;
	vec m_prior;
	uint32_t m_prior_size = 0;
	bool m_sym_norm = false; // Whether it is symmetric normalization or regular stochastic.
	mat m;
	uint32_t rm_node = -1; //TODO make this object changeable instead of reconstructing the obj?

public:
	PropGraph(node_map nodes, bool sym_norm = false, uint32_t rm_node = -1) : m_nodes(nodes), m_sym_norm(sym_norm), m_prior(nodes.size(), 0), m(nodes.size(), nodes.size()), rm_node(rm_node) {
		m_node_cnt = nodes.size();
		m_prior_size = 0;
		for (uint32_t i = 0; i < m_node_cnt; i++) {
			m_prior[i] = 0.0;
		}

		if (rm_node != -1) {
			cout << "removing node" << rm_node << endl;
			m_nodes.erase(rm_node);
		}

		node_sums = new double[m_node_cnt]();
		for (auto& it = m_nodes.begin(); it != m_nodes.end(); it++) {
			double s = 0.0;
			Node* node = it->second;
			list<Edge*>& edges = node->get_edges();
			for (auto& it2 = edges.begin(); it2 != edges.end(); it2++) {
				Edge* e = *it2;
				if (e->node1() == rm_node || e->node2() == rm_node) {
					e->set_valid(false);
					e->rev_edge()->set_valid(false);
					continue;
				}
				s += (*it2)->w();
			}
			node_sums[node->get_node()] = s;
		}
		//cout << "summing nodes done." << endl;

		uint32_t cnt_nodes = 0;
		for (auto& it = m_nodes.begin(); it != m_nodes.end(); it++) {
			Node* node = it->second;
			//cout << "node id" << cnt_nodes << " out of " << m_node_cnt << endl;
			cnt_nodes++;
			list<Edge*>& edges = node->get_edges();
			for (auto& it2 = edges.begin(); it2 != edges.end(); it2++) {
				Edge* edge = *it2;
				uint32_t n1 = edge->node1();
				uint32_t n2 = edge->node2();
				if (n1 == rm_node || n2 == rm_node || !edge->valid() || (node_sums[n1] < 1e-10) || (node_sums[n2] < 1e-10))
					continue;
				if (m_sym_norm) {
					m(n1, n2) = (edge->w() / sqrt(node_sums[n1] * node_sums[n2]));
				}
				else {
					m(n1, n2) = (edge->w() / node_sums[n2]);
				}
			}
		}

		//cout << "initing matrix done." << endl;
	};

	void PropGraph::print_vec(const vec& v);
	void set_prior(list<Node*>& prior);
	vec propagate(double* scores);
	const double* del_node(uint32_t node, double* scores);
	static double vec_sub_sum(const vec& a, const vec& b);
	static void sum_v(const vec& v);
};


#endif 