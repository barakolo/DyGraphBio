#include "propgraph.hpp"
#include <assert.h>

void PropGraph::set_prior(const list<Node*>& prior)
{
	m_prior_size = prior.size();
	for (const auto& it : prior)
	{
		Node* node_obj = it;
		uint32_t n = node_obj->get_node();
        assert(n < m_nodes.size() && n >= 0);
		m_prior[n] = 1.0 / m_prior_size;
	}
    sum_v(m_prior);
}
void PropGraph::clear_prior() {
    for (uint32_t i=0; i<m_nodes.size(); i++)
        m_prior[i] = 0.0;
}

double PropGraph::vec_sub_sum(const vec& a, const vec& b) {
	assert(a.size() != b.size());
	double diff = 0.0;
	for (uint32_t i = 0; i < a.size(); i++) {
		diff += abs(a[i] - b[i]);
	}
	return diff;
}

void PropGraph::sum_v(const vec& v) {
	double s = 0.0;
	for (int i = 0; i < v.size(); i++)
		s += v[i];
	cout << "PropGraph: scores sum is:" << s << endl;
}

void PropGraph::print_vec(const vec& v) {
	double s = 0.0;
	cout << "Printing scores vec for PropGraph:" << endl;
	for (int i = 0; i < v.size(); i++)
		cout << "g" << i << ": " << v[i] << endl;
	//cout << "PropGraph: scores sum is:" << s << endl;
}

vec PropGraph::del_node(uint32_t n, double* scores) {
    m_removed_node = n;
    Node* node = m_nodes[n];
    list<Edge*>& edges = node->get_edges();
    for (const auto& it : edges ) {
        Edge* edge = it;
        
        // First, remove the edge from the graph.
        m_removed[edge->node2()] = m(edge->node1(), edge->node2()); // backup the edge removed.
        m(edge->node1(), edge->node2()) = 0.0;
        m(edge->node2(), edge->node1()) = 0.0; 

        // Second, renormalize neighbor.
        list<Edge*>& n2_edges = m_nodes[edge->node2()]->get_edges();
        for (const auto& n2_edge : n2_edges) {
            uint32_t cur_node = n2_edge->node1();
            uint32_t neigh_node = n2_edge->node2();
            // If this is the edge to the removed node, continue, we already removed it.
            if (neigh_node == edge->node1()) { 
                continue;
            }
            m(cur_node, neigh_node) *= sqrt(node_sums[cur_node] / (node_sums[cur_node] - edge->w()));
            m(neigh_node, cur_node) *= sqrt(node_sums[cur_node] / (node_sums[cur_node] - edge->w()));
        }
    }
    return propagate(scores);
}

void PropGraph::restore() {
    uint32_t n = m_removed_node;
    Node* node = m_nodes[n];
    list<Edge*>& edges = node->get_edges();
    for (const auto& it : edges ) {
        Edge* edge = it;
        m(edge->node1(), edge->node2()) = m_removed[edge->node2()];
        m(edge->node2(), edge->node1()) = m_removed[edge->node2()];

        // Second, renormalize neighbor.
        list<Edge*>& n2_edges = m_nodes[edge->node2()]->get_edges();
        for (const auto& n2_edge : n2_edges) {
            uint32_t cur_node = n2_edge->node1();
            uint32_t neigh_node = n2_edge->node2();
            // If this is the edge to the removed node, continue, we already removed it.
            if (neigh_node == edge->node1()) { 
                continue;
            }
            m(cur_node, neigh_node) /= sqrt(node_sums[cur_node] / (node_sums[cur_node] - edge->w()));
            m(neigh_node, cur_node) /= sqrt(node_sums[cur_node] / (node_sums[cur_node] - edge->w()));
        }
    }
}

vec PropGraph::propagate(double* scores)
{
	// Multiply matrix "-log(alpha * eps)" times.
	// return scores.
	m_nt = 0;
	int times = (int)std::round((log(eps) / log(1.0 - alpha)));
	if (times <= 1) {
		times = 2;
	}
    vec prior_start(alpha * m_prior);
	vec res(m_prior);
	// old: vec res2 = boost::numeric::ublas::prod(m, (1 - alpha) * m_prior) + (alpha * m_prior);
    // new:
    vec res2(prior_start);
    res2 += ((1-alpha) * m) * m_prior;

	uint32_t ctr = 1;
	while (ctr < times) {
		if (ctr % 2 == 0) {
            // Advancing last vector, res, into res2.
            res2 = prior_start;
            res2 += ((1 - alpha) * m) * res;
			// old: res2 = boost::numeric::ublas::prod(m, (1 - alpha) * res) + (alpha * m_prior);
		}
		else {
            // Advancing last vector, res, into res2.
            res = prior_start;
            res += ((1 - alpha) * m) * res2;
			// old: res = boost::numeric::ublas::prod(m, (1 - alpha) * res2) + (alpha * m_prior);
		}
		ctr++;
	}

	m_nt = ctr * m_edges_cnt;
    //cout << "size of res vector is " << res.size() << "size of res2 vector is " << res2.size() << endl;
	if (ctr % 2) {
		//PropGraph::sum_v(res2);
		for (uint32_t i = 0; i < res2.size(); i++)
			scores[i] = res2[i];
		return res2;
	}
	else {
		//PropGraph::sum_v(res);
		for (uint32_t i = 0; i < res.size(); i++)
			scores[i] = res[i];
		return res;
	}
}

uint32_t PropGraph::get_nt() {
	return m_nt;
}

double PropGraph::get_alpha() {
	return alpha;	
}

void PropGraph::set_alpha(double alpha_arg) {
	alpha = alpha_arg;
}

double PropGraph::get_eps() {
	return eps;
}

void PropGraph::set_eps(double eps_arg) {
	eps = eps_arg;
}

uint32_t PropGraph::get_node_cnt() { return m_nodes.size(); }

