#include "dygraph.hpp"
#include "propgraph.hpp"
#include <assert.h>

double DyGraph::sum_node(Node* node)
{
	double s = 0.0;
	list<Edge*> edges = node->get_edges();
	for (auto it = edges.begin(); it != edges.end(); ++it)
	{
		s += (*it)->w();
	}
	return s;
}

void DyGraph::_add_to_hlist(Node* node_obj)
{
	//uint32_t node_id = node_obj->get_node();
	
	// pushing pos-res nodes at the end of the list.
	//if (m_R[node_id] > 0) {
		// TODO: implement
	//}

	if (m_hlist)
	{
		node_obj->set_next(m_hlist);
		m_hlist->set_prev(node_obj);
	}
	node_obj->set_in_hlist(true);
	m_hlist = node_obj;
	node_obj->set_prev(nullptr);
}

void DyGraph::_del_from_hlist(Node* node_obj)
{
	if (node_obj == m_hlist) {
		_pop_hlist();
		return;
	}

	Node* next, *prev;
	next = node_obj->get_next();
	prev = node_obj->get_prev();

	node_obj->set_in_hlist(false);
	node_obj->set_prev(nullptr);
	node_obj->set_next(nullptr);

	if (prev && next)
	{
		prev->set_next(next);
		next->set_prev(prev);
	}
	else if (prev)
	{
		prev->set_next(nullptr);
	}
	else if (next)
	{
		next->set_prev(nullptr);
	}
}

Node* DyGraph::_pop_hlist()
{
	if (m_hlist)
	{
		Node* pnode = m_hlist;
		m_hlist = m_hlist->get_next();
		pnode->set_in_hlist(false);
		pnode->set_next(nullptr);
		pnode->set_prev(nullptr);
		if (m_hlist) {
			m_hlist->set_prev(nullptr);
		}
		return pnode;
	}
	return nullptr;
}

void DyGraph::sum_r() {
	double s = 0.0;
	for (int i = 0; i < m_node_cnt; i++)
		s += m_R[i];
	cout << "DyGraph: m_R sum is:" << s << endl;
}

void DyGraph::set_prior(list<Node*>& prior)
{
	m_prior_size = prior.size();
	for (auto& it = prior.begin(); it != prior.end(); ++it)
	{
		Node* node_obj = *it;
		uint32_t n = node_obj->get_node();
		m_prior[n] = true;
		m_R[n] = 1.0 / m_prior_size;
		if (m_sym_norm)
			m_R[n] *= sqrt(node_sums[n]);
		_add_to_hlist(m_nodes[n]);
	}
	sum_r();
}

const double* DyGraph::propagate(double* scores)
{
	Node* pnode;
	while (m_hlist)
	{
		// pop hlist head.
		pnode = _pop_hlist();

		// Push current hlist head.
		_forward_push(pnode);

		//sum_r();
		//cout << "nt: " << m_nt << endl;
	}

	cout << "nt: " << m_nt << endl;
	sum_r();
	m_nt = 0;
	// Setup and return scores.
	double* out_scores;
	if (scores)
		out_scores = scores;
	else
		out_scores = new double[m_node_cnt];

	for (uint32_t i = 0; i < m_node_cnt; i++)
		if (node_sums[i] > 1e-10) {
			out_scores[i] = m_P[i] * (m_sym_norm ? 1.0 / sqrt(node_sums[i]) : 1);
		}

	/*
	double s = 0.0;
	for (uint32_t i = 0; i < m_node_cnt; i++)
		s += out_scores[i];
	cout << "DyGraph: sum of s:" << s << endl;
	*/
	return out_scores;
	
}

void DyGraph::set_R(double* R) {
	for (uint32_t i = 0; i < m_node_cnt; i++)
		m_R[i] = R[i];
}

void DyGraph::set_P(double* P) {
	for (uint32_t i = 0; i < m_node_cnt; i++)
		m_P[i] = P[i];
}

double* DyGraph::get_R() {
	return m_R;
}

double* DyGraph::get_P() {
	return m_P;
}

int DyGraph::get_length(Node* node_l) {
	uint32_t cnt = 0;
	while (node_l) {
		cnt++;
		node_l = node_l->get_next();
	}
	return cnt;
}

void DyGraph::_forward_push(Node* node)
{
	uint32_t other_node_id;
	uint32_t node_id = node->get_node();
	list<Edge*>& edges = node->get_edges();

	m_P[node_id] += alpha * m_R[node_id];
	//cout << "pushing from " << node_id << ", ttl_sum friends= " << node_sums[node_id] << endl;
	for (auto& it = edges.begin(); it != edges.end(); ++it)
	{
		m_nt++;
		Edge* edge = *it;
		if (!edge->valid() || node_sums[node_id] < 1e-10)
			continue;
		other_node_id = edge->node2();
		//cout << "pushing to " << other_node_id << " w = " << edge->w() << endl;
		m_R[other_node_id] += (1 - alpha) * edge->w() * m_R[node_id] * (1.0 / node_sums[node_id]);
		uint32_t deg = m_nodes[other_node_id]->get_valid_edges(); // TODO: count only valid edges?
		bool above_eps = (m_R[other_node_id] > eps * deg)  || (m_R[other_node_id] < -eps * deg);
		if (above_eps && !m_nodes[other_node_id]->in_hlist())
		{
			_add_to_hlist(m_nodes[other_node_id]);
		}
		else if (!above_eps && m_nodes[other_node_id]->in_hlist()) {
			_del_from_hlist(m_nodes[other_node_id]);
		}
	}
	m_R[node_id] = 0.0;
}

const double* DyGraph::del_node(uint32_t node, double* scores)
{
	Node* rm_node = m_nodes[node];
	list<Edge*> edges = rm_node->get_edges();

	if (m_prior[node])
	{
		cout << "removing prior not fully implemented yet." << endl;
		getchar();
		exit(0);
		m_prior_size--;
		m_prior[node] = false;
		assert(m_prior_size > 0); // TODO: handle the case this is the last in the prior set, go to usual RW.
	}

	// Make all edges to this node invalid state.
	for (auto& it = edges.begin(); it != edges.end(); ++it)
	{
		Edge* edge = *it;
		// 1. invalidate edges.
		edge->set_valid(false);
		edge->rev_edge()->set_valid(false);

		// 2. normalize all its neighbors.
		uint32_t node2 = edge->node2();
		double old_sum = node_sums[node2];
		double new_sum = node_sums[node2] - edge->w();
		double f_uv = old_sum / new_sum;
		node_sums[node2] = new_sum;

		// 3. Fix invariants' R + P (and redo propagate?)
		// R(v) += (1/alpha)*(1 - (1/F_u(v)))P(v) - ((1-alpha) / alpha)*P(u)*W_uv, P(v) /= F_u(v)
		double w_old = edge->w() / node_sums[node];
		m_R[node2] += (1.0 / alpha) * (1 - (1.0 / f_uv)) * m_P[node2];
		m_R[node2] -= ((1.0 - alpha) / alpha) * m_P[node] * w_old;
		m_P[node2] /= f_uv;

		// 3* neighbor node becomes isolated, hence clear its P and R.
		if (new_sum < 1e-10) {
			m_R[node2] = 0.0;
			m_P[node2] = 0.0;
			continue;
		}

		// + If its symmetric, add the diff in the (L^-1) * prior.
		if (m_sym_norm && m_prior[node2] && (new_sum > 1e-10))
		{
			m_R[node2] += ((sqrt(new_sum)) - (sqrt(old_sum))) * (1.0 / m_prior_size);
		}

		// 4: update valid edges count for neighbor.
		Node* n2 = m_nodes[node2];
		n2->set_valid_edges(n2->get_valid_edges() - 1);
		uint32_t deg_node2 = n2->get_valid_edges();

		// 5: push the node to updated nodes list.
		if (abs(m_R[node2]) > eps * deg_node2)
			_add_to_hlist(m_nodes[node2]);
	}

	// clear deleted node.
	m_P[node] = 0.0;
	m_R[node] = 0.0;

	// return propagated scores.
	return propagate(scores);
}


void DyGraph::enable_node(uint32_t node) {
	Node* rm_node = m_nodes[node];
	list<Edge*> edges = rm_node->get_edges();
	// TODO: handle removed prior node, need to re-add it to the prior here somehow.
	double ttl_sum = 0.0;

	for (auto& it = edges.begin(); it != edges.end(); ++it) {
		Edge* edge = *it;
		uint32_t node2 = edge->node2();
		// 1. valid edges again.
		edge->set_valid(true);
		edge->rev_edge()->set_valid(true);

		// 2. add the node_sums accordingly.
		double new_sum = node_sums[node2] + edge->w();
		node_sums[node2] = new_sum;

		// 3: update valid edges count for neighbor.
		Node* n2 = m_nodes[node2];
		n2->set_valid_edges(n2->get_valid_edges() + 1);

		ttl_sum += edge->w();
	}

	node_sums[node] = ttl_sum;
}