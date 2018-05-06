#include "dygraph.hpp"
#include "propgraph.hpp"
#include <assert.h>

double DyGraph::sum_node(Node* node)
{
	double s = 0.0;
	list<Edge*> edges = node->get_edges();
	for (const auto& it : edges)
	{
		s += it->w();
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

void DyGraph::set_prior(const list<Node*>& prior)
{
	m_prior_size = prior.size();
	for (const auto& it : prior)
	{
		Node* node_obj = it;
        //cout << "setting prior node also as" << node_obj->get_node() << endl;
		uint32_t n = node_obj->get_node();
		m_prior[n] = true;
		m_R[n] = 1.0 / m_prior_size;
		if (m_sym_norm)
			m_R[n] *= sqrt(node_sums[n]);
		_add_to_hlist(m_nodes[n]);
	}
	sum_r();
}

void DyGraph::init_prior(const list<Node*>& prior) {
    while (m_hlist) {
        _pop_hlist(); // clear hlist.
    }

    for (uint i=0; i<m_node_cnt; i++) {
        m_R[i] = 0.0;
        m_P[i] = 0.0;
    }

    set_prior(prior);
}

const double* DyGraph::propagate(double* scores)
{
	m_nt = 0;
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

	//cout << "nt: " << m_nt << endl;
	//sum_r();
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

void DyGraph::set_R(const double* R) {
	for (uint32_t i = 0; i < m_node_cnt; i++)
		m_R[i] = R[i];
}

void DyGraph::set_P(const double* P) {
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
	for (const auto& it : edges) {
        m_nt++;
		Edge* edge = it;
		if (!edge->valid() || node_sums[node_id] < 1e-10)
			continue;
		other_node_id = edge->node2();
		//cout << "pushing to " << other_node_id << " w = " << edge->w() << endl;
        Node* other_node = m_nodes[other_node_id];
		m_R[other_node_id] += (1 - alpha) * edge->w() * m_R[node_id] * (1.0 / node_sums[node_id]);
		uint32_t deg = other_node->get_valid_edges(); 
		bool above_eps = (m_R[other_node_id] > eps * deg)  || (m_R[other_node_id] < -eps * deg);
		if (above_eps && !other_node->in_hlist()) {
			_add_to_hlist(other_node);
		}
		else if (!above_eps && other_node->in_hlist()) {
			_del_from_hlist(other_node);
		}
	}
	m_R[node_id] = 0.0;
}

const double* DyGraph::del_node(uint32_t node, double* scores)
{
	Node* rm_node = m_nodes[node];
	list<Edge*>& edges = rm_node->get_edges();

	if (m_prior[node])
	{
		cout << "removing prior not fully implemented yet." << endl;
		getchar();
		exit(0);
		m_prior_size--;
		m_prior[node] = false;
		assert(m_prior_size > 0); // TODO: handle the case this is the last in the prior set, go to usual RW.
	}

    // Start counting total "node-touches".
	uint32_t start_nt = 0;
    start_nt += rm_node->get_valid_edges();

    // Make all edges to this node invalid state.
	for (const auto& it : edges)
	{
		Edge* edge = it;
		// 1. invalidate edges.
		edge->set_valid(false);
		edge->rev_edge()->set_valid(false);


        // 2. neighbor node becomes isolated, hence clear its P and R.
		uint32_t node2 = edge->node2();
		if (m_nodes[node2]->get_valid_edges() == 1) {
            m_nodes[node2]->set_valid_edges(0); // become isolated now.
			node_sums[node2] = 0.0;
            m_R[node2] = 0.0;
			m_P[node2] = m_prior[node2] / (double)m_prior_size;
			continue;
		}

		// 3.a. normalize all its neighbors.
        double old_sum = node_sums[node2];
		double new_sum = node_sums[node2] - edge->w();
		double f_uv = old_sum / new_sum;
		node_sums[node2] = new_sum;
		
		// 3.b. Fix invariants' R + P (and redo propagate?)
		// R(v) += (1/alpha)*(1 - (1/F_u(v)))P(v) - ((1-alpha) / alpha)*P(u)*W_uv, P(negihbor) /= F_u(v)
		double w_old = edge->w() / node_sums[node];
		m_R[node2] += (1.0 / alpha) * (1 - (1.0 / f_uv)) * m_P[node2];
		m_R[node2] -= ((1.0 - alpha) / alpha) * m_P[node] * w_old;
		m_P[node2] /= f_uv;

		// + If its symmetric, add the diff in the (L^-1) * prior.
		if (m_sym_norm && m_prior[node2] && (new_sum > 1e-10))
		{
			m_R[node2] += ((sqrt(new_sum)) - (sqrt(old_sum))) * (1.0 / m_prior_size);
		}

		// 4: update valid edges count for neighbor.
		Node* n2 = m_nodes[node2];
		uint32_t deg_node2 = n2->get_valid_edges() - 1;
		n2->set_valid_edges(deg_node2);

		// 5: push the node to updated nodes list.
		if (abs(m_R[node2]) > eps * deg_node2)
			_add_to_hlist(m_nodes[node2]);
	}

	// clear deleted node.
	m_P[node] = m_prior[node] / (double)m_prior_size;
	m_R[node] = 0.0;

	// return propagated scores.
	const double* res_scores = propagate(scores);
    m_nt += start_nt;
    return res_scores;
}

/*
** Restoring graph structure as if the node "node" were'nt removed. this doesn't restore R,P.
** It restores normalization of weights, valid edges count, and "node"'s edges + weights as well.
** R,P should be restored as well, use the backup_state/restore_state for more guidance about it.
** This doesn't calculate propagation again, so one should apply porpagate if you want the scores after the restore.
*/
void DyGraph::enable_node(uint32_t node) {
	Node* rm_node = m_nodes[node];
	list<Edge*> edges = rm_node->get_edges();
	// TODO: handle removed prior node, need to re-add it to the prior here somehow.
	double ttl_sum = 0.0;

	for (const auto& it : edges) {
		Edge* edge = it;
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

uint32_t DyGraph::disable_node(uint32_t node)
{
	Node* rm_node = m_nodes[node];
	list<Edge*>& edges = rm_node->get_edges();

	if (m_prior[node])
	{
		cout << "removing prior not fully implemented yet." << endl;
		getchar();
		exit(0);
		m_prior_size--;
		m_prior[node] = false;
		assert(m_prior_size > 0); // TODO: handle the case this is the last in the prior set, go to usual RW.
	}

    // Start counting total "node-touches".
	uint32_t start_nt = 0;
    start_nt += rm_node->get_valid_edges();

    // Make all edges to this node invalid state.
	for (const auto& it : edges)
	{
		Edge* edge = it;
		// 1. invalidate edges.
		edge->set_valid(false);
		edge->rev_edge()->set_valid(false);


        // 2. neighbor node becomes isolated, hence clear its P and R.
		uint32_t node2 = edge->node2();
		if (m_nodes[node2]->get_valid_edges() == 1) {
            m_nodes[node2]->set_valid_edges(0); // become isolated now.
			node_sums[node2] = 0.0;
            m_R[node2] = 0.0;
			m_P[node2] = m_prior[node2] / (double)m_prior_size;
			continue;
		}

		// 3.a. normalize all its neighbors.
        double old_sum = node_sums[node2];
		double new_sum = node_sums[node2] - edge->w();
		double f_uv = old_sum / new_sum;
		node_sums[node2] = new_sum;
		
        // + If its symmetric, add the diff in the (L^-1) * prior.
		if (m_sym_norm && m_prior[node2] && (new_sum > 1e-10))
		{
			m_R[node2] += ((sqrt(new_sum)) - (sqrt(old_sum))) * (1.0 / m_prior_size);
		}

		// 4: update valid edges count for neighbor.
		Node* n2 = m_nodes[node2];
		uint32_t deg_node2 = n2->get_valid_edges() - 1;
		n2->set_valid_edges(deg_node2);
    }

	// clear deleted node.
	m_P[node] = m_prior[node] / (double)m_prior_size;
	m_R[node] = 0.0;

	// return removed node degree.
    return edges.size();
}


void DyGraph::set_alpha(double alpha) {
	this->alpha = alpha;
}

double DyGraph::get_alpha() { return this->alpha; }

void DyGraph::set_eps(double eps) {
	this->eps = eps;
}

double DyGraph::get_eps() { return this->eps; }
uint32_t DyGraph::get_nt() { return m_nt; }
uint32_t DyGraph::get_node_cnt () { return m_node_cnt; }

void DyGraph::backup_state() {
	for (uint32_t i=0; i < m_node_cnt; i++) {
			m_R_bak[i] = m_R[i];
			m_P_bak[i] = m_P[i];
	}
}

void DyGraph::restore_state() {
	for (uint32_t i=0; i < m_node_cnt; i++) {
			m_R[i] = m_R_bak[i];
			m_P[i] = m_P_bak[i];
	}
}


double DyGraph::get_res_add(uint32_t node) {
	Node* rm_node = m_nodes[node];
	list<Edge*> edges = rm_node->get_edges();
	double ttl_node_change = 0.0;

	// Make all edges to this node invalid state.
	for (const auto& it : edges)
	{
		Edge* edge = it;
		// 1. invalidate edges.
		edge->set_valid(false);
		edge->rev_edge()->set_valid(false);

		// 2. normalize all its neighbors.
		uint32_t node2 = edge->node2();
		double old_sum = node_sums[node2];
		double new_sum = node_sums[node2] - edge->w();
		double f_uv = old_sum / new_sum;

		// 3. Fix invariants' R + P (and redo propagate?)
		// R(v) += (1/alpha)*(1 - (1/F_u(v)))P(v) - ((1-alpha) / alpha)*P(u)*W_uv, P(negihbor) /= F_u(v)
		double w_old = edge->w() / node_sums[node];
		ttl_node_change += abs(((1.0 / alpha) * (1 - (1.0 / f_uv)) * m_P[node2]) - (((1.0 - alpha) / alpha) * m_P[node] * w_old));

		// 3* neighbor node becomes isolated, hence clear its P and R.
		if (new_sum < 1e-10) {
			ttl_node_change += m_R[node2];
			continue;
		}

		// + If its symmetric, add the diff in the (L^-1) * prior.
		if (m_sym_norm && m_prior[node2] && (new_sum > 1e-10))
		{
			ttl_node_change += ((sqrt(new_sum)) - (sqrt(old_sum))) * (1.0 / m_prior_size);
		}
	}

	// return the added residuals for this deletion.
	return ttl_node_change;
}
