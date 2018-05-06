#include "propgraph.hpp"
#include <assert.h>

void PropGraph::set_prior(list<Node*>& prior)
{
	m_prior_size = prior.size();
	for (auto& it = prior.begin(); it != prior.end(); ++it)
	{
		Node* node_obj = *it;
		uint32_t n = node_obj->get_node();
		m_prior[n] = 1.0 / m_prior_size;
	}
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

vec PropGraph::propagate(double* scores)
{
	// Multiply matrix "-log(alpha * eps)" times.
	// return scores.
	int times = (int)std::round((log(eps) / log(1.0 - alpha)));
	if (times <= 1) {
		times = 2;
	}
	vec res = alpha * m_prior;
	vec res2 = boost::numeric::ublas::prod(m, (1 - alpha) * m_prior) + (alpha * m_prior);

	uint32_t ctr = 1;
	while (ctr < times) {
		//cout << "prop graph cur times" << ctr << endl;
		//cout << "reses diff:" << vec_sub_sum(res, res2) << endl;
		//cout << "sum of res2";
		//PropGraph::sum_v(res2);
		if (ctr % 2 == 0) {
			res2 = boost::numeric::ublas::prod(m, (1 - alpha) * res) + (alpha * m_prior);
		}
		else {
			res = boost::numeric::ublas::prod(m, (1 - alpha) * res2) + (alpha * m_prior);
		}
		ctr++;
	}

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

const double* PropGraph::del_node(uint32_t node, double* scores)
{
	/*list<Edge*> edges = m_nodes[node]->get_edges();
	for (auto& it = edges.begin();  */
	return nullptr;
}
