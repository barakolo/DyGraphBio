// This is the hpp file.
#include "node.hpp"
using namespace std;

list<Edge*>& Node::get_edges() {
	return m_edges;
}

uint32_t Node::get_node() {
	return m_node;
}

string Node::get_name() {
	return m_node_name;
}

Node* Node::get_next() {
	return m_next;
}

Node* Node::get_prev() {
	return m_prev;
}

void Node::set_next(Node* next) {
	m_next = next;
}

void Node::set_prev(Node* prev) {
	m_prev = prev;
}

Edge* Node::add_edge(uint32_t node1, uint32_t node2, double w) {
	Edge* edge = new Edge(node1, node2, w);
	m_edges.push_back(edge);
	m_valid_edges++;
	return edge;

}

bool Node::in_hlist() {
	return m_in_hlist;
}

void Node::set_in_hlist(bool v) {
	m_in_hlist = v;
}

void Node::set_valid_edges(uint32_t ve) {
	m_valid_edges = ve;
}

uint32_t Node::get_valid_edges() {
	return m_valid_edges;
}