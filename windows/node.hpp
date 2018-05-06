// This is the hpp file.
using namespace std;
#include <stdint.h>
#include <unordered_map>
#include <map>
#include "edge.hpp"

#ifndef NODE_CLASS
#define NODE_CLASS 1

class Node {
    uint32_t m_node;
	string m_node_name;
	list<Edge*> m_edges;
	uint32_t m_valid_edges;
	Node* m_next;
	Node* m_prev;
	bool m_in_hlist; // the high residual list bool.

public:
	Node(uint32_t node = -1, string node_name = "", Node* next = NULL, Node* prev = NULL, bool in_hlist = false, uint32_t valid_edges=0) :
		m_node(node), m_node_name(node_name), m_edges(), m_next(next), m_prev(prev), m_in_hlist(in_hlist), m_valid_edges(valid_edges) {};

	Node(uint32_t node, string node_name, list<Edge*> edges, Node* next=NULL, Node* prev=NULL)
    {
		m_node = node;
		m_node_name = node_name;
		m_edges = edges;
		m_next = next;
		m_prev = prev;
	};

    list<Edge*>& get_edges();
    uint32_t get_node();
	string get_name();
	Node* get_next();
	Node* get_prev();
	void set_next(Node* next);
	void set_prev(Node* next);
	bool in_hlist();
	void set_in_hlist(bool v);
	void set_valid_edges(uint32_t ve);
	uint32_t get_valid_edges();
	Edge* add_edge(uint32_t node1, uint32_t node2, double w);
};

#endif

typedef unordered_map<uint32_t, Node*> node_map;
typedef unordered_map<uint32_t, Node*>::iterator node_map_iter;
