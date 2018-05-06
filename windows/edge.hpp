// This is the hpp file.
#include <iostream>
#include <string>
#include <stdint.h>
#include <list>
using namespace std;

#ifndef EDGE_CLASS
#define EDGE_CLASS 1
class Edge {
	uint32_t m_node1;
	uint32_t m_node2;
	double m_weight;
	bool m_valid;
	Edge* m_rev_edge; // optional: its suitable reverse edge on the opposite direction. used in del_node for eficiency.

public:
	Edge(uint32_t node1 = -1, uint32_t node2 = -1, double weight = 0.0, double valid = true, Edge* m_rev_edge = NULL) :
		m_node1(node1), m_node2(node2), m_weight(weight), m_valid(valid), m_rev_edge(NULL) {};
	double w();
	uint32_t node1();
	uint32_t node2();
	Edge* rev_edge();
	void set_rev_edge(Edge* e);
	bool valid();
	void Edge::set_valid(bool v);
};
#endif // ! EDGE_CLASS

