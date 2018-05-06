#include "edge.hpp"
using namespace std;

double Edge::w() { return m_weight; }

uint32_t Edge::node1() { return m_node1; }

uint32_t Edge::node2() { return m_node2; }

bool Edge::valid() { return m_valid; }

// Note: if we invalidate edge, for DyGraph, we should call set_valid_edges, to update valid edges count.
void Edge::set_valid(bool v) { m_valid = v; }

Edge* Edge::rev_edge() { return m_rev_edge;  }

void Edge::set_rev_edge(Edge* e) { m_rev_edge = e;  };