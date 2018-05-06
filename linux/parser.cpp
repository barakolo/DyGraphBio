#include "parser.hpp"
#include "string.h"
using namespace std;


/* 
Parsing the graph by name.
Returns: list of nodes along with their allocated edges.
*/
node_map parse_graph(string gname) {
	cout << "Start parsing " << gname << endl;
	ifstream gfile(gname);
	uint32_t gind = 0;
	uint32_t edges_cnt = 0;
	string line;
	unordered_map<uint32_t, uint32_t> node_to_ind;
	node_map nm;
	while (getline(gfile, line)){
		// Example line: "ITA7_HUMAN	0	ACHA_HUMAN	1	0.73	mumble"
		istringstream iss(line);
		string g1_name, g2_name;
		uint32_t g1=0, g2=0;
		double w=0;
		if (!(iss >> g1_name >> g1 >> g2_name >> g2 >> w)) {
			cout << "Parsing failed on line " << line << endl;
			continue;
		}

		
		if (w < 0.2 /*|| (g1 == g2)*/) { 
		// || (g1_name.find(g2_name) != string::npos) || (g2_name.find(g1_name) != string::npos)
			continue; // the weight is zero.
		}

		edges_cnt++;
		
		// Setting the nodes' indexes.
		if (node_to_ind.count(g1) > 0) {
			g1 = node_to_ind[g1];
		}
		else {
			node_to_ind[g1] = gind;
			g1 = gind;
			gind++;
		}
		

		if (node_to_ind.count(g2) > 0) {
			g2 = node_to_ind[g2];
		}
		else {
			node_to_ind[g2] = gind;
			g2 = gind;
			gind++;
		}
		
		// Adding nodes.
		Node* n1, *n2;
		if (nm.count(g1) == 0) {
			n1 = new Node(g1, g1_name, NULL);
			nm[g1] = n1;
		}
		else {
			n1 = nm[g1];
		}

		if (nm.count(g2) == 0) {
			n2 = new Node(g2, g2_name, NULL);
			nm[g2] = n2;
		}
		else {
			n2 = nm[g2];
		}

		// Adding the relevant edge.
		Edge* e1 = n1->add_edge(g1, g2, w);
		Edge* e2 = n2->add_edge(g2, g1, w);

		// Set up pointers to the reverse direction edges.
		e1->set_rev_edge(e2);
		e2->set_rev_edge(e1);
	}
	cout << "Done parsing, edges cnt=" << edges_cnt << endl;
	return nm;
}
