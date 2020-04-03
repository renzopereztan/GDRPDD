// general preamble
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <bitset>
#include <map>
#include <chrono>
// TdZdd-specific preamble
#include "DdSpecOp.hpp"
#include "DdStructure.hpp"
#include "util/Graph.hpp"
#include "util/IntSubset.hpp"
#include "spec/DegreeConstraint.hpp"
#include "spec/FrontierBasedSearch.hpp"

#define MAX_EDGES 108
#define MAX_LINES 10

using namespace std;

// number of edges and lines
int n=0, L=0;
// line and weight assignments
vector<int> l, d;

// stores the data to represent the results of the ZDD
struct DistData {
	int val, s, t;
	bitset<MAX_EDGES> mask;
	
	DistData(){
		val = 0;
		s = -1;
		t = -1;
	}
	
	DistData(int v){
		val = v;
		s = -1;
		t = -1;
	}
	
	void print(){
		cout << "Path between (" << s << ", " << t << ")" << endl;
		cout << "Distance = " << val << endl;
		cout << "Yielding Path = (";
		for(int i=0, j=0; i<n; i++){
			if(mask[i]){
				if(j != 0) cout << " ";
				cout << n-i-1;
				j++;
			}
		}
		std::cout << ")" << std::endl;
	}
};

class MaxDist: public tdzdd::DdEval<MaxDist,DistData> {
public:
	MaxDist(){}
	void evalTerminal(DistData& data, bool one) const {
		// initialize
		data.val = one ? 0 : INT_MIN;
	}
	void evalNode(DistData& data, int level, tdzdd::DdValues<DistData,2> const& values) const {
		// principal declaration
		const DistData& data0 = values.get(0);
		const DistData& data1 = values.get(1);
		// general logic for finding the maximum
		if(data0.val >= data1.val + d[n-level]){
			data.val = data0.val;
			data.mask = data0.mask;
		}else{
			data.val = data1.val + d[n-level];
			data.mask = data1.mask;
			data.mask[level-1] = 1;
		}
	}
};

class MinDist: public tdzdd::DdEval<MinDist,DistData> {
public:
	MinDist(){}
	void evalTerminal(DistData& data, bool one) const {
		// initialize
		data.val = one ? 0 : 100000;
	}
	void evalNode(DistData& data, int level, tdzdd::DdValues<DistData,2> const& values) const {
		// principal declaration
		const DistData& data0 = values.get(0);
		const DistData& data1 = values.get(1);
		// general logic for finding the minimum
		if(data0.val <= data1.val + d[n-level]){
			data.val = data0.val;
			data.mask = data0.mask;
		}else{
			data.val = data1.val + d[n-level];
			data.mask = data1.mask;
			data.mask[level-1] = 1;
		}
	}
};

// setting up the constraint of using each line at least once
class Lines: public tdzdd::DdSpec<Lines, int, 2> {
public:
  	Lines(){}
  	int getRoot(int& state) const{
    	state = 0;
    	return n;
  	}
  	int getChild(int& state, int level, int value) const{
		if(value == 1) state |= (1<<l[n-level]);
		level--;
		
		if(level == 0){
			if(state == ((1<<L)-1)<<1){
				return -1;
			}else{
				return 0;
			}
		}
    	return level;
  	}
};

int main(int argc, char *argv[]) {
	auto start_time = chrono::high_resolution_clock::now();
	
	if (argc <= 1) {
        cout << "usage: " << argv[0] << " <graph_filename>" << endl;
        return 0;
    }
	// to display computation time
	tdzdd::MessageHandler::showMessages();
	
    tdzdd::Graph graph;
	
	// read graph file
	string line;
	ifstream file;
	file.open(string(argv[1]));
	
	while(getline(file, line)){
		stringstream ss(line);
		string ta, tb, td, tl;
		getline(ss, ta, ' ');
		getline(ss, tb, ' ');
		getline(ss, td, ' ');
		getline(ss, tl, ' ');
		
		n++;
		graph.addEdge(ta, tb);
		d.push_back(stoi(td));
		L = max(L, stoi(tl));
		l.push_back(stoi(tl));
	}

	file.close();
	graph.update();

    // the diagram representing all paths
    tdzdd::FrontierBasedSearch fbs(graph, 1, false, false);
	
	tdzdd::IntRange zeroOrTwo(0, 2, 2);
	tdzdd::IntRange justOne(1, 1);
	// setting up degree constraint
	tdzdd::DegreeConstraint dc(graph);

	map<int, int> vertex_mapping;
	for (int v = 1; v <= graph.vertexSize(); ++v) {
		vertex_mapping[stoi(graph.vertexName(v))] = v;
		dc.setConstraint(v, &zeroOrTwo);
	}
	
	vector<DistData> max_data;
	vector<DistData> min_data;
	max_data.push_back(DistData(INT_MIN));
	min_data.push_back(DistData(INT_MAX));
	
	for(int s = 1; s <= graph.vertexSize(); ++s){
		for(int t = s+1; t <= graph.vertexSize(); ++t){
			cout << "Computing solution for (" << s << ", " << t << ")" << endl;
			
			dc.setConstraint(vertex_mapping[s], &justOne);
			dc.setConstraint(vertex_mapping[t], &justOne);

			// the ZDD for all paths in general
			tdzdd::DdStructure<2> dd1(tdzdd::zddIntersection(dc, fbs));

			// the ZDD for all feasible paths
			tdzdd::DdStructure<2> dd(tdzdd::zddIntersection(dd1, Lines()));

			cout << "Number of ZDD nodes = " << dd.size() << endl;
			cout << "Number of elements = " << dd.evaluate(tdzdd::ZddCardinality<>()) << endl;
			// maximum distance output
			MaxDist dist1 = MaxDist();
			DistData ans1 = dd.evaluate(dist1);
			cout << "Maximum Distance = " << ans1.val << endl;
			
			if(ans1.val > max_data[0].val) max_data.clear();
			if(max_data.size() == 0 || ans1.val == max_data[0].val){
				for(tdzdd::DdStructure<2>::const_iterator it1 = dd.begin(); it1 != dd.end(); ++it1){
					DistData data = DistData();
					for(set<int>::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2){
						data.val += d[n-(*it2)];
						data.mask[(*it2)-1] = 1;
					}
					if(data.val == ans1.val){
						data.s = s;
						data.t = t;
						max_data.push_back(data);
					}
				}
			}
			
			// minimum distance output
			MinDist dist2 = MinDist();
			DistData ans2 = dd.evaluate(dist2);
			cout << "Minimum Distance = " << ans2.val << endl;
			
			if(ans2.val < min_data[0].val) min_data.clear();
			if(min_data.size() == 0 || ans2.val == min_data[0].val){
				for(tdzdd::DdStructure<2>::const_iterator it1 = dd.begin(); it1 != dd.end(); ++it1){
					DistData data = DistData();
					for(set<int>::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2){
						data.val += d[n-(*it2)];
						if(data.val > ans2.val) break;
						data.mask[(*it2)-1] = 1;
					}
					if(data.val == ans2.val){
						data.s = s;
						data.t = t;
						min_data.push_back(data);
					}
				}
			}
			
			dc.setConstraint(vertex_mapping[s], &zeroOrTwo);
			dc.setConstraint(vertex_mapping[t], &zeroOrTwo);
		}
	}
	
	cout << endl << endl;
	cout << "Maximum distance paths found:" << endl;
	for(int i=0; i<max_data.size(); i++) max_data[i].print();
	
	cout << "Minimum distance paths found:" << endl;
	for(int i=0; i<min_data.size(); i++) min_data[i].print();

	auto end_time = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count(); 
  
    time_taken *= 1e-9;
    cout << "Total time taken: " << time_taken << setprecision(9) << "s" << endl;
	
    return 0;
}