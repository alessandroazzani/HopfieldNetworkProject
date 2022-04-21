#include<iostream>
#include<vector>
#include<cassert>
#include<cmath>
#include<random>
#include<fstream>
#include<string>
#include<ctime>
#include <stdio.h>
#include <stdlib.h>




class Neuron {
	int clock = 0;
	int state = 0;

public:
	void flowing_time() {
		if (clock != 0) { clock = clock - 1; }
	}

	void set_clock(int time) { clock = time; }

	int get_clock() { return clock; }

	double get_state() {
		return state;
	}

	void set_state(int state_) {
		state = state_;
	}
};

double operator*(std::vector<double> l, std::vector<Neuron> r) {
	assert(l.size() == r.size());
	double result = 0;
	auto itleft = l.begin();
	auto itright = r.begin();
	for (;itleft != l.end();itleft++, itright++) {
		result = result + (*itleft * (*itright).get_state());
	}
	return result;
}

class Graph {

	void generate_adj() {
		std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis_int(0, n_nodes - 1);
		std::uniform_real_distribution<> dis_real(-0.5, 1);

		auto itr = adj.begin();
		auto const endr = adj.end();
		std::vector<double> row(n_nodes);
		int N = 0;
		int r = 0;

		for (;itr != endr; ++itr, ++r) {
			for (int i = 0; i != in_degree; ++i) {
				N = dis_int(gen);
				if (row[N] == 0 && (adj[N])[r] == 0 && N != r) {
					row[N] = dis_real(gen);
				}
				else { --i; }
			}
			*itr = row;
			for (auto it = row.begin(); it != row.end(); it++) {
				*it = 0;
			}

		}
	}

protected:

	int n_nodes = 0;
	int in_degree = 0;
	int time_active = 0;
	int time_passive = 0;
	int retard = 0;
	std::vector<std::vector<double>> memory;
	std::vector<std::vector<double>> adj{ 0 }; //elemento 1 2 è diretto da 2 ad 1(è la trasposta...)
	std::vector<Neuron> state;
	int time = 0;

	auto get_row(int row) {
		assert(row < n_nodes);
		auto it_row = adj.begin();
		for (int i = 0; i != row; ++i) {
			++it_row;
		}
		return it_row;
	}

public:

	Graph(int n_nodes_, int in_degree_, int time_active_, int time_passive_, int retard_, bool is_cluster = false) {
		n_nodes = n_nodes_;
		in_degree = in_degree_;
		time_active = time_active_ - 1;
		time_passive = time_passive_ - 1;
		retard = retard_ - 1;
		assert(in_degree < (n_nodes - 1));
		state.resize(n_nodes);
		if(!is_cluster){
			generate_adj();
		}
		memory.resize(retard);
		for (int i = 0; i != memory.size(); ++i) {
			memory[i].resize(n_nodes);
		}
		adj.resize(n_nodes);
		for (int i = 0; i != n_nodes; ++i) {
			adj[i].resize(n_nodes);
		}
	}

	void bipartite_adj(){
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis_int(0, n_nodes/2 - 1);
		std::uniform_real_distribution<> dis_real(5, 10);

		//Clear adjacency matrix
		auto itr = adj.begin();
		auto const endr = adj.end();
		auto itc = itr->begin();
		auto endc = itr->end();
		for(; itr != endr; ++itr){
			itc = itr->begin();
			endc = itr->end();
			for(; itc != endc; ++itc){
				(*itc) = 0;
			}
		}

		//Create a bipartite adjacency matrix
		itr = adj.begin();
		std::vector<double> row(n_nodes);
		int N = 0;
		int r = 0;
		int half = 1;
		for (;itr != endr; ++itr, ++r) {
			if (r >= n_nodes / 2) {half = 0;}
			for (int i = 0; i != in_degree; ++i) {
				N = dis_int(gen);
				if (row[N + half * n_nodes / 2] == 0 && (adj[N + half * n_nodes / 2])[r] == 0) {
					row[N + half * n_nodes / 2] = dis_real(gen);
				}
				else { --i; }
			}
			*itr = row;
			for (auto it = row.begin(); it != row.end(); it++) {
				*it = 0;
			}

		}
	}

	int Heaviside(double x) {
		if (x > 0) { return 1; }
		else { return 0; }
	};

	void Heaviside(std::vector<double> v) {
		auto it_state = state.begin();
		for (auto itv = v.begin(); itv != v.end(); itv++) {
			if ((*it_state).get_clock() == 0 && (*it_state).get_state() == 0) {
				(*it_state).set_state(Heaviside(*itv));
				(*it_state).set_clock(time_active);
			}
			else if ((*it_state).get_clock() == 0 && (*it_state).get_state() == 1) {
				(*it_state).set_state(0);
				(*it_state).set_clock(time_passive);
			}
			else {
				(*it_state).flowing_time();
			}
			++it_state;
		}
	}

	void next_step() {
		time++;
		std::vector<double> next(n_nodes);
		auto it_next = next.begin();
		auto end = adj.end();
		for (auto it = adj.begin(); it != end; ++it) {
			*it_next = ((*it) * (state));
			it_next++;
		}
		memory.push_back(next);
		Heaviside(memory[0]);
		memory.erase(memory.begin());
	}

	int get_time() { return time; }

	std::vector<int> get_state() {
		std::vector<int> int_state(state.size());
		auto its = int_state.begin();
		for (auto it = state.begin(); it != state.end(); ++it, ++its) {
			*its = (*it).get_state();
		}
		return int_state;
	}

	void random_init() {
		//srand((unsigned) std::time(0));
		//std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen((unsigned) std::time(0));
		std::uniform_int_distribution<> dis_binary(0, 1);

		for (int i = 0; i != n_nodes; ++i) {
			int random = dis_binary(gen);
			if (random == 1) {
				state[i].set_state(1);
				state[i].set_clock(time_active);
			}
			else { state[i].set_state(0); }
		}
	}

	void all_firing() {
		for (int i = 0; i != n_nodes; ++i) {
			state[i].set_state(1);
		}
	}

	void print_state() {
		for (int i = 0; i != n_nodes; ++i) {
			std::cout << state[i].get_state() << " ";
		}
		std::cout << get_sync() << '\n';
	}

	//basic synchronization
	double get_sync() {
		double sync = 0;
		int s;
		for (auto it = state.begin(); it != state.end(); ++it) {
			s = (*it).get_state();
			switch (s)
			{
			case 0:
				sync += -1;
				break;
			case 1:
				sync += 1;
				break;
			default:
				break;
			}
		}
		sync = sync / n_nodes;
		return sync;
	}

	//Synchronization coefficient for bipartite graphs
	double bipartite_sync() {
		double sync_2 = 0;
		double sync_1 = 0;
		int s;
		for (int i = 0; i != n_nodes; ++i) {
			s = state[i].get_state();
			if(i == n_nodes / 2) {
				sync_1 = sync_2;
				sync_2 = 0;
			}
			switch (s)
			{
			case 0:
				sync_2 += -1;
				break;
			case 1:
				sync_2 += 1;
				break;
			default:
				break;
			}
		}
		return (std::abs(sync_1) + std::abs(sync_2)) / n_nodes;
	}

	void print_adj() {
		auto itr = adj.begin();
		auto const endr = adj.end();
		auto endc = itr->end();
		auto itc = itr->begin();

		for (int i = 1;itr != endr;++itr, i++) {
			itc = itr->begin();
			endc = itr->end();
			std::cout << '\n';
			for (;itc != endc;++itc) {
				std::cout << *itc << "       ";
			}
		}
		std::cout << '\n';
	}

	void print_adj_txt() {
		auto itr = adj.begin();
		auto const endr = adj.end();
		auto endc = itr->end();
		auto itc = itr->begin();
		std::ofstream SaveFile("Matrix.txt");
		SaveFile << "Adjacency Matrix" << '\n';
		for (;itr != endr;++itr) {
			itc = itr->begin();
			endc = itr->end();
			for (;itc != endc;++itc) {
				SaveFile << *itc << '\n';  //this function prints the adj matrix in a txt file
			}
		}
		SaveFile << '\n';
		SaveFile.close();
	}

	void write_adj() {
		auto itr = adj.begin();
		auto const endr = adj.end();
		auto endc = itr->end();
		auto itc = itr->begin();
		std::string val;

		std::ifstream inFile;
		inFile.open("Matrix.txt");
		while (getline(inFile, val)) {
			for (;itr != endr; ++itr) {
				itc = itr->begin();
				endc = itr->end();
				for (;itc != endc;++itc) {
					inFile >> val;
					(*itc) = stod(val);
				}
			}
		}
		inFile.close();
	}

};

class Cluster : public Graph{
	int num_clusters = 0;
	int nodes_in_cluster = 0;
	int inter_connections = 0;
	
	void cluster_adj(){
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis_int(0, nodes_in_cluster - 1);
		std::uniform_real_distribution<> dis_real(5, 10);	//intracluster
		std::uniform_real_distribution<> dis_inhibitor(-150, -140); //intercluster


		//Create connections inside the clusters
		auto itr = adj.begin();
		auto endr = adj.end();
		auto itc = itr->begin();
		auto endc = itr->end();
		std::vector<double> row(n_nodes);
		int N = 0;
		int r = 0;
		int j = 0;
		for(int index_cluster = 0; index_cluster != num_clusters; ++index_cluster){
			for (int j = 0 ;j != nodes_in_cluster; ++j, ++itr, ++r) {
				for (int i = 0; i != in_degree; ++i) {
					N = dis_int(gen);
					if (row[N + index_cluster * nodes_in_cluster] == 0 && (adj[N + index_cluster * nodes_in_cluster])[r] == 0 && (N + index_cluster * nodes_in_cluster) != r) {
						row[N + index_cluster * nodes_in_cluster] = dis_real(gen);
					}
					else { --i; }
				}
				*itr = row;
				for (auto it = row.begin(); it != row.end(); it++) {
					*it = 0;
				}
			}
			//Create inhibitory connections between the clusters
			for(int end_cluster = 0; end_cluster != num_clusters; ++end_cluster){
				if(end_cluster == index_cluster){
					//Do nothing, skip.
				} else {
					for (int i = 0; i != inter_connections; ++i) {
						N = dis_int(gen);
						int C = dis_int(gen);
						if ((adj[C + index_cluster * nodes_in_cluster])[N + end_cluster * nodes_in_cluster] == 0 && (adj[N + end_cluster * nodes_in_cluster])[C + index_cluster * nodes_in_cluster] == 0 ) {
							(adj[C + index_cluster * nodes_in_cluster])[N + end_cluster * nodes_in_cluster] = dis_inhibitor(gen);
						}
						else { --i; }
					}
				}
			}
		}

	}

	public:
	Cluster(int num_clusters_, int nodes_in_cluster_, int inter_connections_, int in_degree_, int time_active_, int time_passive_, int retard_) 
	: Graph(num_clusters_ * nodes_in_cluster_, in_degree_, time_active_, time_passive_, retard_, true) {
		inter_connections = inter_connections_;
		num_clusters = num_clusters_;
		nodes_in_cluster = nodes_in_cluster_;
		cluster_adj();
	}
};


int main() {
	//Parameters used for bipartite graph simulation
	/*
	int time_active = 2;
	int time_passive = 1;
	int retard = 1;
	int number_neurons = 200;
	int in_degree = 40;
	/**/
	
	int time_active = 2;
	int time_passive = 3;
	int retard = 4;
	int number_neurons = 6;
	int in_degree = 1;
	//Graph G(number_neurons, in_degree, time_active, time_passive, retard);

	int num_cluster = 2;
	int nodes_in_cluster = 5;
	int intercluster = 1;
	Cluster G(num_cluster, nodes_in_cluster, intercluster, in_degree, time_active, time_passive, retard);
	//G.random_init();
	G.print_adj();
	//G.print_adj_txt();
	//G.write_adj();
	int steps = 200;
	for (int i = 0; i != steps; ++i) {
		//G.print_state();
		//G.next_step();
		//std::cout << G.bipartite_sync() << '\n';
	}

}