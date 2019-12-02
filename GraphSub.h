#include "Node.h"

// get_clock method for time keeping
void get_time_clock(struct timespec * ts) {
	#ifdef __MACH__
		clock_serv_t cclock; 
		mach_timespec_t mts; 
		host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock); 
		clock_get_time(cclock, &mts); 
		mach_port_deallocate(mach_task_self(), cclock); 
		ts->tv_sec = mts.tv_sec; 
		ts->tv_nsec = mts.tv_nsec; 
	#else
	clock_gettime( CLOCK_MONOTONIC, ts); 
	#endif
}

// Main Class Grpaph
class Graph {

public:
	float c_walkLength, c_numWalks; 

	Graph(char *graphfilename, char *featfilename, int dir_control); 

	void readGraph(char *graphfilename); 
	void readFeat(char *featfilename); 

	void readAndRunBFS(char *queryFileName);
	void readAndRunRR(char *queryFileName, int c_numWalk, int numWalks);
	void dynamicQuery(char *dynamicFileName);

	void addEdge(int src, int dst, int label, int dir_control); 
	int diameter(); 
	void initializeRRParams(int diameter); 

	int rr(int src, int dst, int query_type, vector<int> query_labels, int c_numWalk, int c_walkLength); 
	int bbfs(int src, int dst, int query_type, vector<int> query_labels); 

	int edge_bbfs(int src, int dst, int query_type, vector<int> query_labels);
	int edge_rr(int src, int dst, int query_type, vector<int> query_labels);
	
private:

	struct timespec start, finish;
	double initializationTime; 

	int rand_index; 
	vector<int> rand_vector; 

	vector<Node*> nodes; 
	int numNodes; 
	int numEdges; 

	int numQueries;

	bool rrParamsInitialized; 
	int numWalks, numStops, walkLength; 

}; 

Graph::Graph(char *graphfilename, char *featfilename, int dir_control) {

	string line; 
	initializationTime = 0; 
	double elapsed;

	get_time_clock( &start); 

    numNodes = 0; 

    // Read and add labels and nodes from featfile
	ifstream myfile3(featfilename); 
	while (getline(myfile3, line)) {

		if (line[0] == '#') continue; 
		char *cstr = &line[0u]; 

		char *t = strtok(cstr, " "); 
		int u = atoi(t); 

        for (int j = numNodes; j <= u; j ++ ) {
            nodes.push_back(new Node(j)); 
			numNodes ++ ; 
        }

		t = strtok(NULL, " "); 
		int l = atoi(t); 

		nodes[u]->add_label(l); 
	}
	myfile3.close(); 

	// Sort labels, possible improvement-> use set for maintaining labels
	for (int i = 0; i < numNodes; i ++ ) {
		nodes[i]->sort_labels(); 
	}
	
	// Add edges between Nodes, directed or undirected controlled 
	ifstream myfile2(graphfilename); 
	while (getline(myfile2, line)) {

		if (line[0] == '#') continue; 
		char *cstr = &line[0u]; 

		char *t = strtok(cstr, " "); 
		int u = atoi(t); 

		t = strtok(NULL, " "); 
		int v = atoi(t); 

		int l; 
		t = strtok(NULL, " "); 
		if (t == NULL) {
			l = 0; 
		}
		else 
			l = atoi(t); 

		addEdge(u, v, l, dir_control); 
	}
	myfile2.close();

	// Create random number array
	rand_vector.reserve(numNodes+10);
	for (int i = 0; i < numNodes; i ++) {
		rand_vector.push_back(rand() % numNodes); 
	}
	rand_index = 0; 

	// Estimate diameter of graph by sampling 10 random nodes
	int dia = diameter();
	for (int i = 0; i < 4; i ++ ) {
		int d = diameter(); 
		if (d > dia) dia = d; 
	}
	// Use diameter to initialize the parameters
	initializeRRParams(dia); 

	cout << "NumWalks: " << numWalks << ", WalkLength: " << walkLength << endl; 

	get_time_clock(&finish); 
	initializationTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9); 

	numQueries = 0;
	cout << "1. Initialization Time = " << initializationTime << endl; 

}

void Graph::initializeRRParams(int diameter) {
	numStops = (float)(floor(pow((numNodes), 2.0/3) * pow(log(numNodes), 1.0/3))); 
	walkLength = diameter;
	numWalks = (float)numStops / (2 * walkLength); 

	if (numWalks == 0) numWalks = 1; 
	rrParamsInitialized = true; 

}

// Simple bfs to estimate the depth
int Graph::diameter() {
	vector<int> color;
	color.reserve(numNodes);
	for (int i = 0; i < numNodes; i++) color.push_back(0);
	int src = rand()%numNodes;
	color[src] = 1;
	int dia = 1;
	vector<int> queue;
	int u = src;
	while (u >= 0){
		Node *n = nodes[u];
		int numE = n->numFwdEdges[0];
        vector<int> fwdedges = n -> fwd_labelled_edges[0];
		for (int i = 0; i < numE; i++){
			int v =fwdedges[i];
			if (color[v] >= 1) continue;
			color[v] = color[u] + 1;
			queue.push_back(v);
		}
		if (queue.size() == 0) break;
		u = queue[0];
		queue.erase(queue.begin());
	}
	for (int i = 0; i < numNodes; i++)
		if (dia < (color[i]))
			dia = color[i];
	return dia;
}

// Adding edges constrained on whether directional edge or not
void Graph::addEdge(int src, int dst, int l, int dir_control = 1) {
	if ((src >= numNodes) || (dst >= numNodes)) return; 
	nodes[src]->addEdge(dst, true, l); 
	nodes[dst]->addEdge(src, false, l); 
	if (dir_control == 0) {
		nodes[src]->addEdge(dst, false, l); 
		nodes[dst]->addEdge(src, true, l); 
	}
	numEdges ++ ; 
}

void Graph::readAndRunBFS(char *queryFileName) {

	int total_queries = 0; 

	ifstream myfile(queryFileName); 
	string line; 
	outFile1<<"Result Time"<<endl;

	// query parsing
	while (getline(myfile, line)) {

		total_queries ++; 

		char *cstr = &line[0u];
		char *t = strtok(cstr, " "); 

		t = strtok(NULL, " "); 
		int u = atoi(t); 

		t = strtok(NULL, " "); 
		int v = atoi(t); 

		t = strtok(NULL, " "); 
		int query_type = atoi(t); 

		t = strtok(NULL, " "); 
		int num_labels = atoi(t); 

		vector<int> query_labels; 
		while(num_labels--) {
			t = strtok(NULL, " "); 
			query_labels.push_back(atoi(t)); 
		}

		// bbfs calls
		int x1 = 0; 
        get_time_clock( &start); 
        x1 = bbfs(u, v, query_type, query_labels); 
        get_time_clock( &finish); 

		outFile1 << x1 <<" "<< (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9) << endl; 
	}
	myfile.close();
}

void Graph::readAndRunRR(char *queryFileName, int c_numWalk, int c_walkLength) {

	int total_queries = 0; 
	
	ifstream myfile(queryFileName); 
	string line; 
	outFile2<<"------------"<<endl;
	outFile2<<queryFileName<<" nw="<<numWalks*c_numWalk<<" wl="<<walkLength*c_walkLength<<endl;
	outFile2<<"Result Time"<<endl;

	// query parsing
	while (getline(myfile, line)) {

		total_queries ++; 

		char *cstr = &line[0u];
		char *t = strtok(cstr, " "); 
		int u = atoi(t); 

		t = strtok(NULL, " "); 
		int v = atoi(t); 

		t = strtok(NULL, " "); 
		int query_type = atoi(t); 

		t = strtok(NULL, " "); 
		int num_labels = atoi(t); 

		vector<int> query_labels; 
		while(num_labels--) {
			t = strtok(NULL, " "); 
			query_labels.push_back(atoi(t)); 
		}

		int x2 = 0; 
        get_time_clock(&start); 
        x2 = rr(u, v, query_type, query_labels, c_numWalk, c_walkLength); 
        get_time_clock(&finish); 

		outFile2 << x2 << " " << (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9) << endl; 
	}
	myfile.close();
}


int Graph::bbfs(int src, int dst, int query_type, vector<int> query_labels) {

	// Contains whatever Node has been visited to reach in the current walk
	unordered_set < int > soFar; 

	// Queue contains: Node, label, set of soFar nodes
	vector < pair < int, pair < int, unordered_set<int> > > > queue_F; 
	vector < pair < int, pair < int, unordered_set<int> > > > queue_B; 

	// Current Node -> label number -> vector of paths(vector) to node
	unordered_map < int, unordered_map < int, vector < unordered_set < int > > > > fwd; 
	unordered_map < int, unordered_map < int, vector < unordered_set < int > > > > bwd; 

	// Initialize the queues
	queue_F.push_back(make_pair(src, make_pair(-1, soFar))); 
	queue_B.push_back(make_pair(dst, make_pair(query_labels.size(), soFar))); 

	// If any of the queues becomes empty we have not reachable condition
	while(queue_F.size() > 0 && queue_B.size()>0) {

		// Time limit exceeded
		get_time_clock( &finish); 
		if ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9) > 50) return 2; 

		// Forward direction basic initialization: u->Node, l-> label, soFar-> Nodes already reached
		int u = queue_F[0].first; 
		int l = queue_F[0].second.first; 
		soFar = queue_F[0].second.second; 

		// Pop the first element
		queue_F.erase(queue_F.begin()); 

		// Get where the edges from the Node in the given direction
		vector<int> edges = nodes[u]->fwd_labelled_edges[0]; 
		vector<int> labels_no; 

		// Insert into the Nodes visited so far the current Node
		soFar.insert(u); 

		// change label according to query type
		if (query_type == 2) {
			for (int i = 0; i < query_labels.size(); i ++ ) {
				labels_no.push_back(i); 
			}
		}

		else if (query_type == 3) {
			l ++ ; 
			if (l == query_labels.size())
				l = 0; 
			labels_no.push_back(l); 
		}

		else if (query_type == 4) {
			if (l == -1) {
				labels_no.push_back(0); 
			}
			else {
				labels_no.push_back(l); 
				if (l != query_labels.size() - 1) {
					labels_no.push_back(l + 1); 
				}
			}
		}

		// look at all the possible edges out of the node
		for (int i = 0; i<edges.size(); i ++ ) {
			int e = edges[i]; 

			// Condition path has a cycle
			if (soFar.find(e) != soFar.end()) {
				continue; 
			}

			// Search for query label in the labels of node under question
			vector<int> labels = nodes[e]->labels; 
			for (int i = 0; i < labels_no.size(); i ++ ) {
				l = labels_no[i]; 
				if (binary_search(labels.begin(), labels.end(), query_labels[l])) {

					// No meeting with backward BFS
					if (bwd.find(e) == bwd.end()) {
						fwd[e][l].push_back(soFar); 
						queue_F.push_back(make_pair(e, make_pair(l, soFar))); 
					}

					// Meeting with backward BFS
					else {

						// If same label then do a part by part matching to check whether 
						if (bwd[e].find(l) != bwd[e].end()) {

							// Set to check in
							vector < unordered_set < int > > matches = bwd[e][l]; 
							for (int i = 0; i < matches.size(); i ++) {

								//Iterate over the current set and find in the set to check in
								unordered_set < int > ::iterator it; 
								bool flag = true; 
								for (it = soFar.begin(); it != soFar.end(); it ++ ) {
									if (matches[i].find(*it) != matches[i].end()) {
										flag = false; 
										break; 
									}
								}
								if (flag) {
									return 1; 
								}
							}
						}
						
						fwd[e][l].push_back(soFar); 
						queue_F.push_back(make_pair(e, make_pair(l, soFar))); 
					}
				}
			}
		}

		labels_no.clear(); 

		// Backward walk initialization
		int v = queue_B[0].first; 
		l = queue_B[0].second.first; 
		soFar = queue_B[0].second.second; 

		// Pop the first element
		queue_B.erase(queue_B.begin()); 

		// Get where the edges from the Node in the given direction
		edges = nodes[v]->bwd_labelled_edges[0]; 

		// Insert into the Nodes visited so far the current Node
		soFar.insert(v); 

		// change label according to query type
		if (query_type == 2) {
			for (int i = 0; i < query_labels.size(); i ++ ) {
				labels_no.push_back(i); 
			}
		}

		else if (query_type == 3) {
			l--; 
			if (l == -1)
				l = query_labels.size() - 1; 
			labels_no.push_back(l); 
		}
		
		else if (query_type == 4) {
			if (l == query_labels.size()) {
				labels_no.push_back(l-1); 
			}
			else {
				labels_no.push_back(l); 
				if (l != 0) {
					labels_no.push_back(l - 1); 
				}
			}
		}

		// look at all the possible edges out of the node
		for (int i = 0; i < edges.size(); i ++ ) {
			int e = edges[i]; 

			// Condition path has a cycle
			if (soFar.find(e) != soFar.end()) {
				continue; 
			}

			// Search for query label in the labels of node under question
			vector<int> labels = nodes[e]->labels; 
			for (int i = 0; i < labels_no.size(); i ++ ) {
				l = labels_no[i]; 
				if (binary_search(labels.begin(), labels.end(), query_labels[l])) {

					// No meeting with forward BFS
					if (fwd.find(e) == fwd.end()) {
						bwd[e][l].push_back(soFar); 
						queue_B.push_back(make_pair(e, make_pair(l, soFar))); 
					}

					// Meeting with forward BFS
					else {
						// If same label then do a part by part matching to check whether 
						if (fwd[e].find(l) != fwd[e].end()) {

							// Set to check in
							vector < unordered_set < int > > matches = fwd[e][l]; 
							for (int i = 0; i < matches.size(); i ++ ) {

								//Iterate over the current set and find in the set to check in
								unordered_set < int > ::iterator it; 
								bool flag = true; 
								for (it = soFar.begin(); it != soFar.end(); it ++ ) {
									if (matches[i].find(*it) != matches[i].end()) {
										flag = false; 
										break; 
									}
								}
								if (flag) {
									return 1; 
								}
							}
						}
						
						bwd[e][l].push_back(soFar); 
						queue_B.push_back(make_pair(e, make_pair(l, soFar))); 
					}
				}
			}
		}
	}	
	return 0; 
}

int Graph::rr(int src, int dst, int query_type, vector<int> query_labels, int c_numWalk, int c_walkLength) {
	
	// walk number -> nodes visited in order
	vector<vector<int> > fwd_walk, bwd_walk; 

	 // current walk content
	unordered_set<int> fwd_set, bwd_set;

	// node -> label -> penalty
	unordered_map<int, unordered_map<int, int> > fwdCntr, bwdCntr; 
	
	// ordered vector of current walk
	vector<int> fwd, bwd;

	// node ->label-> walk numbers by which it has been reached
	unordered_map<int, unordered_map<int, vector<int> > > walkNumber_F, walkNumber_B; 

	 // walk counter
	int i_F = 0, i_B = 0;

	// walklength counter
	int j_F = 0, j_B = 0; 

	// Parameter initialization
	int e_F = src; 
	int label_F = -1; 

	int e_B = dst; 
	int label_B = query_labels.size(); 
	 

	Node* prevNode; 
	int prev, prevLabel;
	int ind; 
	int numChild;

	while((i_F < (c_numWalk*numWalks)/10) && (i_B < (c_numWalk*numWalks)/10)) { 
		prevNode = nodes[e_F]; 
		prev = e_F; 
		prevLabel = label_F; 
		fwd.push_back(e_F); 
		fwd_set.insert(e_F); 
		if (query_type == 3) {
			label_F ++ ; 
			if (label_F == query_labels.size()) {
				label_F = 0; 
			}
		}
		// Number of edges out of the previous node
		numChild = (prevNode->fwd_labelled_edges[0]).size(); 
		if (numChild) {
			ind = (rand_vector[rand_index]) % numChild; 
			rand_index = (rand_index + 1) % numNodes; 
			bool foundFlag = false; 
			int label_id = -1;
			if (label_F != -1) 
			 	label_id = query_labels[label_F]; 
			int k; 
			
			 // Finding next in line
			for (k = 0; k < numChild; k ++ ) {
				e_F = prevNode->fwd_labelled_edges[0][ind]; 
				if (query_type == 2) {
					label_F = 0; 
					Node* currNode = nodes[e_F]; 

					 //loop in the forward walk itself, this has been removed
					if (fwd_set.find(e_F) != fwd_set.end()) {
						ind = (ind + 1)%numChild; 
						continue; 
					}
					
					for (int i = 0; i < query_labels.size(); i ++ ) {
						if (binary_search(currNode->labels.begin(), currNode->labels.end(), query_labels[i])) {
							foundFlag = true; 
							if (walkNumber_B.find(e_F)!= walkNumber_B.end()) {
								if (walkNumber_B[e_F].find(label_F) != walkNumber_B[e_F].end()) {
									vector<int> bwd_matches = walkNumber_B[e_F][label_F]; 
									for (int m = 0; m < bwd_matches.size(); m ++ ) {
										bool flag = false; 
										vector<int> inQuestion; 
										vector<int>::iterator bwd_it; 
										if (bwd_matches[m] == i_B) {
											inQuestion = bwd; 
										}
										else
											inQuestion = bwd_walk[bwd_matches[m]]; 
										for (bwd_it = inQuestion.begin(); bwd_it != inQuestion.end(); ++ bwd_it) {											
											if (*bwd_it == e_F) {
												flag = true; 
												break; 
											}
											if (fwd_set.find(*bwd_it) != fwd_set.end()) {
												break; 
											}
										}
										if (flag) {
											return 1; 
										}
									}
								}
							}
							//Fill forward cycle
							walkNumber_F[e_F][label_F].push_back(i_F); 
							break; 
						}
					}
					if (foundFlag) break; 
				}

				else if (query_type == 3) {
					Node* currNode = nodes[e_F]; 
					if (binary_search(currNode->labels.begin(), currNode->labels.end(), label_id)) {
						if (fwdCntr[e_F][label_F] >= (c_numWalk*numWalks)/10) {
							ind = (ind + 1)%numChild; 
							continue; 
						}
						//loop in the forward walk itself, this has been removed
						else if (fwd_set.find(e_F) != fwd_set.end()) { 

							ind = (ind + 1)%numChild; 
							continue; 
						}

						else {
							foundFlag = true; 
							if (walkNumber_B.find(e_F)!= walkNumber_B.end()) {
								if (walkNumber_B[e_F].find(label_F) != walkNumber_B[e_F].end()) {
									vector<int> bwd_matches = walkNumber_B[e_F][label_F]; 
									
									for (int m = 0; m < bwd_matches.size(); m ++ ) {
										bool flag = false; 
										vector<int> inQuestion; 
										vector<int>::iterator bwd_it; 
										if (bwd_matches[m] == i_B) {
											inQuestion = bwd; 
										}
										else
											inQuestion = bwd_walk[bwd_matches[m]]; 
										for (bwd_it = inQuestion.begin(); bwd_it != inQuestion.end(); ++ bwd_it) {											
											if (*bwd_it == e_F) {
												flag = true; 
												break; 
											}
											if (fwd_set.find(*bwd_it) != fwd_set.end()) {
												break; 
											}
										}
										if (flag) {
											return 1; 
										}
									}
								}
							}
							//Fill forward cycle
							walkNumber_F[e_F][label_F].push_back(i_F); 
							break; 
						}
					}
				}

				else if (query_type == 4) {
					int backup_label = label_F;

					int label_id_next = -1;
					if (label_F+1 != query_labels.size()) label_id_next = query_labels[label_F + 1]; 
					Node* currNode = nodes[e_F]; 
					bool prevL = binary_search(currNode->labels.begin(), currNode->labels.end(), label_id); 
					bool nextL = binary_search(currNode->labels.begin(), currNode->labels.end(), label_id_next); 
					if (label_F == -1) prevL = false; 
					if (prevL && nextL) {
						if (rand() % 2 == 0) {
							label_F ++ ; 
							label_id = label_id_next; 
						}
					}

					else if (nextL) {
						label_F ++ ; 
						label_id = label_id_next; 
					}

					else if (!prevL && !nextL) {
						label_F = backup_label; 
						ind = (ind + 1)%numChild; 
						continue; 
					}

					if (fwdCntr[e_F][label_F] >= (c_numWalk*numWalks)/10) {
						label_F = backup_label; 
						ind = (ind + 1)%numChild; 
						continue; 
					}
					//loop in the forward walk itself, this has been removed
					else if (fwd_set.find(e_F) != fwd_set.end()) { 
						ind = (ind + 1)%numChild; 
						continue; 
					}

					else {
						foundFlag = true; 
						if (walkNumber_B.find(e_F)!= walkNumber_B.end()) {
							if (walkNumber_B[e_F].find(label_F) != walkNumber_B[e_F].end()) {
								vector<int> bwd_matches = walkNumber_B[e_F][label_F]; 
								
								for (int m = 0; m < bwd_matches.size(); m ++ ) {
									bool flag = false; 
									vector<int> inQuestion; 
									vector<int>::iterator bwd_it; 
									if (bwd_matches[m] == i_B) {
										inQuestion = bwd; 
									}
									else
										inQuestion = bwd_walk[bwd_matches[m]]; 
									for (bwd_it = inQuestion.begin(); bwd_it != inQuestion.end(); ++ bwd_it) {											
										if (*bwd_it == e_F) {
											flag = true; 
											break; 
										}
										if (fwd_set.find(*bwd_it) != fwd_set.end()) {
											break; 
										}
									}
									if (flag) {
										return 1; 
									}
								}
							}
						}
						//Fill forward cycle
						walkNumber_F[e_F][label_F].push_back(i_F); 
						break; 
					}
				}
				ind = (ind + 1) % numChild; 
			}
			if (!foundFlag) {
				if (k == numChild) {
					fwdCntr[prev][prevLabel] += (c_numWalk*numWalks)/20; 
				}
				e_F = -1; 
			}
		}
		else if (numChild == 0) {
			e_F = -1; 
		}
		if (e_F == -1) {
			e_F = src; 
			j_F = 0; 
			label_F = -1; 
			i_F ++ ; 
			if (fwdCntr[e_F][label_F] >= (c_numWalk*numWalks)/10) {
				return 0; 
			}
			fwd_walk.push_back(fwd); 
			fwd.clear(); 
			fwd_set.clear(); 
		}
		else {
			j_F ++ ; 
			if (j_F == (c_walkLength*walkLength)/10) {	
				e_F = src; 
				j_F = 0; 
				label_F = -1; 
				i_F ++ ; 
				fwd_walk.push_back(fwd); 
				fwd.clear(); 
				fwd_set.clear(); 
			}
		}

		prevNode = nodes[e_B]; 
		prev = e_B; 	
		prevLabel = label_B; 
		bwd.push_back(e_B); 
		bwd_set.insert(e_B); 
		numChild = (prevNode->bwd_labelled_edges[0]).size(); 
		if (numChild) {
			ind = (rand_vector[rand_index])%numChild; 
			rand_index = (rand_index + 1)%numNodes; 
			if (query_type == 3) {
				if (label_B == 0) {
					label_B = query_labels.size(); 
				}
				label_B--; 
			}
			bool foundFlag = false; 
			int label_id = -1;
			if (label_B != query_labels.size() && label_B != -1)
			 	label_id = query_labels[label_B]; 
			int k; 

			// Finding next in line
			for (k = 0; k < numChild; k ++ ) { 
				e_B = prevNode->bwd_labelled_edges[0][ind]; 
				if (query_type == 2) {
					label_B = 0; 
					Node* currNode = nodes[e_B]; 
					bool flag2 = false; 
					int l1 = 0, l2 = 0; 

					//loop in the backward walk itself, this has been removed
					if (bwd_set.find(e_B) != bwd_set.end()) { 
						ind = (ind + 1)%numChild; 
						continue; 
					}

					
					for (int i = 0; i < query_labels.size(); i ++ ) {
						if (binary_search(currNode->labels.begin(), currNode->labels.end(), query_labels[i])) {
							foundFlag = true; 

							if (walkNumber_F.find(e_B)!= walkNumber_F.end()) {
								if (walkNumber_F[e_B].find(label_B) != walkNumber_F[e_B].end()) {
									vector<int> fwd_matches = walkNumber_F[e_B][label_B]; 
									for (int m = 0; m < fwd_matches.size(); m ++ ) {
										bool flag = false; 
										vector<int> inQuestion; 
										vector<int>::iterator fwd_it; 
										if (fwd_matches[m] == i_F)
											inQuestion = fwd; 
										else
											inQuestion = fwd_walk[fwd_matches[m]]; 
										for (fwd_it = inQuestion.begin(); fwd_it != inQuestion.end(); fwd_it ++ ) {
											if (*fwd_it == e_B) {
												flag = true; 
												break; 
											}
											if (bwd_set.find(*fwd_it) != bwd_set.end()) {
												break; 
											}
										}
										if (flag) {
											return 1; 
										}
									}
								}
							}
							walkNumber_B[e_B][label_B].push_back(i_B); 
							//Fill backward cycle
							break; 
						}
					}
					if (foundFlag) break; 
				}
				else if (query_type == 3) {

					Node* currNode = nodes[e_B]; 
					if (binary_search(currNode->labels.begin(), currNode->labels.end(), label_id)) {
						if (bwdCntr[e_B][label_B] >= (c_numWalk*numWalks)/10) {
							ind = (ind + 1)%numChild; 
							continue; 
						}
						//loop in the backward walk itself, this has been removed
						else if (bwd_set.find(e_B) != bwd_set.end()) { 
							ind = (ind + 1)%numChild; 
							continue; 
						}
						else{
							foundFlag = true; 
							if (walkNumber_F.find(e_B)!= walkNumber_F.end()) {
								if (walkNumber_F[e_B].find(label_B) != walkNumber_F[e_B].end()) {
									vector<int> fwd_matches = walkNumber_F[e_B][label_B]; 
									for (int m = 0; m < fwd_matches.size(); m ++ ) {
										bool flag = false; 
										vector<int> inQuestion; 
										vector<int>::iterator fwd_it; 
										if (fwd_matches[m] == i_F)
											inQuestion = fwd; 
										else
											inQuestion = fwd_walk[fwd_matches[m]]; 
										for (fwd_it = inQuestion.begin(); fwd_it != inQuestion.end(); fwd_it ++ ) {
											if (*fwd_it == e_B) {
												flag = true; 
												break; 
											}
											if (bwd_set.find(*fwd_it) != bwd_set.end()) {
												break; 
											}
										}
										if (flag) {
											return 1; 
										}
									}
								}
							}
							walkNumber_B[e_B][label_B].push_back(i_B); 
							//Fill backward cycle
							break; 
						}
					}
				}
				else if (query_type == 4) {
					int restore_label = label_B;
					int label_id_next = 0;
					if (label_B != 0)
						label_id_next = query_labels[label_B - 1]; 
					Node* currNode = nodes[e_B]; 
					bool prevL = binary_search(currNode->labels.begin(), currNode->labels.end(), label_id);
					bool nextL;
					if (label_B == -1) nextL = false; 
					else  nextL = binary_search(currNode->labels.begin(), currNode->labels.end(), label_id_next); 
					if (prevL && nextL) {
						if (rand() % 2 == 0) {
							label_B--; 
							label_id = label_id_next; 
						}
					}

					else if (nextL) {
						label_B--; 
						label_id = label_id_next; 
					}

					else if (!prevL && !nextL) {
						ind = (ind + 1)%numChild; 
						continue; 
					}

					if (bwdCntr[e_B][label_B] >= (c_numWalk*numWalks)/10) {
						label_B = restore_label; 
						ind = (ind + 1)%numChild; 
						continue; 
					}

					//loop in the backawrd walk itself, this has been removed
					else if (bwd_set.find(e_B) != bwd_set.end()) { 
						label_B = restore_label; 
						ind = (ind + 1)%numChild; 
						continue; 
					}

					else{
						foundFlag = true; 
						if (walkNumber_F.find(e_B)!= walkNumber_F.end()) {
							if (walkNumber_F[e_B].find(label_B) != walkNumber_F[e_B].end()) {
								vector<int> fwd_matches = walkNumber_F[e_B][label_B]; 
								for (int m = 0; m < fwd_matches.size(); m ++ ) {
									bool flag = false; 
									vector<int> inQuestion; 
									vector<int>::iterator fwd_it; 
									if (fwd_matches[m] == i_F)
										inQuestion = fwd; 
									else
										inQuestion = fwd_walk[fwd_matches[m]]; 
									for (fwd_it = inQuestion.begin(); fwd_it != inQuestion.end(); ++fwd_it) {
										if (*fwd_it == e_B) {
											flag = true; 
											break; 
										}
										if (bwd_set.find(*fwd_it) != bwd_set.end()) {
											break; 
										}
									}
									if (flag) { 
										return 1; 
									}
								}
							}
						}
						walkNumber_B[e_B][label_B].push_back(i_B); 
						//Fill backward cycle
						break; 
					}
				}
				ind = (ind + 1)%numChild; 
			}
			if (!foundFlag) {
				if (k == numChild) {
					bwdCntr[prev][prevLabel] += (c_numWalk*numWalks)/20; 
				}
				e_B = -1; 
			}
		}
		else if (numChild == 0) {
			e_B = -1; 
		}

		if (e_B == -1) {
			e_B = dst; 
			j_B = 0; 
			label_B = query_labels.size(); 
			i_B ++ ; 	
			if (bwdCntr[e_B][label_B] >= (c_numWalk*numWalks)/10) {
				return 0; 
			}
			bwd_walk.push_back(bwd); 
			bwd.clear(); 
			bwd_set.clear(); 
		}
		else {
			j_B ++ ; 
			if (j_B == (c_walkLength*walkLength)/10) {	
				e_B = dst; 
				j_B = 0; 
				label_B = query_labels.size(); 
				i_B ++ ; 
				bwd_walk.push_back(bwd); 
				bwd.clear(); 
				bwd_set.clear(); 
			}
		}
	}
	return 0; 
}