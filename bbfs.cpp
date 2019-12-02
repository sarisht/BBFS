#include "GraphSub.h"

int main(int argc, char *argv[]){
	string graphname = argv[1];
	char *graphFile = argv[2];
	char *featFile = argv[3];
	int dir_control = atoi(argv[4]);
	int numQueries = 1;
	float c_walkLength = 2;
	float c_numWalk = 1;
	int c = 1, c_budget = 1;
	bool fractional = false;

	srand(time(NULL));
	Graph *tg = new Graph(graphFile, featFile, dir_control);
	outFile1.open(graphname+"BBFS.txt");
	char *queryFile = argv[5];
	tg->readAndRunBFS(queryFile);
}
