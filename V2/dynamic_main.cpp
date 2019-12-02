#include "GraphDyn.h"

int main(int argc, char *argv[]){
	char *graphFile = argv[1];
	char *queryFile = argv[2];
	int c = 1, c_walkLength = 1, c_numWalks = 1, c_budget = 1;
	bool fractional = false;

	srand(time(NULL));
	Graph *tg = new Graph(graphFile, 1);
	tg->c_numWalks =  c_numWalks;
	tg->c_walkLength = c_walkLength;

	int numWalks_input =c_numWalks;
	int walkLength_input = c_walkLength;

	tg->dynamicQuery(queryFile);
}
