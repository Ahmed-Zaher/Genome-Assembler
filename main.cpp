#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

#include "Assembler.h"

int main() {
	ios::sync_with_stdio(false);
	Assembler assembler;
	assembler.read();
	assembler.assemble();
	vector<string> contigs = assembler.getContigs();
	for (int i = 0; i < (int) contigs.size(); ++i)
		cout << "Contig " << i + 1 << ": " << contigs[i] << '\n';
	return 0;
}

