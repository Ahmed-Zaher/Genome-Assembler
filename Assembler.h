class Assembler {

private:

	map<string, int> kMerToInt;
	map<int, string> intToKMer;
	set<pair<int, int>> edges, visEdges;
	map<pair<int, int>, int> coverage;
	vector<string> reads;
	vector<vector<vector<int>>> adj;
	vector<int> inDeg, outDeg, vis;
	vector<vector<int>> paths;
	vector<set<pair<int, int>>> inOutDegSet;
	int k, t, n, visCnt;

	void getKMers() {
		k = 20;
		n = 0;
		for (auto& read : reads) {
			set<pair<int, int>> kMers;
			for (int i = 0; i <= (int) read.size() - k; ++i) {
				string kMer = read.substr(i, k);
				string kMerPref = kMer.substr(0, k - 1);
				string kMerSuff = kMer.substr(1);
				if (kMerToInt.find(kMerPref) == kMerToInt.end())
					kMerToInt[kMerPref] = n, intToKMer[n++] = kMerPref;
				if (kMerToInt.find(kMerSuff) == kMerToInt.end())
					kMerToInt[kMerSuff] = n, intToKMer[n++] = kMerSuff;
				int u = kMerToInt[kMerPref];
				int v = kMerToInt[kMerSuff];
				kMers.insert(make_pair(u, v));
				if (gotEdge(u, v))
					continue;
				edges.insert(make_pair(u, v));
			}
			for (auto& kMer : kMers)
				++coverage[kMer];
		}
		adj.assign(2, vector<vector<int>>(n, vector<int>()));
		inDeg.assign(n, 0);
		outDeg.assign(n, 0);
		for (auto& edge : edges) {
			int u = edge.first, v = edge.second;
			adj[0][u].push_back(v);
			adj[1][v].push_back(u);
			++outDeg[u];
			++inDeg[v];
		}
		inOutDegSet.assign(2, set<pair<int, int>>());
		for (int i = 0; i < n; ++i) {
			inOutDegSet[0].insert(make_pair(inDeg[i], i));
			inOutDegSet[1].insert(make_pair(outDeg[i], i));
		}
	}

	void removeTips() {
		vis.assign(n, 0);
		while ((!inOutDegSet[0].empty() && !inOutDegSet[0].begin()->first)
			|| (!inOutDegSet[1].empty() && !inOutDegSet[1].begin()->first))
			for (int i = 0; i < 2; ++i)
				while(!inOutDegSet[i].empty()) {
					int currDeg = inOutDegSet[i].begin()->first;
					int u = inOutDegSet[i].begin()->second;
					if (vis[u]) {
						inOutDegSet[i].erase(inOutDegSet[i].begin());
						continue;
					}
					if (currDeg)
						break;
					vis[u] = true;
					for (auto& v : adj[i][u])
						if (!vis[v]) {
							if (!i) {
								removeEdge(u, v);
								inOutDegSet[i].insert(make_pair(inDeg[v], v));
							} else {
								removeEdge(v, u);
								inOutDegSet[i].insert(make_pair(outDeg[v], v));
							}
						}
				}
	}

	int DFS(int node, bool type, int depth, int u, int v) {
		if (depth > t)
			return 0;
		if (type && node == v) {
			double avg[2] = {0, 0};
			for (int i = 0; i < 2; ++i) {
				for (int j = 0; j < (int) paths[i].size() - 1; ++j)
					avg[i] += coverage[make_pair(paths[i][j], paths[i][j + 1])];
				avg[i] /= (int) paths[i].size() - 1;
			}
			int del = (avg[0] < avg[1] ? 0 : 1);
			for (int j = 0; j < (int) paths[del].size() - 1; ++j)
				removeEdge(paths[del][j], paths[del][j + 1]);
			return 1;
		}
		vis[node] = visCnt;
		if (!type && node != u) {
			if (DFS(u, 1, 0, u, node))
				return 1;
		}
		for (auto& child : adj[0][node])
			if (gotEdge(node, child) && (vis[child] != visCnt || (depth && child == v))) {
				paths[type].push_back(child);
				if (DFS(child, type, depth + 1, u, v))
					return 1;
				paths[type].pop_back();
			}
		vis[node] = 0;
		return 0;
	}

	void removeBubbles() {
		vis.assign(n, 0);
		paths.assign(2, vector<int>());
		t = k;
		for (int i = 0; i < n; ++i) {
			do {
				paths[0] = paths[1] = {i};
				++visCnt;
			} while(DFS(i, 0, 0, i, -1));
		}
	}

	void removeEdge(int u, int v) {
		--outDeg[u];
		--inDeg[v];
		edges.erase(make_pair(u, v));
	}

	bool gotEdge(int u, int v) {
		return edges.find(make_pair(u, v)) != edges.end();
	}

	bool visitedEdge(int u, int v) {
		return visEdges.find(make_pair(u, v)) != visEdges.end();
	}

	int getNextNode(int u) {
		for (auto& v : adj[0][u])
			if (gotEdge(u, v))
				return v;
	}


	vector<vector<int>> getContigsPaths() {
		vector<vector<int>> contigsPaths;
		// non-braching paths
		for (int u = 0; u < n; ++u) {
			if (inDeg[u] != 1 || outDeg[u] != 1)
				for (auto& v : adj[0][u])
					if (gotEdge(u, v)) {
						vector<int> contigPath = {u, v};
						visEdges.insert(make_pair(u, v));
						while (inDeg[v] == 1 && outDeg[v] == 1) {
							int w = getNextNode(v);
							visEdges.insert(make_pair(v, w));
							v = w;
							contigPath.push_back(v);
						}
						contigsPaths.push_back(contigPath);
					}
		}
		// isolated cycles
		for (auto& edge : edges) {
			int u = edge.first, v = edge.second;
			if (!visitedEdge(u, v)) {
				vector<int> contigPath = {u, v};
				visEdges.insert(make_pair(u, v));
				int nxtNode = getNextNode(v);
				while (nxtNode != u) {
					visEdges.insert(make_pair(v, nxtNode));
					v = nxtNode;
					contigPath.push_back(v);
					nxtNode = getNextNode(v);
				}
				visEdges.insert(make_pair(v, nxtNode));
				contigPath.push_back(u);
				contigsPaths.push_back(contigPath);
			}
		}
		return contigsPaths;
	}

public:

	void read() {
		string read;
		while (cin >> read)
			reads.push_back(read);
	}

	void assemble() {
		getKMers();
		removeTips();
		removeBubbles();
	}

	vector<string> getContigs() {
		vector<string> contigs;
		vector<vector<int>> contigsPaths = getContigsPaths();
		for (auto& contigPath : contigsPaths) {
			string contig;
			if (contigPath.front() == contigPath.back()) {
				for (int i = 0; i < (int) contigPath.size() - 1; ++i)
					contig.push_back(intToKMer[contigPath[i]].back());
			} else {
				contig = intToKMer[contigPath[0]];
				for (int i = 1; i < (int) contigPath.size(); ++i)
					contig.push_back(intToKMer[contigPath[i]].back());
			}
			contigs.push_back(contig);
		}
		return contigs;
	}
};
