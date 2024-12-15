#ifndef _FATTREE_H
#define _FATTREE_H

#include <vector>
#include <string>
#include <algorithm>
using namespace std;

class FatTree {
   public:
    int k;
    int pod;
    int t;
    int iCore;
    int iAgg;
    int iEdge;
    int iHost;
    vector<vector<int>> CoreSwitch;
    vector<vector<int>> AggSwitch;
    vector<vector<int>> EdgeSwitch;
    vector<int> CoreId;
    vector<int> AggId;
    vector<int> EdgeId;
    int switches_num;
    vector<vector<int>> switches;
    vector<vector<int>> path;
    int i_1_hop;
    int i_3_hop;
    int i_5_hop;
    vector<int> betweeness_centrality;

   public:
    FatTree(int input_k = 4) {
        k = input_k;
        pod = k;
        t = k / 2;
        iCore = (k / 2) * (k / 2);
        iAgg = (k / 2) * k;
        iEdge = (k / 2) * k;
        iHost = (k / 2) * (k / 2) * (k / 2);
        switches_num = iCore + iAgg + iEdge + 1;

        switches = vector<vector<int>>(iCore + iAgg + iEdge + 1);
        CoreSwitch = vector<vector<int>>(iCore);
        AggSwitch = vector<vector<int>>(iAgg);
        EdgeSwitch = vector<vector<int>>(iEdge);
        CoreId = vector<int>(iCore);
        for (int i = 0; i < iCore; ++i) {
            CoreId[i] = i + 1;
        }
        AggId = vector<int>(iAgg);
        for (int i = 0; i < iAgg; ++i) {
            AggId[i] = i + iCore + 1;
        }
        EdgeId = vector<int>(iEdge);
        for (int i = 0; i < iEdge; ++i) {
            EdgeId[i] = i + iCore + iAgg + 1;
        }

        i_1_hop = k * k / 2;
        i_3_hop = (k * k * k * k - 2 * k * k * k) / 16;
        i_5_hop = (k * k * k * k * k * k - k * k * k * k * k) / 32;

        build_edges();
        build_paths();

        betweeness_centrality = vector<int>(iCore + iAgg + iEdge + 1, 0);
        get_betweenness_centrality();
    }

    ~FatTree() {
    }

    void build_edges() {
        for (int i = 0; i < iEdge; ++i) {
            int pod = i / t;
            for (int j = 0; j < t; ++j) {
                EdgeSwitch[i].push_back(AggId[pod * t + j]);
                AggSwitch[pod * t + j].push_back(EdgeId[i]);
            }
        }
        for (int p = 0; p < pod; ++p) {
            int core_idx = 0;
            for (int i = 0; i < t; ++i) {
                for (int j = 0; j < t; ++j) {
                    AggSwitch[p * t + i].push_back(CoreId[core_idx]);
                    CoreSwitch[core_idx].push_back(AggId[p * t + i]);
                    ++core_idx;
                }
            }
        }
        for (int i = 0; i < iCore; ++i) {
            switches[CoreId[i]] = CoreSwitch[i];
        }
        for (int i = 0; i < iAgg; ++i) {
            switches[AggId[i]] = AggSwitch[i];
        }
        for (int i = 0; i < iEdge; ++i) {
            switches[EdgeId[i]] = EdgeSwitch[i];
        }
    }

    void dfs(int u, int fa, vector<int> current_path, int f) {
        current_path.push_back(u);
        if (u >= iCore + iAgg + 1) {
            path.push_back(current_path);
            return;
        }
        for (int v : switches[u]) {
            if (v == fa) {
                continue;
            }
            if (u <= iCore) {
                dfs(v, u, current_path, 1);
            } else {
                if (f == 0) {
                    if (v > u) {
                        dfs(v, u, current_path, 1);
                    } else {
                        dfs(v, u, current_path, 0);
                    }
                } else if (f == 1 && v > u) {
                    dfs(v, u, current_path, 0);
                }
            }
        }
    }

    void build_paths() {
        int maxSwitch = iEdge;
        for (int i = 0; i < maxSwitch; ++i) {
            path.push_back({EdgeId[i]});
            for (int j : EdgeSwitch[i]) {
                vector<int> current_path;
                current_path.push_back(EdgeId[i]);
                dfs(j, EdgeId[i], current_path, 0);
            }
        }
        vector<vector<int>> tp = path;
        path.clear();
        for (vector<int> i : tp) {
            vector<int> cur_no = i;
            sort(i.begin(), i.end());
            string cur = "";
            for (int j : i) {
                cur += to_string(j) + " ";
            }
            int f = 1;
            for (vector<int> j : path) {
                sort(j.begin(), j.end());
                string str_j = "";
                for (int k : j) {
                    str_j += to_string(k) + " ";
                }
                if (cur == str_j) {
                    f = 0;
                    break;
                }
            }
            if (f == 1) {
                path.push_back(cur_no);
            }
        }
        tp = path;
        for (vector<int> i : path) {
            if (i.size() == 1) continue;
            vector<int> cur_no = i;
            reverse(cur_no.begin(), cur_no.end());
            tp.push_back(cur_no);
        }
        path = tp;
        sort(path.begin(), path.end(), [](const vector<int>& x, const vector<int>& y) {
            return x.size() < y.size();
        });
        // for (vector<int> i : path) {
        //     for (int j : i) {
        //         cout << j << " ";
        //     }
        //     cout << endl;
        // }
    }

    void get_betweenness_centrality() {
        for (vector<int> p : path) {
            for (int node : p) {
                betweeness_centrality[node]++;
            }
        }
        for (int i = 1;i <= iCore+iAgg;i++){
            betweeness_centrality[i] = betweeness_centrality[i] / 2;
        }
        for (int i = iCore+iAgg+1;i <= iCore+iAgg+iEdge;i++){
            betweeness_centrality[i] = (betweeness_centrality[i]+1)/ 2;
        }
    }

};

#endif  // _FATTREE_H