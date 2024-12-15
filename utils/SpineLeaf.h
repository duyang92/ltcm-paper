#ifndef _SPINELEAF_H
#define _SPINELEAF_H

#include <algorithm>
#include <string>
#include <vector>
using namespace std;

class SpineLeaf {
   public:
    int n;
    int m;
    int iSpine;
    int iLeaf;
    int switches_num;
    vector<vector<int>> SpineSwitch;
    vector<vector<int>> LeafSwitch;
    vector<vector<int>> switches;
    vector<vector<int>> path;
    int i_1_hop;
    int i_3_hop;
    vector<int> betwenness_centrality;

    SpineLeaf(int n = 4, int m = 8) {
        this->n = n;
        this->m = m;
        this->iSpine = n;
        this->iLeaf = m;
        this->switches_num = iSpine + iLeaf + 1;
        this->SpineSwitch = vector<vector<int>>(iSpine, vector<int>());
        this->LeafSwitch = vector<vector<int>>(iLeaf, vector<int>());
        this->switches = vector<vector<int>>(iSpine + iLeaf + 1, vector<int>());
        this->betwenness_centrality = vector<int>(iSpine + iLeaf + 1, 0);
        this->i_1_hop = m;
        this->i_3_hop = n * m * (m - 1) / 2;
        build_paths();
        get_betweenness_centrality();
    }

    void build_paths() {
        for (int i = 1; i <= iSpine; ++i) {
            for (int j = iSpine + 1; j <= iSpine + iLeaf; ++j) {
                for (int k = j + 1; k <= iSpine + iLeaf; ++k) {
                    path.push_back({j, i, k});
                }
            }
        }
        for (int i = iSpine + 1; i <= iSpine + iLeaf; ++i) {
            path.push_back({i});
        }

        vector<vector<int>> temp_p;
        for (auto& i : path) {
            temp_p.push_back(i);
            if (i.size() > 1) {
                vector<int> reversed_i = i;
                reverse(reversed_i.begin(), reversed_i.end());
                temp_p.push_back(reversed_i);
            }
        }
        sort(temp_p.begin(), temp_p.end(), [](const vector<int>& a, const vector<int>& b) {
            return a.size() < b.size();
        });
        path = temp_p;
    }

    void get_betweenness_centrality() {
        for (auto& i : path) {
            for (auto& j : i) {
                betwenness_centrality[j]++;
            }
        }
        for (int i = 1; i <= iSpine; ++i) {
            betwenness_centrality[i] /= 2;
        }
        for (int i = iSpine + 1; i <= iSpine + iLeaf; ++i) {
            betwenness_centrality[i] = (betwenness_centrality[i] + 1) / 2;
        }
    }
};

#endif  // _SPINELEAF_H