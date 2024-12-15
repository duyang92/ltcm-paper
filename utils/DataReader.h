#ifndef _DATAREADER_H_
#define _DATAREADER_H_

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

using namespace std;

class DataReader {
    typedef pair<uint32_t, uint32_t> pii;
    typedef pair<uint64_t, uint64_t> pll;

   private:
    string caida_dir_path;
    vector<pii> data_set, cur_data;
    unordered_map<uint32_t, int> actual_flow_size, cur_size;
    unordered_map<uint32_t, int> actual_flow_spread, cur_spread;
    unordered_map<uint32_t, unordered_set<uint32_t>> actual_flow_set;
    vector<uint32_t> keys;
    int is_spread;

    pii* packets;
    uint32_t* flows;

   public:
    DataReader() {
    }
    DataReader(string caida_dir_path, int test, int test_num, int is_spread = 0) {
        this->caida_dir_path = caida_dir_path;
        this->is_spread = is_spread;
        data_set.clear();
        actual_flow_size.clear();
        actual_flow_spread.clear();
        read_file(test, test_num);
        read_size_and_spread();
    }
    void read_file(int test, int test_num) {
        if (filesystem::exists(caida_dir_path) == false) {
            cout << "Directory not found!" << endl;
            exit(0);
        }
        ifstream inf;
        int count = 0;
        pii temp;
        for (const auto& entry : filesystem::directory_iterator(caida_dir_path)) {
            inf.open(entry.path(), ios::in);
            cout << "Reading " << entry.path() << endl;
            while (inf >> temp.first && inf >> temp.second) {
                count++;
                if (is_spread == 0) temp.second = count;
                else {
                    uint32_t tptp = temp.second;
                    temp.second = temp.first;
                    temp.first = tptp;
                }
                data_set.push_back(temp);
                if (test && count == test_num) break;
            }
            inf.close();
            // break;
        }
        cur_data = data_set;
    }
    void read_size_and_spread() {
        for (int i = 0; i < data_set.size(); i++) {
            if (actual_flow_size.find(data_set[i].first) == actual_flow_size.end()) {
                if (is_spread == 1) {
                    actual_flow_set[data_set[i].first].insert(data_set[i].second);
                }
                actual_flow_size[data_set[i].first] = 1;
                keys.push_back(data_set[i].first);
            } else {
                if (is_spread == 1) {
                    actual_flow_set[data_set[i].first].insert(data_set[i].second);
                }
                actual_flow_size[data_set[i].first]++;
            }
        }
        int max_spread = 0, max_size = 0;
        for (const auto& flow : actual_flow_set) {
            actual_flow_spread[flow.first] = actual_flow_set[flow.first].size();
            max_spread = max(max_spread, actual_flow_spread[flow.first]);
            max_size = max(max_size, actual_flow_size[flow.first]);
        }
        cur_size = actual_flow_size;
        cur_spread = actual_flow_spread;
    }
    unordered_map<uint32_t,int> add_abnormal_flow(int num, int size) {
        unordered_map<uint32_t, int> res;
        while (num) {
            uint32_t flow_id_left = rand() % (1 << 16);
            uint32_t flow_id_right = rand() % (1 << 16);
            uint32_t flow_id = (flow_id_left << 16) | flow_id_right;
            if (actual_flow_size.find(flow_id) == actual_flow_size.end()) {
                actual_flow_size[flow_id] = size;
                res[flow_id] = size;
                num--;
            }
        }
        return res;
    }
    vector<pii> get_data_set(int flows = -1) {
        if (flows == -1 or is_spread == 1) {
            return cur_data;
        }
        srand(time(0));
        random_shuffle(keys.begin(), keys.end());
        cur_data.clear();
        cur_size.clear();
        int count = 0;
        for (int i = 0; i < flows; i++) {
            uint32_t f = keys[i];
            cur_size[f] = actual_flow_size[f];
            for (int j = 0; j < actual_flow_size[f]; j++) {
                cur_data.push_back({f, ++count});
            }
        }
        random_shuffle(cur_data.begin(), cur_data.end());
        return cur_data;
    }
    unordered_map<uint32_t, int> get_actual_flow_size() {
        return cur_size;
    }
    unordered_map<uint32_t, int> get_actual_flow_spread() {
        return cur_spread;
    }
};

#endif