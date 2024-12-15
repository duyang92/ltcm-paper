#ifndef _CHAIN_OUR_
#define _CHAIN_OUR_

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <bitset>
#include "../../utils/FatTree.h"
#include "../../utils/SpineLeaf.h"
#include "../../utils/MurmurHash3.h"

class CHAIN_OUR {
#define MAXI 4294967295

#define MAX_FILTER_COUNTER 100000
#define FILTER_W 1024 * 2
#define FILTER_D 1
#define FILTER_COUNTER_SIZE 4
#define FILTER_FINGERPRINT_SIZE 16
#define MEM_QRY_BIT_LEN 65536

    // #define RE_DIS_THREOLD 10000

#define HASH_SEED 419
#define KEY_SIZE 4

    typedef pair<uint32_t, uint32_t> pii;

   public:
    uint32_t H(uint32_t t, uint32_t s) {
        uint32_t hash_val = 0;
        char hash_input_str[5] = {0};
        memcpy(hash_input_str, &t, sizeof(uint32_t));
        MurmurHash3_x86_32(hash_input_str, 4, s, &hash_val);
        return hash_val;
    }

    uint16_t FP(uint32_t hash) {
        hash ^= hash >> 16;
        hash *= 0x85ebca6b;
        hash ^= hash >> 13;
        hash *= 0xc2b2ae35;
        hash ^= hash >> 16;
        return hash & 65535;
    }

    struct Bucket {
        uint32_t key = 0;
        uint32_t C = 0;
    };

    struct Cell {
        uint16_t flow_id;
        uint8_t size;  // actually use 4 bits
        Cell() {
            flow_id = 0;
            size = 0;
        }
        ~Cell() {}
    };

    class hg_node {
       public:
        Bucket** bucketArray;
        uint32_t ROW_NUM;
        uint32_t COL_NUM;
        uint32_t MAX_CHAIN_OUR_LENGTH = 3;
        mt19937 rng = mt19937(time(0));
        hg_node() {}
        hg_node(uint32_t row, uint32_t col) {
            ROW_NUM = row;
            COL_NUM = col;
            bucketArray = new Bucket*[ROW_NUM];
            for (uint32_t i = 0; i < ROW_NUM; i++) {
                bucketArray[i] = new Bucket[COL_NUM];
                memset(bucketArray[i], 0, sizeof(Bucket) * COL_NUM);
            }
        }
        uint32_t H(uint32_t t, uint32_t s) {
            uint32_t hash_val = 0;
            char hash_input_str[5] = {0};
            memcpy(hash_input_str, &t, sizeof(uint32_t));
            MurmurHash3_x86_32(hash_input_str, 4, s, &hash_val);
            return hash_val;
        }
        void insert(uint32_t key) {
            uint32_t min_counter = UINT32_MAX;
            uint32_t min_rowIdx = -1;
            uint32_t min_colIdx = -1;

            for (uint32_t rowIndex = 0; rowIndex < ROW_NUM; rowIndex++) {
                // unsigned int hashValue = 0;
                // MurmurHash3_x86_32(key, KEY_SIZE, HASH_SEED + 6156137 * rowIndex, &hashValue);
                uint32_t hashValue = H(key, HASH_SEED + 6156137 * rowIndex);

                uint32_t colIndex = hashValue % COL_NUM;
                if (bucketArray[rowIndex][colIndex].C == 0) {
                    bucketArray[rowIndex][colIndex].key = key;
                    bucketArray[rowIndex][colIndex].C = 1;
                    return;
                } else if (bucketArray[rowIndex][colIndex].key == key) {
                    bucketArray[rowIndex][colIndex].C++;
                    return;
                }

                if (bucketArray[rowIndex][colIndex].C < min_counter) {
                    min_rowIdx = rowIndex;
                    min_colIdx = colIndex;
                    min_counter = bucketArray[rowIndex][colIndex].C;
                }
            }

            double randomVal = (double)rng() / 4294967296;
            if (randomVal <= 1.0 / (min_counter + 1)) {
                // the threshold is set to 511
                if (min_counter >= 512) {
                    // unsigned int hashValue = 0;
                    // MurmurHash3_x86_32(key, KEY_SIZE, HASH_SEED + 21151121, &hashValue);
                    uint32_t hashValue = H(key, HASH_SEED + 21151121);
                    uint32_t p2 = (min_colIdx + hashValue) % COL_NUM;

                    unordered_set<uint32_t> tempSet;
                    tempSet.insert(min_colIdx);

                    uint32_t l = 2;
                    uint32_t bb = bucketArray[min_rowIdx][p2].C;
                    while (l <= MAX_CHAIN_OUR_LENGTH && tempSet.find(p2) == tempSet.end() && min_counter > bucketArray[min_rowIdx][p2].C) {
                        tempSet.insert(p2);

                        double randomVal = (double)rng() / 4294967296;
                        if (randomVal <= pow(((double)min_counter / (min_counter + bb)), l)) {
                            bucketArray[min_rowIdx][p2].C = min_counter;
                            bucketArray[min_rowIdx][p2].key = bucketArray[min_rowIdx][min_colIdx].key;
                            break;
                        } else {
                            p2 = (p2 + hashValue) % COL_NUM;
                            bb = bucketArray[min_rowIdx][p2].C;
                            l++;
                        }
                    }
                }

                bucketArray[min_rowIdx][min_colIdx].C++;
                bucketArray[min_rowIdx][min_colIdx].key = key;
            }
        }

        void getEstimatedFlowSizes(unordered_map<uint32_t, uint32_t>& estimatedFlowSizes) {
            for (uint32_t rowIndex = 0; rowIndex < ROW_NUM; rowIndex++) {
                for (uint32_t colIndex = 0; colIndex < COL_NUM; colIndex++) {
                    if (estimatedFlowSizes.find(bucketArray[rowIndex][colIndex].key) == estimatedFlowSizes.end()) {
                        uint32_t key = bucketArray[rowIndex][colIndex].key;
                        estimatedFlowSizes[key] = bucketArray[rowIndex][colIndex].C;
                    } else {
                        estimatedFlowSizes[bucketArray[rowIndex][colIndex].key] += bucketArray[rowIndex][colIndex].C;
                    }
                }
            }
        }

        uint32_t query(uint32_t f) {
            uint32_t size = 1e9;
            for (uint32_t rowIndex = 0; rowIndex < ROW_NUM; rowIndex++) {
                uint32_t hashValue = H(f, HASH_SEED + 6156137 * rowIndex);
                uint32_t colIndex = hashValue % COL_NUM;
                size = min(size, bucketArray[rowIndex][colIndex].C);
            }
            return size;
        }
    };

    // memory config
    uint32_t totalMem;
    uint32_t col_num;
    uint32_t row_num;

    // topology
    int k = 0;
    FatTree ft;
    int n, m;
    SpineLeaf sl;
    int switch_num, path_num;
    vector<vector<int>> path;
    hg_node* swtich_hg;
    unordered_map<uint32_t, uint32_t> per_flow_size;

    // filter sketch
    uint32_t reset;
    uint32_t threshold;
    double filter_ratio;
    uint32_t filter_w, filter_d;
    double filter_memory;
    Cell** filter_cm;
    uint32_t* counter;

    // flow distribution
    pii** interval;
    uint32_t RE_DIS_THREOLD;

    // flow path
    unordered_map<uint32_t, int> flow_path_id;

    // seed
    uint32_t* seed;
    uint32_t index_seed;

    // load
    vector<int> packet_load;
    vector<unordered_set<uint32_t>> flow_load;

    CHAIN_OUR(int k, double total_memory, uint32_t* seeds, uint32_t re_dis_th = 10000) {
        srand((int)time(0));

        this->k = k;
        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();

        // filter sketch
        this->reset = MAX_FILTER_COUNTER;
        this->filter_ratio = 0.001;
        this->threshold = this->reset * filter_ratio;
        this->filter_w = FILTER_W;
        this->filter_d = FILTER_D;
        this->filter_memory = (FILTER_COUNTER_SIZE + FILTER_FINGERPRINT_SIZE) * filter_d * filter_w / 1024 / 8;
        this->filter_cm = new Cell*[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            filter_cm[i] = new Cell[filter_w];
            memset(filter_cm[i], 0, sizeof(Cell) * filter_w);
        }
        this->counter = new uint32_t[switch_num];
        memset(counter, 0, sizeof(uint32_t) * switch_num);

        this->totalMem = total_memory * 8 * 1024;  // the total number of bits
        this->row_num = 4;
        this->col_num = 1.0 * totalMem / row_num / (32 + 24);

        this->swtich_hg = new hg_node[switch_num];
        for (uint32_t i = 0; i < switch_num; ++i) {
            swtich_hg[i] = hg_node(row_num, col_num);
        }

        // flow distribution
        this->interval = new pii*[path_num];
        for (int i = 0; i < path_num; ++i) {
            interval[i] = new pii[path[i].size()];
            memset(interval[i], 0, sizeof(pii) * path[i].size());
        }

        this->RE_DIS_THREOLD = re_dis_th;

        this->seed = seeds;
        this->index_seed = 131;

        // load
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config_FatTree();
    }

    CHAIN_OUR(int n, int m, double total_memory, uint32_t* seeds, uint32_t re_dis_th = 10000) {
        srand((int)time(0));

        this->n = n;
        this->m = m;
        this->sl = SpineLeaf(n, m);
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();

        // filter sketch
        this->reset = MAX_FILTER_COUNTER;
        this->filter_ratio = 0.001;
        this->threshold = this->reset * filter_ratio;
        this->filter_w = FILTER_W;
        this->filter_d = FILTER_D;
        this->filter_memory = (FILTER_COUNTER_SIZE + FILTER_FINGERPRINT_SIZE) * filter_d * filter_w / 1024 / 8;
        this->filter_cm = new Cell*[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            filter_cm[i] = new Cell[filter_w];
            memset(filter_cm[i], 0, sizeof(Cell) * filter_w);
        }
        this->counter = new uint32_t[switch_num];
        memset(counter, 0, sizeof(uint32_t) * switch_num);

        this->totalMem = total_memory * 8 * 1024;  // the total number of bits
        this->row_num = 4;
        this->col_num = 1.0 * totalMem / row_num / (32 + 24);

        this->swtich_hg = new hg_node[switch_num];
        for (uint32_t i = 0; i < switch_num; ++i) {
            swtich_hg[i] = hg_node(row_num, col_num);
        }

        // flow distribution
        this->interval = new pii*[path_num];
        for (int i = 0; i < path_num; ++i) {
            interval[i] = new pii[path[i].size()];
            memset(interval[i], 0, sizeof(pii) * path[i].size());
        }

        this->RE_DIS_THREOLD = re_dis_th;

        this->seed = seeds;
        this->index_seed = 131;

        // load
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config_SpineLeaf();
    }

    ~CHAIN_OUR() {
        delete[] swtich_hg;
        delete[] counter;
        for (int i = 0; i < switch_num; ++i) {
            delete[] filter_cm[i];
        }
        delete[] filter_cm;
        for (int i = 0; i < path_num; ++i) {
            delete[] interval[i];
        }
        delete[] interval;
    }

    void config_FatTree() {

        unordered_map<int, vector<double>> dis;
        dis[1] = {1, 1};
        dis[3] = {1, 8, 1, 10};
        dis[5] = {11, 12, 14, 12, 11, 60};
        unordered_map<int, vector<double>> dis2;
        dis2[1] = {1, 1};
        dis2[3] = {1, 1, 1, 3};
        dis2[5] = {1, 1, 1, 1, 1, 5};

        for (int i = 0; i < path_num; ++i) {
            uint32_t cur1 = 0, cur2 = 0, maxv = MAXI;
            int length = path[i].size();
            vector<pii> tp1, tp2;
            vector<double> cur_dis = dis[length];
            vector<double> cur_dis2 = dis2[length];
            for (int j = 0; j < length; ++j) {
                uint32_t end1 = cur1 + (uint32_t)(cur_dis[j] * 1.0 / cur_dis.back() * maxv);
                uint32_t end2 = cur2 + (uint32_t)(cur_dis2[j] * 1.0 / cur_dis2.back() * maxv);
                tp1.push_back({cur1, end1});
                tp2.push_back({cur2, end2});
                cur1 = end1 + 1;
                cur2 = end2 + 1;
            }
            tp1.back().second = maxv;
            tp2.back().second = maxv;
            // tp2 = tp1;
            for (int j = 0; j < length; ++j) {
                interval[i][j] = tp1[j];
            }
        }
    }

    void config_SpineLeaf() {

        unordered_map<int, vector<double>> dis;
        dis[1] = {1, 1};
        dis[3] = {9, 10, 9, 28};
        unordered_map<int, vector<double>> dis2;
        dis2[1] = {1, 1};
        dis2[3] = {1, 1, 1, 3};

        for (int i = 0; i < path_num; ++i) {
            uint32_t cur1 = 0, cur2 = 0, maxv = MAXI;
            int length = path[i].size();
            vector<pii> tp1, tp2;
            vector<double> cur_dis = dis[length];
            vector<double> cur_dis2 = dis2[length];
            for (int j = 0; j < length; ++j) {
                uint32_t end1 = cur1 + (uint32_t)(cur_dis[j] * 1.0 / cur_dis.back() * maxv);
                uint32_t end2 = cur2 + (uint32_t)(cur_dis2[j] * 1.0 / cur_dis2.back() * maxv);
                tp1.push_back({cur1, end1});
                tp2.push_back({cur2, end2});
                cur1 = end1 + 1;
                cur2 = end2 + 1;
            }
            tp1.back().second = maxv;
            tp2.back().second = maxv;
            // tp2 = tp1;
            for (int j = 0; j < length; ++j) {
                interval[i][j] = tp1[j];
            }
        }
    }

    void show_config() {
        if (k != 0) {
            printf("FatTree: k: %d\n", k);
        } else {
            printf("SpineLeaf: n: %d, m: %d\n", n, m);
        }
        printf("totalMem: %f, row_num: %d, col_num: %d\n", row_num * col_num * (32 + 24.0) / 8 / 1024, row_num, col_num);
    }

    uint16_t QFP(uint32_t hash) {
        hash ^= hash >> 16;
        hash *= 0x85ebca6b;
        hash ^= hash >> 13;
        hash *= 0xc2b2ae35;
        hash ^= hash >> 16;
        return hash & 65535;
    }

    void update_filter_cm(Cell* cm, uint32_t sid, uint32_t f, uint32_t e) {
        uint32_t hash = H(f, seed[0]);
        uint16_t fp = FP(hash);
        uint32_t index = hash % filter_w;
        Cell* cell = &cm[index];
        if (cell->flow_id == 0) {
            cell->flow_id = fp;
            cell->size = 1;
        } else {
            if (cell->flow_id == fp) {
                if (cell->size < (1 << FILTER_COUNTER_SIZE)) {
                    cell->size += 1;
                }
            } else {
                cell->size -= 1;
                if (cell->size == 0) {
                    cell->flow_id = fp;
                    cell->size = 1;
                }
            }
        }
    }

    double query_filter_cm(Cell* cm, uint32_t f) {
        uint32_t hash = H(f, seed[0]);
        uint16_t fp = FP(hash);
        uint32_t index = hash % filter_w;
        Cell* cell = &cm[index];
        if (cell->flow_id == fp) {
            return cell->size / filter_ratio;
        }
        return 0;
    }

    void send(uint32_t f, uint32_t e) {
        int path_id = flow_path_id[f];
        vector<int> current_path = path[path_id];
        int flag = 0, start = current_path[0];
        int switch_id = 0;
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            if (flag == 0) {
                uint32_t hash_value = H(f, index_seed);
                if (interval[path_id][i].first <= hash_value && hash_value <= interval[path_id][i].second) {
                    switch_id = current_switch;
                    swtich_hg[current_switch].insert(f);
                    break;
                }
            }
        }
        packet_load[switch_id]++;
        flow_load[switch_id].insert(f);
    }

    uint32_t query(uint32_t f) {
        if (per_flow_size.find(f) != per_flow_size.end()) {
            return per_flow_size[f];
        }
        int path_id = flow_path_id[f];
        vector<int> current_path = path[path_id];
        uint32_t hash_value = H(f, index_seed);
        uint32_t s1 = 0;
        for (int i = 0; i < current_path.size(); i++) {
            if (interval[path_id][i].first <= hash_value && hash_value <= interval[path_id][i].second) {
                s1 = swtich_hg[current_path[i]].query(f);
            }
        }
        return s1;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int total) {
        unordered_map<uint32_t, uint32_t> res;
        // unordered_set<uint32_t> flows;
        for (int i = 0; i < switch_num; i++) {
            unordered_map<uint32_t, uint32_t> estimatedFlowSizes;
            swtich_hg[i].getEstimatedFlowSizes(estimatedFlowSizes);
            for (auto& it : estimatedFlowSizes) {
                per_flow_size[it.first] += it.second;
                if (per_flow_size[it.first] >= total) {
                    res[it.first] = per_flow_size[it.first];
                }
            }
        }
        return res;
    }
};

#endif