#ifndef _CHAIN_OP_
#define _CHAIN_OP_

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../utils/FatTree.h"
#include "../../utils/MurmurHash3.h"
#include "../../utils/SpineLeaf.h"

class CHAIN_OP {
#define MAXI 4294967295

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

    struct Bucket {
        uint32_t key = 0;
        uint32_t C = 0;
    };

    class hg_node {
       public:
        Bucket** bucketArray;
        uint32_t ROW_NUM;
        uint32_t COL_NUM;
        uint32_t MAX_CHAIN_OP_LENGTH = 3;
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
                    while (l <= MAX_CHAIN_OP_LENGTH && tempSet.find(p2) == tempSet.end() && min_counter > bucketArray[min_rowIdx][p2].C) {
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

    // flow distribution
    int data_id;
    double* prob;

    // flow path
    unordered_map<uint32_t, int> flow_path_id;

    // seed
    uint32_t* seed;
    uint32_t index_seed;

    // load
    uint32_t* packet_pass;
    double* process_time;

    CHAIN_OP(int k, double total_memory, uint32_t* seeds, int data_id = 0) {
        srand((int)time(0));

        this->k = k;
        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();

        this->totalMem = total_memory * 8 * 1024;  // the total number of bits
        this->row_num = 4;
        this->col_num = 1.0 * totalMem / row_num / (32 + 24);

        this->swtich_hg = new hg_node[switch_num];
        for (uint32_t i = 0; i < switch_num; ++i) {
            swtich_hg[i] = hg_node(row_num, col_num);
        }

        // flow distribution
        this->data_id = data_id;
        this->prob = new double[switch_num];

        this->seed = seeds;
        this->index_seed = 13331;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config();
    }

    CHAIN_OP(int n, int m, double total_memory, uint32_t* seeds, int data_id = 0) {
        srand((int)time(0));

        this->sl = SpineLeaf(n, m);
        this->n = n;
        this->m = m;
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();

        this->totalMem = total_memory * 8 * 1024;  // the total number of bits
        this->row_num = 4;
        this->col_num = 1.0 * totalMem / row_num / (32 + 24);

        this->swtich_hg = new hg_node[switch_num];
        for (uint32_t i = 0; i < switch_num; ++i) {
            swtich_hg[i] = hg_node(row_num, col_num);
        }

        // flow distribution
        this->data_id = data_id;
        this->prob = new double[switch_num];

        this->seed = seeds;
        this->index_seed = 13331;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config();
    }

    ~CHAIN_OP() {
        delete[] swtich_hg;
        delete[] prob;
        delete[] packet_pass;
        delete[] process_time;
    }

    void config() {
    }

    void show_config() {
        if (k > 0) {
            printf("FatTree: k=%d\n", k);
        } else {
            printf("SpineLeaf: n=%d, m=%d\n", n, m);
        }
        printf("totalMem: %f, row_num: %d, col_num: %d\n", row_num * col_num * (32 + 24.0) / 8 / 1024, row_num, col_num);
    }

#ifndef _WIN32
    double now() {
        struct timespec tv;
        clock_gettime(CLOCK_MONOTONIC, &tv);
        return (tv.tv_sec + (double)tv.tv_nsec / 1000000000.0);
    }
#endif

    void send(uint32_t f, uint32_t e) {
        int path_id = flow_path_id[f], flag = 0;
        int switch_id = 0;
        vector<int> current_path = path[path_id];
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            packet_pass[current_switch]++;
#ifndef _WIN32
            double start_time = now();
#else
            clock_t start_time = clock();
#endif
            swtich_hg[current_switch].insert(f);
#ifndef _WIN32
            process_time[current_switch] += now() - start_time;
#else
            process_time[current_switch] += (double)(clock() - start_time) / CLOCKS_PER_SEC;
#endif
        }
    }

    uint32_t query(uint32_t f) {
        if (per_flow_size.find(f) != per_flow_size.end()) {
            return per_flow_size[f];
        }
        int path_id = flow_path_id[f];
        vector<int> current_path = path[path_id];
        uint32_t res = 1e9;
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            res = min(res, swtich_hg[current_switch].query(f));
        }
        if (res < 1) res = 1;
        return res;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int total) {
        unordered_map<uint32_t, uint32_t> res;
        // unordered_set<uint32_t> flows;
        for (int i = 0; i < switch_num; i++) {
            unordered_map<uint32_t, uint32_t> estimatedFlowSizes;
            swtich_hg[i].getEstimatedFlowSizes(estimatedFlowSizes);
            for (auto& it : estimatedFlowSizes) {
                per_flow_size[it.first] = it.second;
                if (it.second >= total) {
                    res[it.first] = it.second;
                }
            }
        }
        return res;
    }
};

#endif