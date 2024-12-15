#ifndef _WS_OUR_P_
#define _WS_OUR_P_

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

class WS_OUR_P {
#define MAXI 4294967295

#define MAX_FILTER_COUNTER 100000
#define FILTER_W 1024 * 2
#define FILTER_D 1
#define FILTER_COUNTER_SIZE 4
#define FILTER_FINGERPRINT_SIZE 16

    // #define RE_DIS_THREOLD 10000

#define landa_d 16

#define KEY_SIZE 4

    typedef pair<uint32_t, uint32_t> pii;

   public:
    const int COUNT[2] = {1, -1};
    const int factor = 1;

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
        const int COUNT[2] = {1, -1};
        const int factor = 1;
        uint32_t items[landa_d];
        uint32_t counters[landa_d];
        bitset<landa_d> real;
        uint32_t incast;

        hg_node()
            : items{0}, counters{0} {
            incast = 0;
            real.reset();
        }
        void insert(uint32_t f) {
            uint32_t choice = f & 1;
            uint32_t min_num = UINT32_MAX;
            uint32_t min_pos = -1;

            for (uint32_t i = 0; i < landa_d; ++i) {
                if (counters[i] == 0) {
                    items[i] = f;
                    counters[i] = 1;
                    real[i] = 1;
                    return;
                } else if (items[i] == f) {
                    if (real[i])
                        counters[i]++;
                    else {
                        counters[i]++;
                        incast += COUNT[choice];
                    }
                    return;
                }

                if (counters[i] < min_num) {
                    min_num = counters[i];
                    min_pos = i;
                }
            }

            if (incast * COUNT[choice] >= int(min_num * factor)) {
                // count_type pre_incast = incast;
                if (real[min_pos]) {
                    uint32_t min_choice = items[min_pos] & 1;
                    incast += COUNT[min_choice] * counters[min_pos];
                }
                items[min_pos] = f;
                counters[min_pos] += 1;
                real[min_pos] = 0;
            }
            incast += COUNT[choice];
        }

        uint32_t query(uint32_t f) {
            for (uint32_t i = 0; i < landa_d; ++i) {
                if (items[i] == f) {
                    return counters[i];
                }
            }
            return max(1U, incast * COUNT[f & 1]);
        }
    };

    // memory config
    uint32_t totalMem;
    uint32_t bucketMem;
    uint32_t bucketNum;

    // topology
    int k = 0;
    FatTree ft;
    int n, m;
    SpineLeaf sl;
    int switch_num, path_num;
    vector<vector<int>> path;
    hg_node** swtich_hg;
    unordered_map<uint32_t, uint32_t> per_flow_size;

    // filter sketch
    uint32_t reset;
    uint32_t insert_threshold, query_threshold;
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
    uint32_t* packet_pass;
    double* process_time;

    WS_OUR_P(int k, double total_memory, uint32_t* seeds, int re_dis_th = 10000) {
        srand((int)time(0));

        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();
        this->k = k;

        // filter sketch
        this->reset = MAX_FILTER_COUNTER;
        this->insert_threshold = this->reset * 0.001;
        this->query_threshold = this->reset * 0.5;
        this->filter_ratio = 0.001;
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
        this->bucketMem = totalMem;
        this->bucketNum = bucketMem / ((32 + 24) * landa_d + landa_d + 24);  // this->number of buckets

        this->swtich_hg = new hg_node*[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            this->swtich_hg[i] = new hg_node[bucketNum];
            memset(swtich_hg[i], 0, sizeof(hg_node) * bucketNum);
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
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config_FatTree();
    }

    WS_OUR_P(int n, int m, double total_memory, uint32_t* seeds, int re_dis_th = 10000) {
        srand((int)time(0));

        this->n = n;
        this->m = m;
        this->sl = SpineLeaf(n, m);
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();

        // filter sketch
        this->reset = MAX_FILTER_COUNTER;
        this->insert_threshold = this->reset * 0.001;
        this->filter_ratio = 0.001;
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
        this->bucketMem = totalMem;
        this->bucketNum = bucketMem / ((32 + 24) * landa_d + landa_d + 24);  // this->number of buckets

        this->swtich_hg = new hg_node*[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            this->swtich_hg[i] = new hg_node[bucketNum];
            memset(swtich_hg[i], 0, sizeof(hg_node) * bucketNum);
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
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config_SpineLeaf();
    }

    ~WS_OUR_P() {
        for (int i = 0; i < switch_num; ++i) {
            delete[] swtich_hg[i];
        }
        delete[] swtich_hg;
        for (int i = 0; i < switch_num; ++i) {
            delete[] filter_cm[i];
        }
        delete[] filter_cm;
        delete[] counter;
        for (int i = 0; i < path_num; ++i) {
            delete[] interval[i];
        }
        delete[] interval;
    }

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
        printf("Total Mem: %fKB, Bucket Num: %d\n", 1.0 * totalMem / 1024 / 8, bucketNum);
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

#ifndef _WIN32
    double now() {
        struct timespec tv;
        clock_gettime(CLOCK_MONOTONIC, &tv);
        return (tv.tv_sec + (double)tv.tv_nsec / 1000000000.0);
    }
#endif

    void send(uint32_t f, uint32_t e) {
        int path_id = flow_path_id[f];
        vector<int> current_path = path[path_id];
        int flag = 0, start = current_path[0];
        int switch_id = 0;
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            packet_pass[current_switch]++;
#ifndef _WIN32
            double start_time = now();
#else
            clock_t start_time = clock();
#endif
            if (flag == 0) {
                uint32_t hash_value = H(f, index_seed);
                if (interval[path_id][i].first <= hash_value && hash_value <= interval[path_id][i].second) {
                    switch_id = current_switch;
                    uint32_t hash = hash_value;
                    swtich_hg[current_switch][hash % bucketNum].insert(f);
                    flag = 1;
                }
            }
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
        uint32_t hash_value = H(f, index_seed);
        uint32_t s1 = 0;
        for (int i = 0; i < current_path.size(); i++) {
            if (interval[path_id][i].first <= hash_value && hash_value <= interval[path_id][i].second) {
                s1 = swtich_hg[current_path[i]][hash_value % bucketNum].query(f);
            }
        }
        return s1;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int total) {
        unordered_map<uint32_t, uint32_t> res;
        // unordered_set<uint32_t> flows;
        for (int i = 0; i < switch_num; i++) {
            hg_node* hg = swtich_hg[i];
            for (int j = 0; j < bucketNum; j++) {
                for (int k = 0; k < landa_d; k++) {
                    per_flow_size[hg[j].items[k]] += hg[j].counters[k];
                    if (per_flow_size[hg[j].items[k]] >= total) {
                        res[hg[j].items[k]] = per_flow_size[hg[j].items[k]];
                    }
                }
            }
        }
        return res;
    }
};

#endif