#ifndef _DS_P_
#define _DS_P_

#include <string.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../utils/FatTree.h"
#include "../../utils/MurmurHash3.h"
#include "../../utils/SpineLeaf.h"

using namespace std;

class DS_P {
#define MAXI 4294967295

    typedef pair<double, double> pdd;

   public:
    // count min sketch
    uint32_t cm_w, logical_w;
    uint32_t cm_d;
    double on_chip_size;
    // topology
    int k = 0;
    FatTree ft;
    int n, m;
    SpineLeaf sl;
    int switch_num, path_num;
    vector<vector<int>> path;
    uint32_t*** switch_cm;
    vector<vector<pdd>> switches_interval;
    // heavy hitter list
    uint32_t len_hh;
    unordered_set<uint32_t> offline_hh[21];
    uint32_t thresholdToSaveKeys;
    // flow path
    unordered_map<uint32_t, int> flow_path_id;
    int* betweeness_centrality;
    // load
    uint32_t* packet_pass;
    double* process_time;
    vector<int> packet_load;
    vector<unordered_set<uint32_t>> flow_load;
    // seed
    uint32_t* seed;

    DS_P(int k, double memory, uint32_t len, uint32_t* seeds, int hh_th = 8000) {
        srand((int)time(0));

        this->thresholdToSaveKeys = hh_th;

        this->cm_d = 2;
        this->cm_w = round(memory * 8 * 1024 / 24 / cm_d);
        this->logical_w = (1 << 16);
        this->on_chip_size = cm_w * cm_d * 24 / 8 / 1024;

        this->k = k;
        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();
        this->switch_cm = new uint32_t**[switch_num];
        this->switches_interval = vector<vector<pdd>>();
        this->betweeness_centrality = new int[switch_num];
        memset(betweeness_centrality, 0, sizeof(int) * switch_num);

        this->len_hh = len;

        this->seed = seeds;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config();
    }

    DS_P(int n, int m, double memory, uint32_t len, uint32_t* seeds, int hh_th = 8000) {
        srand((int)time(0));

        this->thresholdToSaveKeys = hh_th;

        this->cm_d = 2;
        this->cm_w = round(memory * 8 * 1024 / 24 / cm_d);
        this->logical_w = (1 << 16);
        this->on_chip_size = cm_w * cm_d * 24 / 8 / 1024;

        this->n = n;
        this->m = m;
        this->sl = SpineLeaf(n, m);
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();
        this->switch_cm = new uint32_t**[switch_num];
        this->switches_interval = vector<vector<pdd>>();
        this->betweeness_centrality = new int[switch_num];
        memset(betweeness_centrality, 0, sizeof(int) * switch_num);

        this->len_hh = len;

        this->seed = seeds;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config();
    }

    ~DS_P() {
        for (int i = 0; i < switch_num; ++i) {
            for (int j = 0; j < cm_d; ++j) {
                delete[] switch_cm[i][j];
            }
            delete[] switch_cm[i];
        }
        delete[] switch_cm;
        delete[] betweeness_centrality;
        delete[] packet_pass;
        delete[] process_time;
    }

    uint32_t H(uint32_t t, uint32_t s) {
        uint32_t hash_val = 0;
        char hash_input_str[5] = {0};
        memcpy(hash_input_str, &t, sizeof(uint32_t));
        MurmurHash3_x86_32(hash_input_str, 4, s, &hash_val);
        return hash_val;
    }

    void config() {
        for (int i = 0; i < switch_num; ++i) {
            switch_cm[i] = new uint32_t*[cm_d];
            for (int j = 0; j < cm_d; ++j) {
                switch_cm[i][j] = new uint32_t[cm_w];
                memset(switch_cm[i][j], 0, sizeof(uint32_t) * cm_w);
            }
        }

        for (int i = 0; i < switch_num; i++) {
            if (k != 0) {
                betweeness_centrality[i] = ft.betweeness_centrality[i];
            } else {
                betweeness_centrality[i] = sl.betwenness_centrality[i];
            }
        }

        for (int i = 0; i < path.size(); ++i) {
            int total = 0;
            int length = path[i].size();
            vector<pdd> tp;
            for (int j = 0; j < length; ++j) {
                total += cm_w / betweeness_centrality[path[i][j]];
            }
            double a = 0, b = 0;
            for (int j = 0; j < length; ++j) {
                if (j > 0) {
                    a += cm_w / betweeness_centrality[path[i][j - 1]];
                }
                b += cm_w / betweeness_centrality[path[i][j]];
                tp.push_back({a / total, b / total});
            }
            switches_interval.push_back(tp);
        }
    }

    void show_config() {
        if (k != 0) {
            printf("FatTree: k: %d\n", k);
        } else {
            printf("SpineLeaf: n: %d, m: %d\n", n, m);
        }
        printf("cm_w: %d, cm_d: %d, on_chip_size: %.2f KB\n", cm_w, cm_d, on_chip_size);
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
        for (int i = 0; i < current_path.size(); ++i) {
            int current_switch = current_path[i];
            packet_pass[current_switch]++;
#ifndef _WIN32
            double start_time = now();
#else
            clock_t start_time = clock();
#endif
            uint32_t** current_cm = switch_cm[current_switch];
            int alpha = (int)round(switches_interval[path_id][i].first * logical_w);
            int beta = (int)round(switches_interval[path_id][i].second * logical_w);
            int process_flag = 0, value = 1e9;
            for (int j = 0; j < cm_d; ++j) {
                int hashV = H(f, seed[j]) % logical_w + 1;
                if (alpha + 1 <= hashV && hashV <= beta) {
                    int index = (int)round((double)(hashV - alpha - 1) / (beta - alpha - 1) * (cm_w - 1));
                    current_cm[j][index] += 1;
                    value = min(value, (int)current_cm[j][index]);
                    process_flag = 1;
                }
            }
            if (value == thresholdToSaveKeys && process_flag == 1) {
                if (offline_hh[current_switch].size() < len_hh) {
                    offline_hh[current_switch].insert(f);
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
        int path_id = flow_path_id[f];
        vector<int> current_path = path[path_id];
        uint32_t res = MAXI;
        for (int i = 0; i < current_path.size(); ++i) {
            uint32_t** current_cm = switch_cm[current_path[i]];
            int alpha = (int)round(switches_interval[path_id][i].first * logical_w);
            int beta = (int)round(switches_interval[path_id][i].second * logical_w);
            for (int j = 0; j < cm_d; ++j) {
                int hashV = H(f, seed[j]) % logical_w + 1;
                if (alpha + 1 <= hashV && hashV <= beta) {
                    int index = (int)round((double)(hashV - alpha - 1) / (beta - alpha - 1) * (cm_w - 1));
                    res = min(res, current_cm[j][index]);
                }
            }
        }
        if (res < 1 or res == MAXI) res = 1;
        return res;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int threshold) {
        unordered_map<uint32_t, uint32_t> hh;
        for (int i = 1; i < switch_num; i++) {
            for (auto f : offline_hh[i]) {
                uint32_t sz = query(f);
                if (sz >= threshold) {
                    hh[f] = sz;
                }
            }
        }
        return hh;
    }
};

#endif