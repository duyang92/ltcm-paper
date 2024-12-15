#ifndef _CMAX_P_
#define _CMAX_P_

#include <string.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../utils/FatTree.h"
#include "../../utils/MurmurHash3.h"
#include "../../utils/SpineLeaf.h"

using namespace std;

class CMAX_P {
#define MAXI 4294967295

    typedef pair<uint32_t, uint32_t> pii;

   public:
    // count min sketch
    uint32_t cm_w;
    uint32_t cm_d;
    double on_chip_size;

    // topology
    int k = 0;
    FatTree ft;
    int n, m;
    SpineLeaf sl;
    int switch_num, path_num;
    vector<vector<int>> path;
    pii*** switch_cm;
    bool* enable;

    // flow path
    unordered_map<uint32_t, int> flow_path_id;

    // seed
    uint32_t* seed;

    // load
    uint32_t* packet_pass;
    double* process_time;

    CMAX_P(int k, double memory, uint32_t* seeds) {
        srand((int)time(0));

        this->cm_w = round(memory * 8 * 1024 / (32 + 24) / 3);
        this->cm_d = 3;
        this->on_chip_size = cm_w * cm_d * (32 + 24) / 8 / 1024;

        this->k = k;
        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();
        this->switch_cm = new pii**[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            switch_cm[i] = new pii*[cm_d];
            for (int j = 0; j < cm_d; ++j) {
                switch_cm[i][j] = new pii[cm_w];
                memset(switch_cm[i][j], 0, sizeof(pii) * cm_w);
            }
        }
        this->enable = new bool[switch_num];
        memset(enable, 0, sizeof(bool) * switch_num);

        this->seed = seeds;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config();
    }

    CMAX_P(int n, int m, double memory, uint32_t* seeds) {
        srand((int)time(0));

        this->cm_w = round(memory * 8 * 1024 / (32 + 24) / 3);
        this->cm_d = 3;
        this->on_chip_size = cm_w * cm_d * (32 + 24) / 8 / 1024;

        this->n = n;
        this->m = m;
        this->sl = SpineLeaf(n, m);
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();
        this->switch_cm = new pii**[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            switch_cm[i] = new pii*[cm_d];
            for (int j = 0; j < cm_d; ++j) {
                switch_cm[i][j] = new pii[cm_w];
                memset(switch_cm[i][j], 0, sizeof(pii) * cm_w);
            }
        }
        this->enable = new bool[switch_num];
        memset(enable, 0, sizeof(bool) * switch_num);

        this->seed = seeds;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config();
    }

    ~CMAX_P() {
        for (int i = 0; i < switch_num; ++i) {
            for (int j = 0; j < cm_d; ++j) {
                delete[] switch_cm[i][j];
            }
            delete[] switch_cm[i];
        }
        delete[] switch_cm;
        delete[] enable;
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
        if (k != 0) {
            for (int i = ft.iCore + ft.iAgg + 1; i < ft.switches_num; ++i) {
                enable[i] = 1;
            }
        } else {
            for (int i = sl.iSpine + 1; i < sl.switches_num; ++i) {
                enable[i] = 1;
            }
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

    void update_cm(pii** cm, uint32_t f, uint32_t e) {
        for (int i = 0; i < cm_d; i++) {
            uint32_t j = H(f, seed[i]) % cm_w;
            if (cm[i][j].first == 0) {
                cm[i][j].first = f;
                cm[i][j].second = 1;
            } else {
                if (cm[i][j].first == f) {
                    cm[i][j].second += 1;
                } else {
                    int tmp = (cm[i][j].second -= 1);
                    if (tmp <= 0) {
                        cm[i][j].first = f;
                        cm[i][j].second = 0;
                    }
                }
            }
        }
    }

    uint32_t query_cm(pii** cm, uint32_t f) {
        uint32_t size = 1;
        for (uint32_t i = 0; i < cm_d; i++) {
            uint32_t j = H(f, seed[i]) % cm_w;
            if (cm[i][j].first == f) {
                size = max(size, cm[i][j].second);
            }
        }
        return size;
    }

#ifndef _WIN32
    double now() {
        struct timespec tv;
        clock_gettime(CLOCK_MONOTONIC, &tv);
        return (tv.tv_sec + (double)tv.tv_nsec / 1000000000.0);
    }
#endif

    void send(uint32_t f, uint32_t e) {
        uint32_t id = flow_path_id[f];
        vector<int> current_path = path[id];
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            if (enable[current_switch] == 0) continue;
            packet_pass[current_switch] += 1;
#ifndef _WIN32
            double start_time = now();
#else
            clock_t start_time = clock();
#endif
            // uint32_t hashV = H(f ^ e, 31);
            uint32_t hashV = H(f, 31);
            if ((hashV % 2 == 0 && i == 0) or (hashV % 2 == 1 && i == current_path.size() - 1)) {
                update_cm(switch_cm[current_switch], f, e);
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
        uint32_t res = 0;
        uint32_t hashV = H(f, 31);
        for (int i = 0; i < current_path.size(); i++) {
            if (enable[current_path[i]] == 0) {
                continue;
            }
            if ((hashV % 2 == 0 && i == 0) or (hashV % 2 == 1 && i == current_path.size() - 1)) {
                uint32_t size1 = query_cm(switch_cm[current_path[i]], f);
                res = size1;
            }
        }
        if (res < 1) res = 1;
        return res;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int total) {
        unordered_map<uint32_t, uint32_t> res;
        for (int i = 0; i < switch_num; i++) {
            if (enable[i] == 0) continue;
            for (int j = 0; j < cm_d; j++) {
                for (int k = 0; k < cm_w; k++) {
                    uint32_t f = switch_cm[i][j][k].first;
                    int ss = query(f);
                    if (ss >= total) {
                        res[f] = ss;
                    }
                }
            }
        }
        return res;
    }
};

#endif