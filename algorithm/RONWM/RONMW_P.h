#ifndef _RONMW_P_
#define _RONMW_P_

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
#include "NWPS.h"
#include "UWRA.h"

using namespace std;

template <typename T>
class RONMW_P {
#define MAXI 4294967295

    typedef pair<uint32_t, uint32_t> pii;

   public:
    // topology
    int k = 0;
    FatTree ft;
    int n, m;
    SpineLeaf sl;
    int switch_num, path_num;
    vector<vector<int>> path;
    int capacity, rt;
    double on_chip_size;
    vector<T> store;

    // flow path
    unordered_map<uint32_t, int> flow_path_id;

    // result
    unordered_map<uint32_t, int> flow_count;

    // load
    uint32_t* packet_pass;
    double* process_time;
    vector<int> packet_load;
    vector<unordered_set<uint32_t>> flow_load;

    RONMW_P(int k, double memory) {
        srand((int)time(0));

        this->capacity = round(memory * 8 * 1024 / (32 + 24 + 32));

        this->k = k;
        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();
        this->on_chip_size = capacity * (32 + 24 + 32) / 8 / 1024;
        this->flow_path_id = unordered_map<uint32_t, int>();

        this->flow_count = unordered_map<uint32_t, int>();

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config();
    }

    RONMW_P(int n, int m, double memory) {
        srand((int)time(0));

        this->capacity = round(memory * 8 * 1024 / (32 + 24 + 32));

        this->n = n;
        this->m = m;
        this->sl = SpineLeaf(n, m);
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();
        this->on_chip_size = capacity * (32 + 24 + 32) / 8 / 1024;
        this->flow_path_id = unordered_map<uint32_t, int>();

        this->flow_count = unordered_map<uint32_t, int>();

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config();
    }

    ~RONMW_P() {
        delete[] packet_pass;
        delete[] process_time;
    }

    void config() {
        store = vector<T>(switch_num, T(0, capacity));
        for (int i = 0; i < switch_num; i++) {
            store[i] = T(i, capacity);
        }
    }

    void show_config() {
        if (k != 0) {
            printf("FatTree: k: %d\n", k);
        } else {
            printf("SpineLeaf: n: %d, m: %d\n", n, m);
        }
        printf("Entry num: %d, size: %.2fKB\n", capacity, on_chip_size);
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
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            packet_pass[current_switch]++;
#ifndef _WIN32
            double start_time = now();
#else
            clock_t start_time = clock();
#endif
            store[current_switch].addFlow(f, e);
#ifndef _WIN32
            process_time[current_switch] += now() - start_time;
#else
            process_time[current_switch] += (double)(clock() - start_time) / CLOCKS_PER_SEC;
#endif
        }
    }

    int query(uint32_t f) {
        if (flow_count.size() == 0) {
            flow_count = store[0].getFlowCount(store);
        }
        int res = flow_count[f];
        if (res < 1) res = 1;
        return res;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int total) {
        flow_count = store[0].getFlowCount(store);
        unordered_map<uint32_t, uint32_t> res;
        for (auto [f, c] : flow_count) {
            if (c >= total) res[f] = c;
        }
        return res;
    }
};

#endif