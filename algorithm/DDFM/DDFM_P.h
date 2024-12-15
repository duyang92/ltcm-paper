#ifndef _DDFM_P_
#define _DDFM_P_

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

class DDFM_P {
#define MAXI 4294967295

#define MAX_FILTER_COUNTER 100000
#define FILTER_W 1024 * 2
#define FILTER_D 1
#define FILTER_COUNTER_SIZE 4
#define FILTER_FINGERPRINT_SIZE 16

#define HH_L_THRESHOLD 1
    // #define HH_THRESHOLD 8000
    // #define RE_DIS_THREOLD 10000

#define BUCKET_NUM_L 8
#define BUCKET_NUM_S 8

    typedef pair<uint32_t, uint32_t> pii;

   public:
    struct Bucket_L {
        uint16_t fingerprints[BUCKET_NUM_L];
        uint32_t size[BUCKET_NUM_L];
        Bucket_L() {
            memset(fingerprints, 0, sizeof(fingerprints));
            memset(size, 0, sizeof(size));
        }
        ~Bucket_L() {}
    };

    struct Bucket_S {
        uint16_t fingerprints[BUCKET_NUM_S];
        uint32_t size[BUCKET_NUM_S];
        Bucket_S() {
            memset(fingerprints, 0, sizeof(fingerprints));
            memset(size, 0, sizeof(size));
        }
        ~Bucket_S() {}
    };

    struct Cell {
        uint16_t flow_id;
        uint8_t size;
        Cell() {
            flow_id = 0;
            size = 0;
        }
        ~Cell() {}
    };

   public:
    // filter sketch
    uint32_t reset;
    uint32_t insert_threshold, query_threshold;
    double filter_ratio;
    uint32_t filter_w, filter_d;
    double filter_memory;
    Cell** filter_cm;
    uint32_t* counter;

    // topology
    int k = 0;
    FatTree ft;
    int n, m;
    SpineLeaf sl;
    int switch_num, path_num;
    vector<vector<int>> path;

    // Bucket_L
    Bucket_L** buckets_L;
    uint32_t w_L;
    double buckets_memory_L;

    // Bucket_S
    Bucket_S** buckets_S;
    uint32_t w_S;
    double buckets_memory_S;

    // heavy hitter list
    uint32_t len_hh;
    unordered_set<uint32_t> offline_hh[21];
    uint32_t HH_THRESHOLD;
    uint32_t RE_DIS_THREOLD;
    unordered_set<uint32_t> re_dis_list;

    // flow distribution
    pii*** interval;

    // flow path
    unordered_map<uint32_t, int> flow_path_id;
    int* betweeness_centrality;

    // seed
    uint32_t* seed;
    uint32_t index_seed;

    // load
    uint32_t* packet_pass;
    double* process_time;
    vector<int> packet_load;
    vector<unordered_set<uint32_t>> flow_load;

    uint32_t l1, l2, h1, h2;

    DDFM_P(int k, double memory, uint32_t len, uint32_t* seeds, int hh_th = 8000, int re_dis_th = 10000) {
        srand((int)time(0));

        this->HH_THRESHOLD = hh_th;
        this->RE_DIS_THREOLD = re_dis_th;

        // fat tree config
        this->k = k;
        this->ft = FatTree(k);
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();

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

        // uint32_t wL = round(memory / 3 * 1024 * 8 / ((24 + 16) * BUCKET_NUM_L));
        uint32_t wL = round(memory * 2 / 9 * 1024 * 8 / ((24 + 16) * BUCKET_NUM_L));
        // uint32_t wS = round(memory * 2 / 3 * 1024 * 8 / ((16 + 16) * BUCKET_NUM_S));
        uint32_t wS = round(memory * 7 / 9 * 1024 * 8 / ((16 + 16) * BUCKET_NUM_S));

        // Bucket_L
        this->buckets_L = new Bucket_L*[switch_num];
        this->w_L = wL;
        for (int i = 0; i < switch_num; ++i) {
            buckets_L[i] = new Bucket_L[wL];
        }
        this->buckets_memory_L = (double)w_L * (24 + 16) * BUCKET_NUM_L / 8 / 1024;

        // Bucket_S
        this->buckets_S = new Bucket_S*[switch_num];
        this->w_S = wS;
        for (int i = 0; i < switch_num; ++i) {
            buckets_S[i] = new Bucket_S[wS];
        }
        this->buckets_memory_S = (double)w_S * (16 + 16) * BUCKET_NUM_S / 8 / 1024;

        // heavy hitter list
        this->len_hh = len;

        // flow distribution
        // this->original_interval = vector<vector<pii>>(path_num);
        // this->overflow_interval = vector<vector<pii>>(path_num);
        this->interval = new pii**[path_num];
        for (int i = 0; i < path_num; ++i) {
            interval[i] = new pii*[2];
            for (int j = 0; j < 2; ++j) {
                interval[i][j] = new pii[path[i].size()];
                memset(interval[i][j], 0, sizeof(pii) * path[i].size());
            }
        }

        // flow path
        this->betweeness_centrality = new int[switch_num];
        memset(betweeness_centrality, 0, sizeof(int) * switch_num);

        // seed
        this->seed = seeds;
        this->index_seed = 131;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config_FatTree();
    }

    DDFM_P(int n, int m, double memory, uint32_t len, uint32_t* seeds, int hh_th = 8000, int re_dis_th = 10000) {
        srand((int)time(0));

        this->HH_THRESHOLD = hh_th;
        this->RE_DIS_THREOLD = re_dis_th;

        // spine leaf config
        this->n = n;
        this->m = m;
        this->sl = SpineLeaf(n, m);
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();

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

        uint32_t wL = round(memory / 3 * 1024 * 8 / ((24 + 16) * BUCKET_NUM_L));
        uint32_t wS = round(memory * 2 / 3 * 1024 * 8 / ((16 + 16) * BUCKET_NUM_S));

        // Bucket_L
        this->buckets_L = new Bucket_L*[switch_num];
        this->w_L = wL;
        for (int i = 0; i < switch_num; ++i) {
            buckets_L[i] = new Bucket_L[wL];
        }
        this->buckets_memory_L = (double)w_L * (24 + 16) * BUCKET_NUM_L / 8 / 1024;

        // Bucket_S
        this->buckets_S = new Bucket_S*[switch_num];
        this->w_S = wS;
        for (int i = 0; i < switch_num; ++i) {
            buckets_S[i] = new Bucket_S[wS];
        }
        this->buckets_memory_S = (double)w_S * (16 + 16) * BUCKET_NUM_S / 8 / 1024;

        // heavy hitter list
        this->len_hh = len;

        // flow distribution
        this->interval = new pii**[path_num];
        for (int i = 0; i < path_num; ++i) {
            interval[i] = new pii*[2];
            for (int j = 0; j < 2; ++j) {
                interval[i][j] = new pii[path[i].size()];
                memset(interval[i][j], 0, sizeof(pii) * path[i].size());
            }
        }

        // flow path
        this->betweeness_centrality = new int[switch_num];
        memset(betweeness_centrality, 0, sizeof(int) * switch_num);

        // seed
        this->seed = seeds;
        this->index_seed = 131;

        // load
        this->packet_pass = new uint32_t[switch_num];
        memset(packet_pass, 0, sizeof(uint32_t) * switch_num);
        this->process_time = new double[switch_num];
        memset(process_time, 0, sizeof(double) * switch_num);

        config_SpineLeaf();
    }

    ~DDFM_P() {
        for (int i = 0; i < switch_num; ++i) {
            delete[] buckets_L[i];
            delete[] buckets_S[i];
            delete[] filter_cm[i];
        }
        delete[] buckets_L;
        delete[] buckets_S;
        delete[] filter_cm;
        delete[] counter;
        delete[] betweeness_centrality;
        for (int i = 0; i < path_num; ++i) {
            for (int j = 0; j < 2; ++j) {
                delete[] interval[i][j];
            }
            delete[] interval[i];
        }
        delete[] interval;
        delete[] packet_pass;
        delete[] process_time;
    }

    void config_FatTree() {
        for (int i = 0; i < switch_num; ++i) {
            betweeness_centrality[i] = ft.betweeness_centrality[i];
        }

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
                interval[i][0][j] = tp1[j];
                interval[i][1][j] = tp2[j];
            }
        }
    }

    void config_SpineLeaf() {
        for (int i = 0; i < switch_num; ++i) {
            betweeness_centrality[i] = sl.betwenness_centrality[i];
        }

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
                interval[i][0][j] = tp1[j];
                interval[i][1][j] = tp2[j];
            }
        }
    }

    void show_config() {
        if (k != 0) {
            printf("FatTree: k: %d\n", k);
        } else {
            printf("SpineLeaf: n: %d, m: %d\n", n, m);
        }
        printf("filter_w: %d, filter_d: %d, filter_memory: %.2fKB\n", filter_w, filter_d, filter_memory);
        printf("w_L: %d, buckets_memory_L: %.2fKB\n", w_L, buckets_memory_L);
        printf("w_S: %d, buckets_memory_S: %.2fKB\n", w_S, buckets_memory_S);
        printf("len_hh: %d, hh_list_size: %fKB\n", len_hh, len_hh * 4.0 / 1024);
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

    void update_buckets_L(int& sid, uint32_t& f, uint32_t& fingerprint, uint32_t& hash_value) {
        uint32_t index = hash_value % w_L;
        Bucket_L* cur_bucket = &buckets_L[sid][index];
        int flag = 1, minv = MAXI, mini = 0;
        for (int i = 0; i < BUCKET_NUM_L; i++) {
            uint16_t* fp = &cur_bucket->fingerprints[i];
            uint32_t* size = &cur_bucket->size[i];
            if (*fp == 0) {
                *fp = fingerprint;
                *size = 1;
                if (l1 <= hash_value && hash_value <= h1 && offline_hh[sid].size() < len_hh) {
                    offline_hh[sid].insert(f);
                }
                flag = 0;
                break;
            } else if (*fp == fingerprint) {
                *size += 1;
                flag = 0;
                break;
            } else if (*size < minv) {
                minv = *size;
                mini = i;
            }
        }
        if (flag) {
            if (l1 <= hash_value && hash_value <= h1 && offline_hh[sid].size() < len_hh) {
                offline_hh[sid].insert(f);
            }
            if ((double)rand() * (minv + 1) < RAND_MAX) {
                cur_bucket->fingerprints[mini] = fingerprint;
            }
        }
    }

    uint32_t query_buckets_L(Bucket_L* buckets, uint32_t& f, uint32_t& fingerprint, uint32_t& hash_value) {
        uint32_t index = hash_value % w_L;
        Bucket_L* cur_bucket = &buckets[index];
        for (int i = 0; i < BUCKET_NUM_L; i++) {
            if (cur_bucket->fingerprints[i] == fingerprint) {
                return cur_bucket->size[i];
            }
        }
        return 0;
    }

    void update_buckets_S(int& sid, uint32_t& f, uint32_t& fingerprint, uint32_t& hash_value) {
        uint32_t index = hash_value % w_S;
        Bucket_S* cur_bucket = &buckets_S[sid][index];
        int flag = 1, minv = MAXI, mini = 0;
        for (int i = 0; i < BUCKET_NUM_S; i++) {
            uint16_t* fp = &cur_bucket->fingerprints[i];
            uint32_t* size = &cur_bucket->size[i];
            if (*fp == 0) {
                *fp = fingerprint;
                *size = 1;
                flag = 0;
                break;
            } else if (*fp == fingerprint) {
                if (*size == 65535) {
                    update_buckets_L(sid, f, fingerprint, hash_value);
                    return;
                }
                if (*size == HH_THRESHOLD) {
                    if (offline_hh[sid].size() < len_hh) {
                        offline_hh[sid].insert(f);
                    }
                }
                *size += 1;
                flag = 0;
                break;
            } else if (*size < minv) {
                minv = *size;
                mini = i;
            }
        }
        if (flag) {
            if ((double)rand() * (minv + 1) < RAND_MAX) {
                cur_bucket->fingerprints[mini] = fingerprint;
            }
        }
    }

    uint32_t query_buckets_S(Bucket_S* buckets, uint32_t& f, uint32_t& fingerprint, uint32_t& hash_value) {
        uint32_t index = hash_value % w_S;
        Bucket_S* cur_bucket = &buckets[index];
        uint32_t size = 1e9;
        for (int i = 0; i < BUCKET_NUM_S; i++) {
            if (cur_bucket->fingerprints[i] == fingerprint) {
                return cur_bucket->size[i];
            }
            size = min(size, (uint32_t)cur_bucket->size[i]);
        }
        return 1;
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
        int flag = 0, start = current_path[0], switch_id = 0;
        uint32_t hash_tp = 0, fp_tp = 0, index_tp = 0;
        packet_pass[start] += 1;
#ifndef _WIN32
        double start_time = now();
#else
        clock_t start_time = clock();
#endif
        l1 = interval[path_id][0][0].first;
        h1 = interval[path_id][0][0].second;
        l2 = interval[path_id][1][0].first;
        h2 = interval[path_id][1][0].second;
        if (counter[start] < insert_threshold) {
            hash_tp = H(f, index_seed);
            fp_tp = FP(hash_tp);
            index_tp = hash_tp % filter_w;
            Cell* cell = &filter_cm[start][index_tp];
            if (cell->flow_id == 0) {
                cell->flow_id = fp_tp;
                cell->size = 1;
            } else {
                if (cell->flow_id == fp_tp) {
                    cell->size += 1;
                } else {
                    if ((cell->size -= 1) == 0) {
                        cell->flow_id = fp_tp;
                        cell->size = 1;
                    }
                }
            }
        }
        if (++counter[start] == reset) counter[start] = 0;
        if (hash_tp == 0) {
            hash_tp = H(f, index_seed);
            fp_tp = FP(hash_tp);
            index_tp = hash_tp % filter_w;
        }
        Cell* cell = &filter_cm[start][index_tp];
        if (cell->flow_id == fp_tp && cell->size >= filter_ratio * RE_DIS_THREOLD) {
            flag = 1;
        }
        if (flag == 0) {
            if (l1 <= hash_tp && hash_tp <= h1) {
                update_buckets_S(start, f, fp_tp, hash_tp);
                flag = 2;
            }
        } else {
            uint32_t hash_value = H(f ^ e, index_seed);
            if (l2 <= hash_value && hash_value <= h2) {
                update_buckets_L(start, f, fp_tp, hash_tp);
                flag = 2;
            }
        }
#ifndef _WIN32
        process_time[start] += now() - start_time;
#else
        process_time[start] += (double)(clock() - start_time) / CLOCKS_PER_SEC;
#endif
        for (int i = 1; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            packet_pass[current_switch] += 1;
            if (flag == 2) continue;
#ifndef _WIN32
            double start_time = now();
#else
            clock_t start_time = clock();
#endif
            l1 = interval[path_id][0][i].first;
            h1 = interval[path_id][0][i].second;
            l2 = interval[path_id][1][i].first;
            h2 = interval[path_id][1][i].second;
            if (flag == 0) {
                uint32_t hash_value = H(f, index_seed);
                if (l1 <= hash_value && hash_value <= h1) {
                    uint32_t fingerprint = FP(hash_value);
                    update_buckets_S(current_switch, f, fingerprint, hash_value);
                    flag = 2;
                }
            } else {
                uint32_t hash_value = H(f ^ e, index_seed);
                if (l2 <= hash_value && hash_value <= h2) {
                    uint32_t hash_v = H(f, index_seed);
                    uint32_t fingerprint = FP(hash_v);
                    update_buckets_L(current_switch, f, fingerprint, hash_v);
                    flag = 2;
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
        uint32_t hash_value = H(f, index_seed);
        uint32_t fingerprint = FP(hash_value);
        uint32_t s1 = 0, s2 = 0;
        uint32_t flag = re_dis_list.find(f) == re_dis_list.end() ? 0 : 1;
        for (int i = 0; i < current_path.size(); i++) {
            if (interval[path_id][0][i].first <= hash_value && hash_value <= interval[path_id][0][i].second) {
                s1 = query_buckets_S(buckets_S[current_path[i]], f, fingerprint, hash_value);
            }
            if (flag == 1) {
                uint32_t size2 = query_buckets_L(buckets_L[current_path[i]], f, fingerprint, hash_value);
                s2 += size2;
            }
        }
        return s1 + s2;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int threshold) {
        unordered_map<uint32_t, uint32_t> hh;
        for (int i = 1; i < switch_num; i++) {
            for (const uint32_t& f : offline_hh[i]) {
                re_dis_list.insert(f);
                uint32_t sz = query(f);
                if (sz >= threshold) {
                    hh[f] = max(hh[f], sz);
                }
            }
        }
        return hh;
    }
};

#endif