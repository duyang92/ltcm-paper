#ifndef _WS_
#define _WS_

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

class WS {
#define MAXI 4294967295

#define landa_d 16

#define KEY_SIZE 4

    typedef pair<uint32_t, uint32_t> pii;

   public:
    const int COUNT[2] = {1, -1};
    const int factor = 1;

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

    // flow distribution
    int data_id;
    double* prob;

    // flow path
    unordered_map<uint32_t, int> flow_path_id;

    // seed
    uint32_t* seed;
    uint32_t index_seed;

    // load
    vector<int> packet_load;
    vector<unordered_set<uint32_t>> flow_load;

    WS(int k, double total_memory, uint32_t* seeds, int data_id = 0) {
        srand((int)time(0));

        this->k = k;
        this->ft = FatTree(k);  
        this->path = ft.path;
        this->switch_num = ft.switches_num;
        this->path_num = path.size();

        this->totalMem = total_memory * 8 * 1024;  // the total number of bits
        this->bucketMem = totalMem;
        this->bucketNum = bucketMem / ((32 + 24) * landa_d + landa_d + 24);  // this->number of buckets

        this->swtich_hg = new hg_node*[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            this->swtich_hg[i] = new hg_node[bucketNum];
            memset(swtich_hg[i], 0, sizeof(hg_node) * bucketNum);
        }

        // flow distribution
        this->data_id = data_id;
        this->prob = new double[switch_num];

        this->seed = seeds;
        this->index_seed = 13331;

        // load
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config();
    }

    WS(int n, int m, double total_memory, uint32_t* seeds, int data_id = 0) {
        srand((int)time(0));

        this->sl = SpineLeaf(n, m);
        this->n = n;
        this->m = m;
        this->path = sl.path;
        this->switch_num = sl.switches_num;
        this->path_num = path.size();

        this->totalMem = total_memory * 8 * 1024;  // the total number of bits
        this->bucketMem = totalMem;
        this->bucketNum = bucketMem / ((32 + 24) * landa_d + landa_d + 24);  // this->number of buckets

        this->swtich_hg = new hg_node*[switch_num];
        for (int i = 0; i < switch_num; ++i) {
            this->swtich_hg[i] = new hg_node[bucketNum];
            memset(swtich_hg[i], 0, sizeof(hg_node) * bucketNum);
        }

        // flow distribution
        this->data_id = data_id;
        this->prob = new double[switch_num];

        this->seed = seeds;
        this->index_seed = 13331;

        // load
        this->packet_load = vector<int>(switch_num, 0);
        this->flow_load = vector<unordered_set<uint32_t>>(switch_num);

        config();
    }

    ~WS() {
        for (int i = 0; i < switch_num; ++i) {
            delete[] swtich_hg[i];
        }
        delete[] swtich_hg;
        delete[] prob;
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

    void config() {
        double caida_ft[] = {0.0, 0.5434393679375896, 0.5434393679375896, 0.5446132883575244, 0.5446132883575244, 0.41391063851283083, 0.4144392600179558, 0.4146621036051138, 0.414984996011598, 0.41432661858680603, 0.41470145382308654, 0.4146358742749429, 0.415019152236538, 0.003142603863073307, 0.0033564768157314812, 0.005160061181215711, 0.0012958902363513676, 0.002707113134044416, 0.002755518993977512, 0.002539109281614459, 0.0};
        double dc_ft[] = {0.0, 0.49559492988133763, 0.49559492988133763, 0.49553689196518824, 0.49553689196518824, 0.4397959462767478, 0.4397706869272814, 0.43949783046909696, 0.4394684014463482, 0.4403113888383059, 0.44028185076572174, 0.43950023301230934, 0.4394605944079356, 0.0011775243326903773, 0.0006489529735050958, 0.0006284612946489328, 0.0008474907237139485, 0.000882239029778619, 0.001694656028177642, 0.0, 0.00029736625161362243};
        double caida_sl[] = {0.0, 0.33398748946864315, 0.33398748946864315, 0.33398748946864315, 0.33398748946864315, 0.0, 0.0002066434655861768, 0.0012634899230677917, 0.0037997888534383553, 0.002920922789378278, 0.003805940502429052, 0.0018866372723737774, 0.0018075666734540643};
        double dc_sl[] = {0.0, 0.33313147229577683, 0.33344340576441484, 0.3334012969384339, 0.3333573325703338, 0.0015699428926792779, 0.0008011638732034876, 0.0, 0.0016151013019250896, 0.000788146298229774, 0.0005671609314951893, 0.0015323257529476263, 0.0008177297828085974};
        for (int i = 0; i < switch_num; i++) {
            if (data_id == 0) {
                if (k > 0) {
                    prob[i] = caida_ft[i];
                } else {
                    prob[i] = caida_sl[i];
                }
            } else {
                if (k > 0) {
                    prob[i] = dc_ft[i];
                } else {
                    prob[i] = dc_sl[i];
                }
            }
        }
    }

    void show_config() {
        if (k > 0) {
            printf("FatTree: k=%d\n", k);
        } else {
            printf("SpineLeaf: n=%d, m=%d\n", n, m);
        }
        printf("Total Mem: %fKB, Bucket Num: %d\n", 1.0 * totalMem / 1024 / 8, bucketNum);
    }

    void send(uint32_t f, uint32_t e) {
        int path_id = flow_path_id[f], flag = 0;
        vector<int> current_path = path[path_id];
        uint32_t hash = H(f, index_seed);
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            if (flag == 1) continue;
            double hashV = H(f, 31);
            hashV /= (double)MAXI;
            if (i != current_path.size() - 1 && hashV >= prob[current_switch]) continue;
            swtich_hg[current_switch][hash % bucketNum].insert(f);
            packet_load[current_switch]++;
            flow_load[current_switch].insert(f);
            flag = 1;
        }
    }

    uint32_t query(uint32_t f) {
        if (per_flow_size.find(f) != per_flow_size.end()) {
            return per_flow_size[f];
        }
        int path_id = flow_path_id[f];
        vector<int> current_path = path[path_id];
        uint32_t res = 0;
        uint32_t hash = H(f, index_seed);
        for (int i = 0; i < current_path.size(); i++) {
            int current_switch = current_path[i];
            double hashV = H(f, 31);
            hashV /= (double)MAXI;
            if (i != current_path.size() - 1 && hashV >= prob[current_switch]) continue;
            res = swtich_hg[current_switch][hash % bucketNum].query(f);
        }
        if (res < 1) res = 1;
        return res;
    }

    unordered_map<uint32_t, uint32_t> query_heavy_hitters(int total) {
        unordered_map<uint32_t, uint32_t> res;
        // unordered_set<uint32_t> flows;
        for (int i = 0; i < switch_num; i++) {
            hg_node* hg = swtich_hg[i];
            for (int j = 0; j < bucketNum; j++) {
                for (int k = 0; k < landa_d; k++) {
                    per_flow_size[hg[j].items[k]] = hg[j].counters[k];
                    if (hg[j].counters[k] >= total) {
                        res[hg[j].items[k]] = hg[j].counters[k];
                    }
                }
            }
        }
        return res;
    }
};

#endif