#ifndef _NWPS_H_
#define _NWPS_H_

#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <vector>

#include "../../utils/MurmurHash3.h"

using namespace std;

class NWPS {
    typedef pair<double, pair<uint32_t, uint32_t>> pdii;

   public:
    class cmp {
       public:
        bool operator()(pdii& a, pdii& b) {
            return a.first > b.first;
        }
    };
    int sid;
    int capacity;
    priority_queue<pdii, vector<pdii>, cmp> pq;

    NWPS(int id, int cap) {
        this->sid = id;
        this->capacity = cap;
    }

    void addFlow(uint32_t flowId, uint32_t packetId) {
        uint64_t u_id = (((uint64_t)flowId) << 32) | packetId;
        uint32_t hashV = 0, seed = 31;
        char key[9] = {0};
        memcpy(key, &u_id, sizeof(uint64_t));
        MurmurHash3_x86_32(key, 8, seed, &hashV);
        double rands = (double)hashV / (double)0xFFFFFFFF;

        while (rands == 0) {
            seed++;
            MurmurHash3_x86_32(key, 8, seed, &hashV);
            rands = (double)hashV / (double)0xFFFFFFFF;
        }

        double p = 1.0 / rands; 

        if (pq.size() < capacity || p > pq.top().first) {
            if (pq.size() == capacity) {
                pq.pop();
            }
            pq.push({p, {flowId, packetId}});
        }
    }

    unordered_map<uint32_t, int> getFlowCount(vector<NWPS>& store) {
        set<double> f_s;

        priority_queue<pdii, vector<pdii>, cmp> temp_pq;
        for (int i = 0; i < store.size(); ++i) {
            while (store[i].pq.size() > 0) {
                auto [p, lst] = store[i].pq.top();
                store[i].pq.pop();
                if (f_s.find(p) == f_s.end()) {
                    f_s.insert(p);
                    temp_pq.push({p, lst});
                }
            }
        }
        double T = temp_pq.top().first;

        unordered_map<uint32_t, double> flow_count;
        while (temp_pq.size() > 0) {
            auto [p, lst] = temp_pq.top();
            temp_pq.pop();
            uint32_t f = lst.first;
            flow_count[f] += max(1.0, T);
        }
        unordered_map<uint32_t, int> res;
        for (auto [f, c] : flow_count) {
            res[f] = (int)round(c);
        }
        return res;
    }
};

#endif