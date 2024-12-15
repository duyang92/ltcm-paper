#ifndef _MINHEAP_H
#define _MINHEAP_H
#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
using namespace std;

class MinHeap {
    typedef pair<uint32_t, uint32_t> node;
   public:
    vector<node> heap;       
    int CurrentSize;  
    int DefaultSize = 500;
    unordered_map<uint32_t, int> flow_idx;

    void ShiftDown(int start, int m) {
        int i = start, j = 2 * i + 1; 
        node temp = heap[i];
        while (j <= m) {
            if (j < m && heap[j].second > heap[j + 1].second)
                j++;  
            if (temp.second <= heap[j].second)
                break;
            else {
                heap[i] = heap[j];
                flow_idx[heap[j].first] = i;
                i = j;
                j = 2 * j + 1;
            }
        }
        heap[i] = temp;
        flow_idx[temp.first] = i;
    } 

    void ShiftUp(int start) {
        int j = start, i = (j - 1) / 2;
        node temp = heap[j];
        while (j > 0) {
            if (heap[i].second <= temp.second)
                break;
            else {
                heap[j] = heap[i];
                flow_idx[heap[i].first] = j;
                j = i;
                i = (j - 1) / 2;
            }
        }
        heap[j] = temp;
        flow_idx[temp.first] = j;
    }  

   public:
    MinHeap(int sz = 500) {
        heap = vector<node>(sz);
        DefaultSize = sz;
        CurrentSize = 0;
    }

    // ~MinHeap() {
    //     if (heap) delete[] heap;
    // }

    bool Insert(node x) {
        if (flow_idx.find(x.first) != flow_idx.end()) {
            int idx = flow_idx[x.first];
            if (x.second < heap[idx].second) return false;
            heap[idx] = x;
            ShiftDown(idx, CurrentSize - 1);
            return true;
        }
        if (CurrentSize == DefaultSize) {  
            return ReplaceMin(x);
        }
        heap[CurrentSize] = x;  
        flow_idx[x.first] = CurrentSize;
        ShiftUp(CurrentSize);   
        CurrentSize++;
        return true;
    } 

    bool ReplaceMin(node x) {
        if (!CurrentSize) {
            std::cout << "no size\n";
            return false;
        }
        if (x.second >= heap[0].second) {
            heap[0] = x;
            flow_idx[x.first] = 0;
            ShiftDown(0, CurrentSize - 1);  
            return true;
        }
        return false;  
    }

    bool empty() { return CurrentSize == 0; }
};

#endif