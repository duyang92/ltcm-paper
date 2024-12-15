#ifdef _WIN32
#include <direct.h> // For _mkdir on Windows
#else
#include <sys/stat.h> // For mkdir on POSIX systems
#include <sys/types.h>
#endif
#include <string.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithm/CMAX/CMAX.h"
#include "algorithm/ChainSketch/CHAIN_O.h"
#include "algorithm/ChainSketch/CHAIN.h"
#include "algorithm/ChainSketch/CHAIN_OUR.h"
#include "algorithm/DDFM/DDFM_T.h"
#include "algorithm/DS/DS.h"
#include "algorithm/RONWM/RONMW.h"
#include "algorithm/WavingSketch/WS_O.h"
#include "algorithm/WavingSketch/WS.h"
#include "algorithm/WavingSketch/WS_OUR.h"

#include "utils/FatTree.h"
#include "utils/SpineLeaf.h"
#include "utils/MurmurHash3.h"

using namespace std;

typedef pair<uint32_t, uint32_t> pii;

// functions
void process_args(int argc, char* argv[]);
void data_preprocess();
void read_packet_data();
void generate_hash_seeds(int len = 8);
bool create_directory(const std::string& dir);
void free();
// core functions
void get_route_FatTree(FatTree& ft);
void get_route_SpineLeaf(SpineLeaf& sl);
void get_heavy_hitter();
template <class T>
void OtherResultAnalysis(T& sketch, string name);

// file
ifstream inf;
ofstream ouf;

// const variable
int data_id = 0, topo_id = 0;
int ft_k = 4, sl_n = 4, sl_m = 8;
string result_dir = "./result/";

// experiment variable
int hh_threshold = -1;
double hh_ratio = 0.0001;
double key_list_size = 512, hh_key_memory = key_list_size * 4 / 1024;
double ds_key_list_size = 512, ds_key_memory = ds_key_list_size * 4 / 1024;
double memory = 72, filter_memory = 5;
double hh_key_th = 8000, re_dis_th = 10000;

// data
int packet_num, flow_num;
uint32_t* seeds;
pii* packets;
uint32_t* flows;
vector<vector<pii>> data_set;
unordered_set<uint32_t> flow_set;
unordered_map<uint32_t, int> actual_flow_size;
vector<pair<int, uint32_t>> actual_freq_vec;
unordered_set<uint32_t> actual_heavy_hitters;
unordered_map<uint32_t, int> routing;
unordered_map<uint32_t, unordered_set<uint32_t>> switch_flow;

int main(int argc, char* argv[]) {
    // process args
    process_args(argc, argv);
    // read data
    read_packet_data();
    // generate hash seeds
    generate_hash_seeds();
    // route
    FatTree ft(ft_k);
    SpineLeaf sl(sl_n, sl_m);
    if (topo_id == 0) {
        printf("------------------ FatTree ------------------\n");
        get_route_FatTree(ft);
        filter_memory = filter_memory * ft.iEdge / (ft.switches_num - 1);
    } else {
        printf("------------------ SpineLeaf ------------------\n");
        get_route_SpineLeaf(sl);
        filter_memory = filter_memory * sl.iLeaf / (sl.switches_num - 1);
    }
    // heavy hitter
    hh_threshold = -1;
    vector<double> hh_ratios = {0.000005, 0.000010, 0.000015, 0.000020, 0.000025, 0.000030};
    vector<vector<double>> memo = {{950, 650, 500, 400, 350, 300}, {650, 300, 200, 150, 150, 150}, {1450, 950, 750, 650, 550, 500}, {1000, 550, 300, 250, 250, 200}};
    vector<vector<double>> ds_memo = {{1600, 1000, 800, 650, 550, 500}, {1600, 700, 400, 250, 200, 200}, {2300, 1500, 1200, 950, 800, 700}, {3800, 1250, 650, 400, 300, 300}};
    vector<double> hh_key_memorys = memo[topo_id * 2 + data_id];
    vector<double> ds_key_memorys = ds_memo[topo_id * 2 + data_id];
    vector<double> re_dis_ths = {15000, 15000, 15000, 15000, 15000, 15000, 15000};
    // var ratio
    for (int index = 0; index < 1; index++) {
        hh_ratio = hh_ratios[index];
        hh_key_th = hh_ratio * packet_num * 0.9;
        key_list_size = hh_key_memorys[index];
        hh_key_memory = key_list_size * 4 / 1024;
        ds_key_list_size = ds_key_memorys[index];
        ds_key_memory = ds_key_list_size * 4 / 1024;
        re_dis_th = re_dis_ths[index];
        get_heavy_hitter();
        printf("hh_ratio: %f, key_list_size: %f, hh_key_memory: %f, hh_key_th: %f, filter_memory: %f\n", hh_ratio, key_list_size, hh_key_memory, hh_key_th, filter_memory);
        string hh_dir = result_dir + to_string(topo_id) + to_string(data_id) + "_HH_" + to_string((hh_threshold == -1 ? hh_ratio : hh_threshold)) + "/";
        create_directory(hh_dir);
        // var memory
        for (memory = 24; memory <= 96; memory += 12) {
            // ====================== Algorithm ======================
            printf("================== Memory: %.6fKB(%.6f)================\n", memory, (hh_threshold == -1 ? hh_ratio : hh_threshold));
            ouf.open(hh_dir + to_string(memory) + ".txt", ios::out);
            if (topo_id == 0) {
                // DDFM
                DDFM_T ddfm(ft_k, memory, key_list_size, seeds, hh_key_th, re_dis_th);
                OtherResultAnalysis(ddfm, "DDFM");
                // CMAX
                CMAX cmax(ft_k, memory + filter_memory + hh_key_memory, seeds);
                OtherResultAnalysis(cmax, "CMAX");
                // RONWM
                RONMW<UWRA> ronwm_uwra(ft_k, memory + filter_memory + hh_key_memory);
                OtherResultAnalysis(ronwm_uwra, "UWRA");
                RONMW<NWPS> ronwm_nwps(ft_k, memory + filter_memory + hh_key_memory);
                OtherResultAnalysis(ronwm_nwps, "NWPS");
                // Distributed Sketch
                DS ds(ft_k, memory + filter_memory + hh_key_memory - ds_key_memory, ds_key_list_size, seeds, hh_ratio * packet_num);
                OtherResultAnalysis(ds, "DS");
                // WavingSketch
                WS_O ws(ft_k, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(ws, "WS");
                WS ws_nspa(ft_k, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(ws_nspa, "WS_NSPA");
                WS_OUR ws_our(ft_k, memory + filter_memory + hh_key_memory, seeds, re_dis_th);
                OtherResultAnalysis(ws_our, "WS_OUR");
                // ChainSketch
                CHAIN_O chain(ft_k, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(chain, "CHAIN");
                CHAIN chain_nspa(ft_k, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(chain_nspa, "CHAIN_NSPA");
                CHAIN_OUR chain_our(ft_k, memory + filter_memory + hh_key_memory, seeds, re_dis_th);
                OtherResultAnalysis(chain_our, "CHAIN_OUR");
            }else{
                // DDFM
                DDFM_T ddfm(sl_n, sl_m, memory, key_list_size, seeds, hh_key_th, re_dis_th);
                OtherResultAnalysis(ddfm, "DDFM");
                // CMAX
                CMAX cmax(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds);
                OtherResultAnalysis(cmax, "CMAX");
                // RONWM
                RONMW<UWRA> ronwm_uwra(sl_n, sl_m, memory + filter_memory + hh_key_memory);
                OtherResultAnalysis(ronwm_uwra, "UWRA");
                RONMW<NWPS> ronwm_nwps(sl_n, sl_m, memory + filter_memory + hh_key_memory);
                OtherResultAnalysis(ronwm_nwps, "NWPS");
                // Distributed Sketch
                DS ds(sl_n, sl_m, memory + filter_memory + hh_key_memory - ds_key_memory, ds_key_list_size, seeds, hh_ratio * packet_num);
                OtherResultAnalysis(ds, "DS");
                // WavingSketch
                WS_O ws(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(ws, "WS");
                WS ws_nspa(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(ws_nspa, "WS_NSPA");
                WS_OUR ws_our(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds, re_dis_th);
                OtherResultAnalysis(ws_our, "WS_OUR");
                // ChainSketch
                CHAIN_O chain(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(chain, "CHAIN");
                CHAIN chain_nspa(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds, data_id);
                OtherResultAnalysis(chain_nspa, "CHAIN_NSPA");
                CHAIN_OUR chain_our(sl_n, sl_m, memory + filter_memory + hh_key_memory, seeds, re_dis_th);
                OtherResultAnalysis(chain_our, "CHAIN_OUR");
            }
            ouf.close();
            // ======================================================
        }
    }
    // free
    free();
    return 0;
}

template <class T>
void OtherResultAnalysis(T& sketch, string name) {
    ouf << name << endl;
    printf("+++++++++++> Start %s\n", name.c_str());
    for (const auto& flow : actual_flow_size) {
        sketch.flow_path_id[flow.first] = routing[flow.first];
    }
    sketch.show_config();
    // Insert
    double time = 0;
    clock_t start = clock();
    for (int i = 0; i < packet_num; i++) {
        sketch.send(packets[i].first, packets[i].second);
    }
    clock_t end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Insert throughputs: %.2f pps\n", packet_num / time);
    double max_flow = 0, max_packet = 0, flow_avg = 0, packet_avg = 0;
    double flow_var = 0, packet_var = 0;
    for (int i = 1; i < sketch.switch_num; i++) {
        double fl = sketch.flow_load[i].size();
        double pl = sketch.packet_load[i];
        max_flow = max(max_flow, fl);
        max_packet = max(max_packet, pl);
        flow_avg += fl;
        packet_avg += pl;
        ouf << fl << " " << pl << " ";
    }
    ouf << endl;
    flow_avg /= sketch.switch_num - 1;
    packet_avg /= sketch.switch_num - 1;
    for (int i = 1; i < sketch.switch_num; i++) {
        double fl = (int)sketch.flow_load[i].size();
        double pl = sketch.packet_load[i];
        flow_var += (fl - flow_avg) * (fl - flow_avg);
        packet_var += (pl - packet_avg) * (pl - packet_avg);
    }
    flow_var /= sketch.switch_num - 1;
    packet_var /= sketch.switch_num - 1;
    ouf << max_flow << " " << max_packet << endl;
    ouf << flow_var << " " << packet_var << endl;
    printf("Flow: %.2f, %.2f, %.2f\n", flow_avg, max_flow, sqrt(flow_var));
    printf("Packet: %.2f, %.2f, %.2f\n", packet_avg, max_packet, sqrt(packet_var));
    // heavy hitter
    unordered_map<uint32_t, uint32_t> heavy_hitters = sketch.query_heavy_hitters(hh_threshold == -1 ? int(hh_ratio * packet_num) : hh_threshold);
    double FP = 0, FN = 0, TP = 0, TN = 0, ARE_HH = 0, AAE_HH = 0;
    for (const auto& f : actual_heavy_hitters) {
        uint32_t est = 0;
        if (heavy_hitters.find(f) == heavy_hitters.end()) {
            FN++;
            est = 0;
        } else {
            TP++;
            est = heavy_hitters[f];
        }
        ARE_HH += abs((int)actual_flow_size[f] - (int)est) / (double)actual_flow_size[f];
        AAE_HH += abs((int)actual_flow_size[f] - (int)est);
    }
    ARE_HH /= actual_heavy_hitters.size();
    AAE_HH /= actual_heavy_hitters.size();
    for (const auto& f : heavy_hitters) {
        if (actual_heavy_hitters.find(f.first) == actual_heavy_hitters.end()) {
            FP++;
        }
    }
    TN = actual_flow_size.size() - TP;
    double PR = TP / (TP + FP);
    double RR = TP / (TP + FN);
    double F1 = 2 * PR * RR / (PR + RR);
    double F2 = 5 * PR * RR / (4 * PR + RR);
    double FNR = FN / (FN + TP);
    double FPR = FP / (FP + TN);
    printf("FNR(%f), FPR(%f), PR: %f, RR: %f, F1: %f, F2: %f, ARE_HH: %f, AAE_HH: %f\n", FNR, FPR, PR, RR, F1, F2, ARE_HH, AAE_HH);
    ouf << FNR << " " << FPR << " " << PR << " " << RR << " " << F1 << " " << F2 << " " << ARE_HH << " " << AAE_HH << endl;
    // per flow size
    double ARE_100 = 0, ARE_1000 = 0, ARE_10000 = 0, ARE_ = 0, ARE_all = 0;
    double AAE_100 = 0, AAE_1000 = 0, AAE_10000 = 0, AAE_ = 0, AAE_all = 0;
    int N_100 = 0, N_1000 = 0, N_10000 = 0, N_ = 0, N_all = 0;
    for (const auto& flow : actual_flow_size) {
        uint32_t est = 0;
        if (heavy_hitters.find(flow.first) != heavy_hitters.end()) {
            est = heavy_hitters[flow.first];
        } else {
            est = sketch.query(flow.first);
        }
        ARE_all += abs((int)flow.second - (int)est) / (double)flow.second;
        N_all++;
        AAE_all += abs((int)flow.second - (int)est);
        if (flow.second <= 100) {
            ARE_100 += abs((int)flow.second - (int)est) / (double)flow.second;
            AAE_100 += abs((int)flow.second - (int)est);
            N_100++;
        } else if (flow.second <= 1000) {
            ARE_1000 += abs((int)flow.second - (int)est) / (double)flow.second;
            AAE_1000 += abs((int)flow.second - (int)est);
            N_1000++;
        } else if (flow.second <= 10000) {
            ARE_10000 += abs((int)flow.second - (int)est) / (double)flow.second;
            AAE_10000 += abs((int)flow.second - (int)est);
            N_10000++;
        } else {
            ARE_ += abs((int)flow.second - (int)est) / (double)flow.second;
            AAE_ += abs((int)flow.second - (int)est);
            N_++;
        }
    }
    ARE_100 /= N_100, ARE_1000 /= N_1000, ARE_10000 /= N_10000, ARE_ /= N_;
    AAE_100 /= N_100, AAE_1000 /= N_1000, AAE_10000 /= N_10000, AAE_ /= N_;
    ARE_all /= N_all;
    AAE_all /= N_all;
    printf("ARE_100(%d): %.4f, ARE_1000(%d): %.4f, ARE_10000(%d): %.4f, ARE_(%d): %.4f\n", N_100, ARE_100, N_1000, ARE_1000, N_10000, ARE_10000, N_, ARE_);
    printf("ARE_all(%d): %.4f\n", N_all, ARE_all);
    printf("AAE_100(%d): %.4f, AAE_1000(%d): %.4f, AAE_10000(%d): %.4f, AAE_(%d): %.4f\n", N_100, AAE_100, N_1000, AAE_1000, N_10000, AAE_10000, N_, AAE_);
    printf("AAE_all(%d): %.4f\n", N_all, AAE_all);
    ouf << ARE_all << " " << AAE_all << endl;
}

void get_heavy_hitter() {
    actual_heavy_hitters.clear();
    for (int i = 0; i < actual_freq_vec.size(); i++) {
        if (hh_threshold != -1) {
            if (actual_freq_vec[i].first >= hh_threshold) {
                actual_heavy_hitters.insert(actual_freq_vec[i].second);
            }
        } else {
            if (actual_freq_vec[i].first >= hh_ratio * packet_num) {
                actual_heavy_hitters.insert(actual_freq_vec[i].second);
            } else {
                break;
            }
        }
    }
    if (hh_threshold != -1) {
        printf("Threshold: %d, Heavy Hitters: %d\n", hh_threshold, actual_heavy_hitters.size());
    } else {
        printf("Heavy Hitters(%f, %f): %d\n", hh_ratio, hh_ratio * packet_num, actual_heavy_hitters.size());
    }
}

void get_route_FatTree(FatTree& ft) {
    unordered_map<int, vector<int>> host2path;
    int host = ft.iEdge * 2;
    for (int i = 0; i < host; i++) {
        for (int j = 0; j < host; j++) {
            vector<int> candi;
            for (int k = 0; k < ft.path.size(); k++) {
                int tp1 = ft.path[k][0] - (ft.iCore + ft.iAgg + 1);
                int tp2 = ft.path[k][ft.path[k].size() - 1] - (ft.iCore + ft.iAgg + 1);
                if (tp1 == i / 2 && tp2 == j / 2) {
                    candi.push_back(k);
                }
            }
            host2path[i * 138 + j] = candi;
        }
    }
    unordered_map<int, int> path_cnt;
    for (const auto& flow : actual_flow_size) {
        int host1 = rand() % host, host2 = rand() % host;
        while (host1 == host2)
            host2 = rand() % host;
        vector<int> routes = host2path[host1 * 138 + host2];
        int index = rand() % routes.size();
        routing[flow.first] = routes[index];
        switch_flow[ft.path[routes[index]][0]].insert(flow.first);
        for (int i = 0; i < ft.path[routing[flow.first]].size(); i++) {
            path_cnt[ft.path[routing[flow.first]][i]]++;
        }
    }
}

void get_route_SpineLeaf(SpineLeaf& ft) {
    unordered_map<int, vector<int>> host2path;
    int host = ft.iLeaf * 2;
    for (int i = 0; i < host; i++) {
        for (int j = 0; j < host; j++) {
            vector<int> candi;
            for (int k = 0; k < ft.path.size(); k++) {
                int tp1 = ft.path[k][0] - (ft.iSpine + 1);
                int tp2 = ft.path[k][ft.path[k].size() - 1] - (ft.iSpine + 1);
                if (tp1 == i / 2 && tp2 == j / 2) {
                    candi.push_back(k);
                }
            }
            host2path[i * 138 + j] = candi;
        }
    }
    unordered_map<int, int> path_cnt;
    for (const auto& flow : actual_flow_size) {
        int host1 = rand() % host, host2 = rand() % host;
        while (host1 == host2)
            host2 = rand() % host;
        vector<int> routes = host2path[host1 * 138 + host2];
        int index = rand() % routes.size();
        routing[flow.first] = routes[index];
        switch_flow[ft.path[routes[index]][0]].insert(flow.first);
        for (int i = 0; i < ft.path[routing[flow.first]].size(); i++) {
            path_cnt[ft.path[routing[flow.first]][i]]++;
        }
    }
}

void read_packet_data() {
    // string path = "../data2019/";
    // // string path = "../data2019_temp/";
    // if (data_id == 1) {
    //     path = "../data_DC/";
    //     // path = "../data_DC_temp/";
    // }
    string path = "./data/";
    data_set.clear();

    pii temp;
    int total = 0;

    vector<string> fs = {"00.txt", "01.txt", "02.txt", "03.txt", "04.txt"};
    // vector<string> fs = {"00.txt", "01.txt", "02.txt", "03.txt", "04.txt", "05.txt", "06.txt", "07.txt", "08.txt", "09.txt"};
    if (data_id == 1) {
        fs = {"dc.txt"};
    }
    random_device rd;
    mt19937 rng(rd());
    for (string f : fs) {
        vector<pii> data;
        string c_path = path + f;
        inf.open(c_path, ios::in);
        cout << "Reading " << c_path << endl;
        int cc = 0;
        while (inf >> temp.first && inf >> temp.second) {
            cc += 1;
            total += 1;
            // temp.second = total;
            temp.second = rng();
            data.push_back({temp.first, temp.second});
            flow_set.insert(temp.first);
            actual_flow_size[temp.first] += 1;
        }
        cout << "Read " << cc << " lines" << endl;
        data_set.push_back(data);
        inf.close();
        // break;
    }
    packet_num = total;
    flow_num = flow_set.size();

    cout << "Total data size: " << total << ", flow size: " << flow_num << endl;

    data_preprocess();
}

void data_preprocess() {
    packets = new pii[packet_num];
    memset(packets, 0, packet_num * sizeof(pii));
    int index = 0;
    for (const auto& data : data_set) {
        for (const auto& d : data) {
            packets[index++] = d;
        }
    }
    flows = new uint32_t[flow_num];
    memset(flows, 0, flow_num * sizeof(uint32_t));
    index = 0;
    for (const auto& flow : flow_set) {
        flows[index++] = flow;
    }
    for (const auto& ele : actual_flow_size) {
        actual_freq_vec.push_back({ele.second, ele.first});
    }
    sort(actual_freq_vec.begin(), actual_freq_vec.end(), greater<pair<int, uint32_t>>());
}

void process_args(int argc, char* argv[]) {
    // if (argc < 3) {
    if (argc < 2) {
        printf("Usage: %s <topo_id> \n", argv[0]);
        exit(1);
    }
    topo_id = atoi(argv[1]);
    // data_id = atoi(argv[2]);
    data_id = 1;
}

void generate_hash_seeds(int len) {
    seeds = new uint32_t[len];
    int count = 0;
    std::unordered_set<int> diff_ele;
    while (count < len) {
        int num = rand();
        if (diff_ele.find(num) == diff_ele.end()) {
            diff_ele.insert(num);
            seeds[count++] = num;
        }
    }
}

bool create_directory(const std::string& dir) {
#ifdef _WIN32
    if (_mkdir(dir.c_str()) == 0 || errno == EEXIST) {
#else
    if (mkdir(dir.c_str(), 0777) == 0 || errno == EEXIST) {
#endif
        return true; 
    } else {
        return false; 
    }
}

void free() {
    delete[] packets;
    delete[] flows;
    delete[] seeds;
}