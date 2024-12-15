# Lightweight Two-level Collaborative Network Traffic Measurement for Data Center Networks

## Introduction

Network traffic measurement is crucial for the effective management of data center networks. Collaborative measurement solutions distribute measurement tasks to switches based on flow-level or packet-level granularity to alleviate the measurement load on each switch. However, flow-level solutions often experience severe imbalances in measurement overhead between switches measuring large or small flows, and face scalability challenges due to the costly optimization of collaborative plans. Additionally, packet-level solutions do not adequately reduce hash collisions in sketches and fail to significantly enhance measurement accuracy. In this paper, we present the Lightweight Two-level Collaborative Measurement (LTCM) that synergies flow-level and packet-level strategy to optimize measurement load balancing, reduce overall measurement overhead, and enhance measurement accuracy. We first design a Lightweight Flow-level Measurement (FCM) framework that balances the number of flows measured by each switch, incorporating a novel interval-matching technique that significantly lowers the computational costs of collaborative strategies. Based on FCM, LTCM implements our measurement load balancing strategy selector at ingress switches to detect flows whose number of packets entering the network exceeds a given threshold and evenly distribute their subsequent packets across all switches for measurement. To further improve measurement accuracy and speed, we design a two-layer collaborative sketch that reduces hash collisions between large and small flows. We implement LTCM on a Tofino-based programmable switch. Experimental results based on real Internet traces show that LTCM achieves highly efficient flow-level and packet-level measurement load balancing, improving accuracy by up to $96.6$\% and throughput by $112.79$\%.

## About this repo

The core LTCM algorithm is implemented in /algorithm/DDFM.

Other baseline methods are also implemented in /algorithm.

The dataset files are placed under the /data.

The general modules are placed under the /utils.

The main function is `main.cpp`, which is used to measure the accuracy of various algorithms. The file `throughput.cpp` is used to measure the throughput of these algorithms.

## Dataset

Please note that the dataset needs to have its IP addresses converted to integers beforehand. Each line should consist of two unsigned integers separated by a space: the first number representing the source IP address and the second number representing the destination IP address.

Additionally, due to the large size of the original dataset, we do not provide the dataset.

## How to build

You can use the following commands to build and run.

```bash
g++ ./main.cpp -o main -std=c++17 -I -w --O3
./main 0
```

> Note: The parameter `0` in the last command is used to specify the topology of the data center network. `0` represents the FatTree topology, and `1` represents the SpineLeaf topology.