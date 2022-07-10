#ifndef CCCPATH_H
#define CCCPAHT_H

#include "../../graph/graph.hpp"
#include "../../tools/type.hpp"
#include "../../tools/hopstotchHash.hpp"
#include "../../tools/linearSet.hpp"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <boost/functional/hash.hpp>

using Pair = std::pair<v_size, v_size>;
using std::vector;
using std::unordered_map;

// constexpr v_size batchSize = 50;

template <typename Container> // we can make this generic for any container [1]
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

struct cccpath {
    v_size sz; //sz is the size of S
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double * experiments;
    u_int64_t * exp;
    double sumW;
    double ** dp;
    double * memoryPool = nullptr;
    u_int64_t suW;
    u_int64_t * c;

    v_size * pEdge = nullptr;
    v_size * pIdx = nullptr;
    v_size * pColor = nullptr;
    v_size * sortByColor = nullptr;
    v_size vCnt, eCnt;

    v_size * clique = nullptr;
    std::default_random_engine e;
    e_size N = 5000000;
    unordered_map<vector<v_size>, u_int64_t, container_hash<vector<v_size>>> dpm;

    void init(v_size sz_, std::vector<v_size> & nodes, e_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];
        exp = new u_int64_t[sz];

        // auto cmp = [&](v_size a, v_size b) {
        //     return g->color[a] > g->color[b];
        // };
        
        for(v_size i = 0; i < sz; i++) { // for each node in s
            v_size u = nodes[i]; // guess nodes represent S
            printf("u: %u\n", u);
            // sortByColor = g->pEdge + g->pIdx2[u];

            // std::sort(sortByColor,
            //     sortByColor + g->pIdx[u + 1] - g->pIdx2[u], cmp);
            // sortGraph(u);
            
            memset(pColor, 0, sizeof(v_size)*(g->cc+1)); // pColor = {0,0,0,0,0,0,0,...,0}
            v_size deg = g->pIdx[u+1] - g->pIdx2[u];
            for(v_size i = 0; i < deg; i++) {
                pColor[g->color[g->pEdge[g->pIdx2[u] + i]] + 1]++; // for each out-neighbour of u
                // set its corresponding color place in pColor 1 
                // e.g. neighbours with colors 0,1,3, then pColor is {0,1,1,0,1,0...}
            }

            for(v_size i = 1; i < g->cc; i++) {
                pColor[i] += pColor[i - 1]; // get the cumulative {0,1,2,2,3,3...}
            }

            for(v_size i = 0; i < deg; i++) {
                v_size v = g->pEdge[g->pIdx2[u] + i]; // v is each out-neighbour of u
                sortByColor[ pColor[g->color[v]]++ ] = v; //sortByColor is sorted list of out-neighours of u by color
                // {1212,9302,142,3802,...} by color
            }
            
            memcpy(g->pEdge + g->pIdx2[u], sortByColor, sizeof(v_size)*deg);
            // important step!!! now g->pEdge outneighbours of u are sorted by color!!!! In order #####
            // g->pEdge has same content as sortByColor
            computeDP(u);

            //double sumD = 0.0;
            uint64_t suD = 0;
            /*for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
                sumD += dp[i][k];
            }*/

            for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
                for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    suD += dpm[{i, pEdge[l], k}];
                }
            }

            //experiments[i] = sumD; // experiments[i] stores the total number of k-paths in the graph sum(dp[i][k])
            exp[i] = suD;
            suW += suD;
            //sumW += sumD; // sumW is the total number of k-paths in S
        }
        printf("finished those loops\n");
        delete [] sortByColor;
    }

    // careful, k here is k - 1
    void initForSingleNode(v_size k_, Graph * g_, hopstotchHash * hashTable_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        sumW = 0.0;
        clique = new v_size[k];
/*
        dp = new double*[g->degeneracy];
        memoryPool = new double[g->degeneracy * (k+1)]();
        v_size p = 0;
        for(v_size i = 0; i < g->degeneracy; i++) {
            dp[i] = memoryPool + p;
            p += k + 1;
        }

        for(v_size i = 0; i < g->degeneracy; i++) {
            dp[i][0] = 0;
            dp[i][1] = 1;
        }
*/
        pEdge = new v_size[g->degeneracy*g->degeneracy];
        pIdx = new v_size[g->degeneracy + 1];
        sortByColor = new v_size[g->degeneracy + 1];
        pColor = new v_size[g->cc + 1];
    }

    // void sortGraph(v_size u) {
        // v_size deg = g->pIdx[u+1] - g->pIdx2[u];
        // memcpy(sortByColor, 
        //     g->pEdge + g->pIdx2[u], sizeof(v_size)*deg);
        // auto cmp = [&](v_size a, v_size b) {
        //     return g->color[a] > g->color[b];
        // };
        // std::sort(sortByColor, sortByColor + deg, cmp);

//         pIdx[0] = 0;
//         for(v_size i = 0; i < deg; i++) {
//             v_size v = sortByColor[i];
//             pIdx[i + 1] = pIdx[i]; 

//             for(v_size j = i + 1; j < deg; j++) {
//                 v_size w = sortByColor[j];
// // assert(g->color[v] > g->color[w]);
//                 if(g->color[v] == g->color[w]) continue;
//                 if(hashTable[v].contain(w)) {
// // assert(g->color[v] > g->color[w]);
//                     pEdge[pIdx[i + 1]++] = j;
//                 }
//             }
//         }
    // }

    ~cccpath() {
        if(experiments != nullptr) delete [] experiments;
        if (exp != nullptr) delete [] exp;
        if(memoryPool != nullptr) delete [] memoryPool;
        if(dp != nullptr) delete [] dp;
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
        // if(sortByColor != nullptr) delete [] sortByColor;
        if(clique != nullptr) delete [] clique;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }
    
    void computeDP(v_size u) {
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
        pIdx[0] = 0;
        for(v_size i = 0; i < outDegree; i++) {
            v_size v = sortByColor[i]; // this is the v_i
            pIdx[i + 1] = pIdx[i]; // step up from previous entry

            for(v_size j = i + 1; j < outDegree; j++) {
                v_size w = sortByColor[j];
// assert(g->color[v] > g->color[w]);
                if(g->color[v] == g->color[w]) continue;
                if(hashTable[v].contain(w)) {
                // if there is an edge from v to w
                
// assert(g->color[v] > g->color[w]);

                    //////////////////////////////Important/////////////////////////
                    // ++ because number of increments represent number of outgoingedges of v, hence the gap in pIdx
                
                    pEdge[pIdx[i + 1]++] = j; // Here constructs the local pEdge!!!!, so the value is a corresponding node in local index
                }
            }
        }

        // After above steps
        // pEdge is populated with nodes in local index
        // pIdx is linked to pEdge
        // pIdx is in the same order as the nodes in SortByColor

        // below are real DP process described in the algorithm

        for (v_size i = 0; i < outDegree; i++) {
            for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    dpm[{i, pEdge[l], 1}] = 1;
                }
        }

        for (v_size j = 2; j <= k; j++) {
            for (v_size i = 0; i < outDegree; i++) {
                for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    v_size x = pEdge[l];
                    for (v_size p = pIdx[x]; p < pIdx[x + 1]; p++) {
                        v_size t = pEdge[p];
                        if (dpm[{i, t, 1}] == 1) {
                            dpm[{i, x, j}] = dpm[{i, x, j}] + dpm[{x, t, j-1}];
                        }
                    }
                }
            }
        }



/*    
        for(v_size j = 2; j <= k; j++) {
            for(v_size i = 0; i < outDegree; i++) {
                dp[i][j] = 0.0;
                for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    // for out-neighbour of v_i 
                    // a out-neighbour is pEdge[l]

                    dp[i][j] += dp[pEdge[l]][j - 1];
                }
            }
        }
*/
    }

    // sortByColor a mapping of rank v_i to actual node id (the value of entry)
    int sampleOneTime(v_size id, v_size u, 
        std::uniform_real_distribution<double> & d) { // uiDistribution(0, 1);
        v_size preId = -1;
        v_size prId = -1;

        double sumD = experiments[id];
        double x = d(e); // get a random double

        uint64_t suD = exp[id];
        
        
        double sumTmp = 0.0;
        uint64_t sumT = 0;
        v_size deg = g->pIdx[u+1] - g->pIdx2[u]; // out degree of node u

        v_size last;
        v_size secLast;

        for (v_size i = 0; i < k; i++) {
            if (i == 0) {
                for (v_size j = 0; j < deg; j++) {
                    sumT += c[j];
                    if (sumT + 1e-10 >= x * suD) {
                        clique[0] = sortByColor[ j ];
                        prId = j;
                        break;
                    }
                }
            } else if (i == 1) {
                sumT = 0;
                suD = c[last];
                for (v_size j = pIdx[last]; j < pIdx[last + 1]; j++) {
                    sumT += dpm[{last, j, k-1}];
                    if (sumT + 1e-10 >= x * suD) {
                        clique[1] = sortByColor[ pEdge[j] ];
                        prId = pEdge[j];
                        break;
                    }
                }
            } else {
                sumT = 0;
                suD = dpm[{secLast, last, k - i + 1}];
                for (v_size j = pIdx[last]; j < pIdx[last + 1]; j++) {
                    sumT += dpm[{last, j, k-i}];
                    if (sumT + 1e-10 >= x * suD) {
                        clique[i] = sortByColor[ pEdge[j] ];
                        prId = pEdge[j];
                        break;
                    }
                }

            }

            for(v_size j = 0; j < i-1; j++) {
                if(!connect(clique[i], clique[j])) {
                    return 0;
                }
            }

            secLast = last;
            last = prId;    
        }



        // This is to select the first node, "the start node" of the clique
        /*for (v_size i = 0; i < deg; i++) {
            sumTmp += dp[i][k];
            if(sumTmp + 1e-10 >= x * sumD) {
                clique[0] = sortByColor[ i ];
                preId = i; //preId is the start node Id, a node in S
                break;
            }
        }

        for(v_size i = 1; i < k; i++) {
            sumTmp = sumD = 0.0;
            for(v_size j = pIdx[preId]; j < pIdx[preId + 1]; j++) {
                sumD += dp[pEdge[j]][k - i]; // sumD is cnt in the pseudocode
            }
// bool f = false;
            x = d(e); // sample again
            for(v_size j = pIdx[preId]; j < pIdx[preId + 1]; j++) { // pEdge[j] is each neighbour node, this iterates over neighbour of preId node
                sumTmp += dp[pEdge[j]][k - i];
                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[i] = sortByColor[ pEdge[j] ]; // sortByColor[ pEdge[j] ] this is the actual id of v_j, retrieves the global node index
                    preId = pEdge[j]; 
// f = true;
                    break;
                }
            }
// assert(f);
            for(v_size j = 0; j < i-1; j++) {
                if(!connect(clique[i], clique[j])) {
                    return 0;
                }
            }
        }*/
        
        return 1;
    }

    double sample(std::vector<v_size> & nodes, e_size sampleTimes, double expectedN) {
        e_size t = 0;
        e_size sampleTotalTimes = 0;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        double ans = 0.0;
        printf("abc\n"); 


        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            //printf("j %u\n", u);
            //e_size expectedSampleTime
                //= std::round(sampleTimes * (experiments[i] / sumW) + 0.000001);
            e_size expectedSampleTime = std::round(sampleTimes * (exp[i] / sumW) + 0.000001);
            // expected SampleTime is the expected sample time for sampling around node u
            printf("exp s t: %u | exp: %u | sumW: %u\n", expectedSampleTime, exp[i], sumW);
            if(expectedSampleTime == 0) continue;

            c = new u_int64_t[g->pIdx[u+1] - g->pIdx2[u]];

            sortByColor = g->pEdge + g->pIdx2[u];
            // sortGraph(u);
            printf("just %u\n", u);
            computeDP(u);

            // first mistake made
            for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
                c[i] = 0;
                for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    c[i] += dpm[{i, pEdge[l], k}];
                }
            }

            printf("si: %u\n", i);
            fflush(stdout);

            v_size tt = 0;
            for(v_size j = 0; j < expectedSampleTime; j++) {
                tt += sampleOneTime(i, u, uiDistribution);
            }
            t += tt;
            //ans += 1.0*tt/expectedSampleTime*experiments[i];
            ans += 1.0*tt/expectedSampleTime*exp[i];
            sampleTotalTimes += expectedSampleTime;
            if(c != nullptr) delete [] c;
        }
        
        if(sampleTotalTimes < sampleTimes) {
            //std::discrete_distribution<int> 
              //udistribution(experiments, experiments + sz);
            std::discrete_distribution<int> udistribution(exp, exp + sz);
printf("|not expected %llu ", sampleTimes - sampleTotalTimes);
              while(sampleTotalTimes < sampleTimes) { // if sample number not met target, pick ones following the distribution
                  int id = udistribution(generator);
                  v_size u = nodes[id];
                  sortByColor = g->pEdge + g->pIdx2[u];
                  // sortGraph(u);
                  c = new u_int64_t[g->pIdx[u+1] - g->pIdx2[u]];
                  computeDP(u);

                  for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
                      c[i] = 0;
                      for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                          c[i] += dpm[{i, pEdge[l], k}];
                      }
                  }

                  t += sampleOneTime(id, u, uiDistribution);
                  sampleTotalTimes++;
                  if(c != nullptr) delete [] c;
              }
             // printf("|small %.6f %u %u", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
             // return 1.0 * t / sampleTotalTimes * sumW;
        }
        

        // printf("sampleTimes %u\n", sampleTimes);
        // printf("sample rate %f\n", 1.0 * t / sampleTimes);
        printf("| %.6f %u %u", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
        // printf("| %.8f", expectedN / sumW);
        if(c != nullptr) delete [] c;
        return 1.0 * t / sampleTotalTimes * sumW;
        //return ans;
    }
};

#endif
