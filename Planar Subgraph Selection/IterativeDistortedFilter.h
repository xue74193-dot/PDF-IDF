
#ifndef PLANARGRAPH_ITERATIVEDISTORTEDFILTER_H
#define PLANARGRAPH_ITERATIVEDISTORTEDFILTER_H
#include "haba_non.h"

Result IterativeDistortedFilter(double eps, int max_edge) {

    double r;
    double p=1.0;
    double beta = (k*(1+eps)+1)/(k*(1+eps));
    double F_max = -99999999999;

    long long int oracle_times = 0;
    int m = 1;

    bernoulli_distribution ran(p);
    default_random_engine e(12345);

    vector<S_class> S_sets;
    S_class S_0 = S_class();
    S_class S_best = S_class();

    cout<<"IterativeDistortedFilter & Max edge: "<<max_edge<<"---------start---------"<<endl;
    cout << "m: " << m << endl << "beta: " << beta << endl;
    cout<<"p: "<<p<<endl;

    //  First pass to find the best singleton

    for(const auto &u:Groundset){

        double fu_value = f_u(u);  // g({u}) - c({u})
        oracle_times++;
        if (fu_value > F_max) {
            F_max = fu_value;
            S_0.replace_with_singleton(fu_value, u);
        }

    }


    S_best = S_0;
    S_sets.resize(m,S_class());



    double min_r = (eps / (edge_num * (1.0 + eps))) * F_max;
    r = F_max;



    vector<int> belongs(edge_num, -1);  //Recall processed elements.


    vector<int> current_version(m, 0); // Increment version number of $S_i$ upon element addition
    vector<vector<double>> last_gain(m, vector<double>(edge_num, 0.0));//Record the latest marginal gain respect to S_i
    vector<vector<int>> cache_version(m, vector<int>(edge_num, -1)); // Record the version of S_i at which the value was cached (initially -1 to denote uncomputed).

    while (r > min_r) {

        for(const auto &u:Groundset){

            if (belongs[u.index] != -1) continue;

            int best_set = -1;
            double best_delta = -1e18;

            for (int i = 0; i < m; i++) {
                if (S_sets[i].is_feasible(u, max_edge)) {

                    double delta;
                    bool needs_recalc = true;

                    if (cache_version[i][u.index] == current_version[i]) {
                        needs_recalc = false;

                        delta = last_gain[i][u.index];
                    }
                    else if (cache_version[i][u.index] != -1 && last_gain[i][u.index] < r) {
                        needs_recalc = false;
                        delta = last_gain[i][u.index];
                    }
                    if (needs_recalc) {
                        delta = S_sets[i].marginal_gain_f_distorted(u, beta);
                        oracle_times++;

                        last_gain[i][u.index] = delta;
                        cache_version[i][u.index] = current_version[i];
                    }

                    if(delta > best_delta){
                        best_delta = delta;
                        best_set = i;
                    }
                }
            }
            if(best_set!=-1 && best_delta >= r){
                if(ran(e)==1) {

                    double add_marginal = S_sets[best_set].marginal_gain_f(u);
                    oracle_times++;
                    S_sets[best_set].add_element(add_marginal, u);
                    current_version[best_set]++;


                }
                belongs[u.index] = best_set;
            }
        }

        r = r / (1.0 + eps);

    }

    //    select best solution
    for (int i = 0; i < m; i++) {
        if (S_sets[i].S_revenue > S_best.S_revenue)
            S_best = S_sets[i];
    }
    if (S_0.S_revenue > S_best.S_revenue)
        S_best = S_0;



//    cout<<"IterativeDistortedFilter & Max Edge: "<<max_edge<<endl;
//    cout << "S*:" <<endl;
//    cout << "  revenue: "<<S_best.S_revenue << "  cost1: " << S_best.S_cost[0] << "  cost2: " << S_best.S_cost[1] << "  size: " << S_best.edge_set.size() << endl;
//    cout<<"Graph Edge: "<<endl;
//    for(const auto &p:S_best.edge_set)
//        cout<<"("<<p.x<<","<<p.y<<")" <<" ";
//    cout<<endl;
//    cout<<"Graph Node: "<<endl;
//    for(const auto &p:S_best.node_set)
//        cout<<p<<" ";
//    cout<<endl;
//    cout << "  Oracle times: " << oracle_times << endl;
//

    cout << "Objective Values: " << S_best.S_revenue << endl;
    cout << "Oracle Queries: " << oracle_times << endl;

    cout<<"IterativeDistortedFilter ---------end--------- "<<endl<<endl;
    return Result(S_best.S_revenue, S_best.S_cost, S_best.edge_set.size(), oracle_times);
}



#endif //PLANARGRAPH_ITERATIVEDISTORTEDFILTER_H
