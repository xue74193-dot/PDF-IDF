
#ifndef COVERAGE_ITERATIVEDISTORTEDFILTER_H
#define COVERAGE_ITERATIVEDISTORTEDFILTER_H
#include "haba_non.h"


Result IterativeDistortedFilter(double eps, const int &max_node) {

    cout<<"IterativeDistortedFilter & Max node: "<<max_node<<"---------start---------"<<endl;

    int m = 1;  // number of candidate solutions

    double p=1.0;
    double beta = (k*(1+eps) + 1)/(k*(1+eps));
    double F_max = -99999999999;

    long long int oracle_times = 0;

    bernoulli_distribution ran(p);
    default_random_engine e(1000);

    vector<S_class> S_sets;
    S_class S_0 = S_class();
    S_class S_best = S_class();

    cout << "m: " << m << endl << "beta: " << beta << endl;


    // First pass to find the best singleton
    for (int u = 0; u < node_num; u++) {

        double fu_value = f_u(u);  // g({u}) - c({u})
        oracle_times++;
        if (fu_value > F_max) {
            F_max = fu_value;
            S_0.replace_with_singleton(fu_value, u);
        }

    }

    S_best = S_0;
    S_sets.resize(m,S_class());

    double r = F_max;

    // double r = (eps / (node_num * (1.0 + eps))) * F_max;

    vector<int> belongs(node_num, -1);  //Recall processed elements.



    vector<int> current_version(m, 0); // Increment version number of $S_i$ upon element addition
    vector<vector<double>> last_gain(m, vector<double>(node_num, 0.0));//Record the latest marginal gain respect to S_i
    vector<vector<int>> cache_version(m, vector<int>(node_num, -1));// Record the version of S_i at which the value was cached (initially -1 to denote uncomputed).


    while (r > (eps / (node_num * (1.0 + eps))) * F_max) {


        for (int u = 0; u < node_num; u++) {

            if (belongs[u] != -1) continue;

            int best_set = -1;
            double best_delta = -1e18;

            for (int i = 0; i < m; i++) {
                if (S_sets[i].is_feasible(u, max_node)) {   //check constraint

                    double delta;
                    bool needs_recalc = true;

                    if (cache_version[i][u] == current_version[i]) {
                        needs_recalc = false;

                        delta = last_gain[i][u];
                    }
                    else if (cache_version[i][u] != -1 && last_gain[i][u] < r) {
                        needs_recalc = false;
                        delta = last_gain[i][u];
                    }
                    if (needs_recalc) {
                        delta = S_sets[i].marginal_gain_distorted(u, beta);
                        oracle_times++;

                        last_gain[i][u] = delta;
                        cache_version[i][u] = current_version[i];
                    }

                    if(delta > best_delta){
                        best_delta = delta;
                        best_set = i;
                    }
                }
            }
            if( best_set != -1 && best_delta >= r ){
                if(ran(e)==1) {

                    double add_marginal = S_sets[best_set].marginal_gain(u);
                    oracle_times++;
                    S_sets[best_set].add_element(add_marginal, u);
                    current_version[best_set]++;

                }
                belongs[u] = best_set;
            }
        }

        r = r / (1.0 + eps);
    }

    //     select best solution
    for (int i = 0; i < m; i++) {
        if (S_sets[i].S_revenue > S_best.S_revenue)
            S_best = S_sets[i];
    }
    if (S_0.S_revenue > S_best.S_revenue)
        S_best = S_0;


//    cout<<"F_max: "<<F_max<<endl;
//    cout<<"min_r: "<<r<<endl;
//    cout << "S*:" <<endl;
//    cout << "  revenue: "<<S_best.S_revenue << "  cost: " << S_best.S_cost << "  size: " << S_best.solution.size() << endl;
//    cout << "  Elements: ";
//    for (const auto &p:S_best.solution)
//        cout<<p<<", ";
//    cout << endl;
//    cout << "  group limit: ";
//    for (int j = 0; j < S_best.group_exist.size(); j++) {
//        cout << S_best.group_exist[j] << " ";
//    }
//    cout << endl;



    cout << "Objective Values: " << S_best.S_revenue << endl;
    cout << "Oracle Queries: " << oracle_times << endl;

    cout<<"IterativeDistortedFilter ---------end--------- "<<endl<<endl;
    return Result(S_best.S_revenue, S_best.S_cost, S_best.solution.size(), oracle_times, 0);
}

#endif //COVERAGE_ITERATIVEDISTORTEDFILTER_H
