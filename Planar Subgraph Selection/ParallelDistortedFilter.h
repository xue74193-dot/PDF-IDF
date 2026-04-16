#ifndef PLA_DISTORTEDFILTERING_H
#define PLA_DISTORTEDFILTERING_H
#include "utilityfunction.h"

using namespace std;
class ALG_known
{
public:
    double tau;
    vector<S_class> E;
    ALG_known(const double t)
    {
        tau=t;
    }
};
class ALG_Distorted
{
public:
    ALG_Distorted(){};
    vector<ALG_known> alg;
    double max_tau=0.0;//max tau in C, help to update tau
    int min_tau_index_in_C=-1;//-1 represent C is empty, else represent min tau index which is still available now
    int max_tau_index_in_C=-1;
};

default_random_engine e(12345);
Result ParallelDistortedFilter(double eps,int  max_edge) {

    bernoulli_distribution ran(1.0);

    cout << "ParallelDistortedFilter & Max edge: " << max_edge << "---------start---------" << endl;

    long long int oracle_times = 0;

    double delta = 3.0;
    double alpha_1 = 1.0 + 1.0 / (k + 1.0);


    cout << "alpha: " << alpha_1 << endl;
    S_class S_best;
    ALG_Distorted Gamma;
    S_class G;
    /*****only used for generate C, its first value seems useless?*****/
    pair<Edge, double> M = make_pair(Edge(), 0.0);
    int counter = 0;


    for (const auto &u: Groundset) {

        {

            /*******k-system constraint***********/
            if (G.is_feasible(u, max_edge)) {
                //marginal is useless for G, so we use a arbitrary value
                G.add_element(0.0, u);
            }
            /*****singleton must satisfy planar, so only check budget******/
            double fu_value = f_u(u);
            if (fu_value > M.second) {
                M.first = u;
                M.second = fu_value;
            }
            oracle_times++;
            if (fu_value > S_best.S_revenue) {
                S_best.replace_with_singleton(fu_value, u);
            }
            if (fu_value <= 0) continue;

            //generate ALG2
            double left_temp = M.second;
            double right_temp = M.second * pow(1.0 + delta, 1.0 + 2.0 * log(k) / log(1.0 + delta) +
                                                            2.0 * log(G.edge_set.size()) / log(1.0 + delta));

            if (left_temp > right_temp) continue;

            double left = ceil(log(left_temp) / log(1.0 + delta));
            double right = floor(log(right_temp) / log(1.0 + delta));

            if (Gamma.min_tau_index_in_C == -1)//first visit C
            {
                for (int t = left; t <= right; t++) {
                    double tau_temp = pow(1.0 + delta, t);
                    Gamma.alg.push_back(ALG_known(tau_temp));
                }
                Gamma.min_tau_index_in_C = 0;//the index of min tau which is still available now
                Gamma.max_tau = pow(1.0 + delta, right);//now max tau
                Gamma.max_tau_index_in_C = Gamma.alg.size() - 1;//max tau index in C now
            } else//not first visit
            {
                double now_min_tau_in_C = pow(1.0 + delta, left);//min tau in C now
                double now_max_tau_in_C = pow(1.0 + delta, right);//max tau in C now
                if (now_min_tau_in_C > Gamma.max_tau)//all old S_array should be removed
                {
                    Gamma.min_tau_index_in_C = Gamma.alg.size();//the index of min tau which is still available now
                    for (int t = left; t <= right; t++)//then add new S gamma pair
                    {
                        double tau_temp = pow(1.0 + delta, t);
                        Gamma.alg.push_back(ALG_known(tau_temp));
                    }
                    Gamma.max_tau_index_in_C = Gamma.alg.size() - 1;//max gamma index in C now
                    Gamma.max_tau = Gamma.alg.back().tau;//the last element always is the max_gamma anyway
                } else//else find where is the min gamma index now, which is equivalent to remove all S whose gamma < left
                {
                    bool need_update = true;//judge S_array should be updated or not
                    for (int z = Gamma.min_tau_index_in_C; z < Gamma.alg.size(); z++) {
                        if (Gamma.alg[z].tau < now_min_tau_in_C)
                            Gamma.min_tau_index_in_C++;
                        if (Gamma.alg[z].tau >= now_max_tau_in_C) {
                            Gamma.max_tau_index_in_C = z;
                            need_update = false;
                            break;
                        }
                    }
                    if (need_update) {
                        //finally, go through all gamma now, put new pair (S,gamma) if needed
                        for (int t = left; t <= right; t++) {
                            double tau_temp = pow(1.0 + delta, t);
                            if (tau_temp > Gamma.max_tau)
                                Gamma.alg.push_back(ALG_known(tau_temp));
                        }
                        Gamma.max_tau_index_in_C = Gamma.alg.size() - 1;//max gamma index in C now
                        Gamma.max_tau = Gamma.alg.back().tau;//the last element always is the max_gamma anyway
                    }
                }
            }

            for (int b = Gamma.min_tau_index_in_C; b <= Gamma.max_tau_index_in_C; b++) {
                //if F_max<=0, then continue;
                if (Gamma.alg[b].tau <= 0) continue;

                int m = ceil(1.0 + 2.0 * log(k) / log(1.0 + delta) + 2.0 * log(G.edge_set.size()) / log(1.0 + delta));
                if (m > Gamma.alg[b].E.size()) {
                    for (int i = Gamma.alg[b].E.size(); i < m; i++) {
                        Gamma.alg[b].E.push_back(S_class());
                    }
                }

                for (int i = 1; i <= m; i++) {
                    //double alpha_i=i*alpha_1;
                    double alpha_i = alpha_1;  //β
                    double threshold_i = Gamma.alg[b].tau / pow(1.0 + delta, i - 1.0);

                    /*****this part should be different for different applications****************/
//                    S_class Cup;
//                    for (int x = 0; x < i; x++)//all sets are disjoint
//                        Cup.edge_set.insert(Cup.edge_set.end(), Gamma.alg[b].E[x].edge_set.begin(),
//                                            Gamma.alg[b].E[x].edge_set.end());

                    S_class Cup;
                    for(int x=0;x<i;x++){

                        for(const auto& edge:Gamma.alg[b].E[x].edge_set){
                            Cup.edge_set.push_back(edge);

                            Cup.node_set.insert(edge.x);
                            Cup.node_set.insert(edge.y);
                        }

                    }

                    double distorted = Cup.marginal_gain_f_distorted(u, alpha_i);
                    oracle_times++;

                    bool flag = false;
                    if (distorted >= threshold_i) {
                        flag = true;
                        if (Gamma.alg[b].E[i - 1].is_feasible(u, max_edge)) {
                            //the marginal gain input g()-c();
                            if (ran(e) == 1) {

                              //  cout<<"threshold_i: "<<threshold_i<<endl;
                                /*******can be further improved********/
                              //  double Before_adding_to_set = Gamma.alg[b].E[i - 1].f_S();
                                double single_set_marginal = Gamma.alg[b].E[i - 1].marginal_gain_f(u);
                                Gamma.alg[b].E[i - 1].add_element(single_set_marginal, u);
                               // cout<<"single_set_marginal:"<<single_set_marginal<<endl;


//                                if(single_set_marginal==Gamma.alg[b].E[i - 1].f_S()-Before_adding_to_set)  cout<<"1"<<" ";
//                                cout<<"single_set_marginal: "<< single_set_marginal<<"after_adding_to_set-Before_adding_to_set: "<< Gamma.alg[b].E[i - 1].f_S()-Before_adding_to_set<<" "<<endl;

                                //update S*
                                if (Gamma.alg[b].E[i - 1].S_revenue > S_best.S_revenue) {
                                    S_best = Gamma.alg[b].E[i - 1];
                                }
                            }
                            break;
                        }
                    }

//                if(flag)
//                    break;
                }

            }
        }
    }
//        cout << "max_tau_index: " << Gamma.max_tau_index_in_C << endl;
//        cout << "min_tau_index: " << Gamma.min_tau_index_in_C << endl;
        int m = ceil(1.0 + 2.0 * log(k) / log(1.0 + delta) + 2.0 * log(G.edge_set.size()) / log(1.0 + delta));
//        cout << "m: " << m << endl;
/*        for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
       {
           for(int i=0;i<m;i++)
           {
               cout<<"  revenue: "<<Gamma.alg[b].E[i].S_revenue<<"  Real revenue: "<<Gamma.alg[b].E[i].f_S(data)<<" cost: "<<Gamma.alg[b].E[i].S_cost<<" size: "<<Gamma.alg[b].E[i].Set.size()<<endl;
               for(const auto &p:Gamma.alg[b].E[i].Set)
                   cout<<p<<" ";
               cout<<endl;
               cout<<" genres_sum: "<<endl;
               for(int j=0;j<Gamma.alg[b].E[i].genres_sum.size();j++)
               {
                   cout<<Gamma.alg[b].E[i].genres_sum[j]<<" ";
               }
               cout<<endl;
           }
       }*/
        long long int saved_element = 0;
        for (int b = Gamma.min_tau_index_in_C; b <= Gamma.max_tau_index_in_C; b++) {
            for (const auto &p: Gamma.alg[b].E) {
                S_class temp_best = p;

                for (const auto &it1: Gamma.alg[b].E) {
                    for (const auto &it2: it1.edge_set) {
                        if (temp_best.is_feasible(it2, max_edge)) {
                            double temp_marginal = temp_best.marginal_gain_f(it2);
                            if (temp_marginal > 0)


                                temp_best.add_element(0.0, it2);
                        }
                    }
                }

                //sum element/memory
                saved_element += temp_best.edge_set.size();

                temp_best.S_revenue = temp_best.f_S();

                if (temp_best.S_revenue > S_best.S_revenue)
                    S_best = temp_best;
            }


        }


        /*****call offline for every gamma *******/
/*    for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
    {
        vector<int> U;
        for(const auto &p:Gamma.alg[b].E) {
            U.insert(U.end(), p.Set.begin(), p.Set.end());
        }
        pair<S_gamma,long long int> temp=Boost(eps,data,U,alpha_1,M);
        oracle_times+=temp.second;
        if (temp.first.S_revenue > S_best.S_revenue)
            S_best = temp.first;

*//*        cout<<"offline result: "<<endl;
        cout<<"  revenue: "<<temp.first.S_revenue<<"  Real revenue: "<<temp.first.f_S(data)<<" cost: "<<temp.first.S_cost<<" size: "<<temp.first.Set.size()<<endl;
        for(const auto &p:temp.first.Set)
            cout<<p<<" ";
        cout<<endl;
        cout<<" genres_sum: "<<endl;
        for(int j=0;j<temp.first.genres_sum.size();j++)
        {
            cout<<temp.first.genres_sum[j]<<" ";
        }
        cout<<endl;*//*
    }*/
        /*****call offline for all existing sets *******/
/*    vector<int> U;
    for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
    {

        for(const auto &p:Gamma.alg[b].S) {
            U.insert(U.end(), p.solution.begin(), p.solution.end());
        }
    }
    pair<S_class,long long int> temp=Boost(eps,U,max_node,alpha_1);
    oracle_times+=temp.second;
    if (temp.first.S_revenue > S_best.S_revenue)
        S_best = temp.first;*/


//        cout << "S*:" << endl;
//        cout << "  revenue: " << S_best.S_revenue << " cost_1: " << S_best.S_cost[0]<<"cost_2: " << S_best.S_cost[1]<< " size: "
//             << S_best.edge_set.size() << endl;
//        for (int i = 0; i < S_best.edge_set.size(); i++)
//            cout << S_best.edge_set[i].x << " " << S_best.edge_set[i].y << " ";
//        cout << endl;
//        cout << "Graph Edge: " << endl;
//        for (const auto &p: S_best.edge_set)
//            cout << "(" << p.x << "," << p.y << ")" << p.cost[0]<<" ";
//        cout << endl;
//        cout << "Graph Node: " << endl;
//        for (const auto &p: S_best.node_set)
//            cout << p << " ";
//        cout << endl;

        cout << "Objective Values: "<< S_best.S_revenue << endl;
        cout << "Oracle Queries: " << oracle_times << endl;

        cout<<"ParallelDistortedFilter ----------------end---------------"<<endl<<endl;
        return Result(S_best.S_revenue, S_best.S_cost, S_best.edge_set.size(), oracle_times);

}
#endif //PLA_DISTORTEDFILTERING_H






