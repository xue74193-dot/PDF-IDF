
#ifndef COVERAGE_HABA_NON_H
#define COVERAGE_HABA_NON_H
#include "RandomSelection.h"
class ALG_known_tau
{
public:
    double tau;
    vector<S_class> E;
    ALG_known_tau(const double t)
    {
        tau=t;
    }
};
class ALG
{
public:
    ALG(){};
    vector<ALG_known_tau> alg;
    double max_tau=0.0;//max tau in C, help to update tau
    int min_tau_index_in_C=-1;//-1 represent C is empty, else represent min tau index which is still available now
    int max_tau_index_in_C=-1;
};
class Haba_non
{
public:
    Haba_non(){};
    ALG Gamma;
    S_class G;
    pair<int,double> M=make_pair(-1,0.0);
    int counter=0;
    int ell=-1;
};
Result Haba(double eps,const int &max_node)
{
    cout<<"Haba & Max node: "<<max_node<<"---------start---------"<<endl;
    //int max_movie=999999999;//+inf
    long long int oracle_times=0;
    int h=ceil(log2(2.0*k+1.0));
    cout<<"k:"<<k<<endl;
    cout<<"h: "<<h<<endl;

    int r=4;
    vector<Haba_non> non;
    for(int i=0;i<r;i++)
        non.emplace_back(Haba_non());
    /*****only used for generate C, its first value seems useless?*****/

    int counter=0;
    for(int u=0;u<node_num;u++)
    {
        if(counter%1000==0)
           // cout<<"Counter: "<<counter<<endl;
        counter++;

        for(int it=0;it<r;it++)
        {
            if (non[it].G.is_feasible(u,  max_node)) {
                //marginal is useless for G, so we use a arbitrary value
                non[it].G.add_element(0.0, u);
            }
            /*****singleton must satisfy planar, so only check budget******/
            //if (node_budget_feasible(u, B, data)) {
            double fu_value = f_u(u);
            if (fu_value > non[it].M.second) {
                non[it].M.first = u;
                non[it].M.second = fu_value;
            }
            oracle_times++;
            if(fu_value<=0) continue;

//            if (fu_value > S_best.S_revenue) {
//                S_best.replace_with_singleton(fu_value, u);
//            }
            //}
            //else
            //    continue;

            //generate ALG2
            double left_temp = non[it].M.second;
            //double right_temp=M.second*pow(2.0,4.0+2.0*log2(k)+1.0+2.0*log2(G.Set.size()));
            double right_temp = non[it].M.second * pow(2.0, 5.0 + 2.0 * log2(k) + 2.0 * log2(non[it].G.solution.size()));

            if (left_temp > right_temp) continue;

            double left = ceil(log(left_temp) / log(2.0));
            double right = floor(log(right_temp) / log(2.0));

            if (non[it].Gamma.min_tau_index_in_C == -1)//first visit C
            {
                for (int t = left; t <= right; t++) {
                    double tau_temp = pow(2.0, t);
                    non[it].Gamma.alg.push_back(ALG_known_tau(tau_temp));
                }
                non[it].Gamma.min_tau_index_in_C = 0;//the index of min tau which is still available now
                non[it].Gamma.max_tau = pow(2.0, right);//now max tau
                non[it].Gamma.max_tau_index_in_C = non[it].Gamma.alg.size() - 1;//max tau index in C now
            } else//not first visit
            {
                double now_min_tau_in_C = pow(2.0, left);//min tau in C now
                double now_max_tau_in_C = pow(2.0, right);//max tau in C now
                if (now_min_tau_in_C > non[it].Gamma.max_tau)//all old S_array should be removed
                {
                    non[it].Gamma.min_tau_index_in_C = non[it].Gamma.alg.size();//the index of min tau which is still available now
                    for (int t = left; t <= right; t++)//then add new S gamma pair
                    {
                        double tau_temp = pow(2.0, t);
                        non[it].Gamma.alg.push_back(ALG_known_tau(tau_temp));
                    }
                    non[it].Gamma.max_tau_index_in_C = non[it].Gamma.alg.size() - 1;//max gamma index in C now
                    non[it].Gamma.max_tau = non[it].Gamma.alg.back().tau;//the last element always is the max_gamma anyway
                } else//else find where is the min gamma index now, which is equivalent to remove all S whose gamma < left
                {
                    bool need_update = true;//judge S_array should be updated or not
                    for (int z = non[it].Gamma.min_tau_index_in_C; z < non[it].Gamma.alg.size(); z++) {
                        if (non[it].Gamma.alg[z].tau < now_min_tau_in_C)
                            non[it].Gamma.min_tau_index_in_C++;
                        if (non[it].Gamma.alg[z].tau >= now_max_tau_in_C) {
                            non[it].Gamma.max_tau_index_in_C = z;
                            need_update = false;
                            break;
                        }
                    }
                    if (need_update) {
                        //finally, go through all gamma now, put new pair (S,gamma) if needed
                        for (int t = left; t <= right; t++) {
                            double tau_temp = pow(2.0, t);
                            if (tau_temp > non[it].Gamma.max_tau)
                                non[it].Gamma.alg.push_back(ALG_known_tau(tau_temp));
                        }
                        non[it].Gamma.max_tau_index_in_C = non[it].Gamma.alg.size() - 1;//max gamma index in C now
                        non[it].Gamma.max_tau = non[it].Gamma.alg.back().tau;//the last element always is the max_gamma anyway
                    }
                }
            }
            bool u_is_choosen=false;
            for (int b = non[it].Gamma.min_tau_index_in_C; b <= non[it].Gamma.max_tau_index_in_C; b++)
            {
                //cout<<"tau: "<<Gamma[a].max_tau_index_in_C-Gamma[a].min_tau_index_in_C<<endl;
                //int ell=ceil(4.0+2.0*log2(k)+1.0+2.0*log2(G.Set.size()));
                non[it].ell = ceil(3.0 + 2.0 * log2(k) + 2.0 * log2(non[it].G.solution.size()));
                /*****start form S_0, so the number of the S is {0,1,2,...,ell}=ell+1******/
                if (non[it].ell + 1 > non[it].Gamma.alg[b].E.size()) {
                    for (int i = non[it].Gamma.alg[b].E.size(); i <= non[it].ell; i++) {
                        non[it].Gamma.alg[b].E.push_back(S_class());
                    }
                }
                /*****calculate cup marginal gain and send to m(u)******/
                double sum_weight=0.0;
                //bool circle=false;
                /****check every end_vertex v in N(u) whether has been covered by any solution******/
                for(const auto &v:Edge[u])
                {
                    bool is_covered=false;
                    for(int x=0;x<non[it].ell;x++)
                    {
                        if(non[it].Gamma.alg[b].E[x].node_covered.count(v))
                        {
                            is_covered=true;
                            break;
                        }
                    }
                    //if not covered, then we can sum weight as marginal gain
                    if(!is_covered)
                        sum_weight+=Groundset[v].weight;

/*                    //mark whether has circle
                    if(v==u)
                        circle=true;*/
                }
/*                *//****check u whether has been covered by any solution******//*
                bool u_is_covered=false;
                for(int x=0;x<non[it].ell;x++)
                {
                    if(non[it].Gamma.alg[b].E[x].node_covered.count(u))
                    {
                        u_is_covered=true;
                        break;
                    }
                }
                if(!u_is_covered&&!circle)
                    sum_weight+=Groundset[u].weight;*/

                double m_u=sum_weight-Groundset[u].cost;
                /*****calculate cup marginal gain******/

                oracle_times++;
                long long int i_u;
                if (m_u > 0) {
                    i_u = floor(log2(non[it].Gamma.alg[b].tau / m_u));
                } else {
                    i_u = 999999999999;
                }
                if (i_u >= 0 && i_u <= non[it].ell)
                {
                    if (non[it].Gamma.alg[b].E[i_u].is_feasible(u, max_node))
                    {
                        non[it].Gamma.alg[b].E[i_u].add_element(0.0, u);
                        u_is_choosen=true;
                    }
                }
            }
            if(u_is_choosen)
                break;
        }
    }

    /****go through all available tau*****/
    S_class S_best;
    for(int it=0;it<r;it++)
    {
//        cout << "max_tau_index: " <<non[it].Gamma.max_tau_index_in_C << endl;
//        cout << "min_tau_index: " << non[it].Gamma.min_tau_index_in_C << endl;
        for (int b = non[it].Gamma.min_tau_index_in_C; b <= non[it].Gamma.max_tau_index_in_C; b++)
        {
            //cout << "now_tau_index: " << b << endl;
            for (int j = 0; j <= h - 1; j++)
            {
                int i = j;
                S_class T_j;
                while (i <= non[it].ell)
                {
                    for (const auto &u:non[it].Gamma.alg[b].E[i].solution)
                    {
                        if (T_j.is_feasible(u, max_node)) {
                            //we only calculate the revenue whenever ends
                            T_j.add_element(0.0, u);
                        }
                    }
                    i += h;
                }
                T_j.S_revenue = T_j.f_S();
                oracle_times++;
                if (T_j.S_revenue > S_best.S_revenue)
                    S_best = T_j;
            }
        }
    }
    /*****call offline for every r and gamma *******/
    long long int saved_element=0;
    for(int it=0;it<r;it++)
    {
      //  cout<<"now r: "<<r<<endl;
       // cout<<"non[it].Gamma.max_tau_index_in_C: "<<non[it].Gamma.max_tau_index_in_C<<endl;
       // cout<<"non[it].Gamma.min_tau_index_in_C: "<<non[it].Gamma.min_tau_index_in_C<<endl;
        for (int b = non[it].Gamma.min_tau_index_in_C; b <= non[it].Gamma.max_tau_index_in_C; b++)
        {
         //   cout<<"now gamma id: "<<b<<endl;
            vector<int> A_it;
            for (int i = 0; i <= non[it].ell; i++)
            {
                //sum element/memory
                /***********/
                saved_element+=non[it].Gamma.alg[b].E[i].solution.size();
                /***********/

                /****set insert to vector for every r and every gamma****/
                A_it.insert(A_it.end(), non[it].Gamma.alg[b].E[i].solution.begin(), non[it].Gamma.alg[b].E[i].solution.end());
            }
        //    cout<<"call offline groundset size: "<<A_it.size()<<endl;

            //sum element/memory
            int offline_memory=0;

            pair<S_class,long long int> temp=multi_random_acc(eps,A_it,max_node,offline_memory);
            oracle_times+=temp.second;
            if (temp.first.S_revenue > S_best.S_revenue)
                S_best = temp.first;

            saved_element+=offline_memory;
        }
    }
    /*****call offline for every r*******/
/*    for(int it=0;it<r;it++)
    {
        vector<int> A_it;
        for (int b = non[it].Gamma.min_tau_index_in_C; b <= non[it].Gamma.max_tau_index_in_C; b++)
        {
            for (int i = 0; i <= non[it].ell; i++)
            {
                A_it.insert(A_it.end(), non[it].Gamma.alg[b].E[i].solution.begin(), non[it].Gamma.alg[b].E[i].solution.end());
            }
        }
//        cout<<"A size: "<<A_it.size()<<endl;
        pair<S_class,long long int> temp=multi_random_acc(eps,A_it,max_node);
        oracle_times+=temp.second;
        if (temp.first.S_revenue > S_best.S_revenue)
            S_best = temp.first;
    }*/
//    cout<<"Haba & Max node: "<<max_node<<endl;
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<S_best.S_revenue<<" cost: "<<S_best.S_cost<<" size: "<<S_best.solution.size()<<endl;
//    cout<<"Node: "<<endl;
//    for(const auto &p:S_best.solution)
//        cout<<p<<'\t';
//    cout<<endl;
//    cout<<" group limit: "<<endl;
//    for(int j=0;j<S_best.group_exist.size();j++)
//    {
//        cout<<S_best.group_exist[j]<<" ";
//    }
//    cout<<endl;
//   cout<<"memory: "<<saved_element<<endl;

    cout<<"Objective Values: "<<S_best.S_revenue<<endl;
    cout<<"Oracle Queries: "<<oracle_times<<endl;

    cout<<"Haba ---------end--------- "<<endl<<endl;
    return Result(S_best.S_revenue,S_best.S_cost,S_best.solution.size(),oracle_times,saved_element);
}


#endif //COVERAGE_HABA_NON_H
