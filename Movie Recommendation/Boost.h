#ifndef MOV_BOOST_H
#define MOV_BOOST_H
#include "utility_functions.h"
#include "time.h"
#include "algorithm"
class Ai_for_boost
{
public:
    Ai_for_boost(int node_temp,double ratio)
    {
        node=node_temp;
        weight=f_u_distorted(node,ratio);
        l=1;
    }
    int node;
    double weight;
    int l;
};
bool cmp_for_boost(const Ai_for_boost &a,const Ai_for_boost &b)
{
    return a.weight>b.weight;
}
pair<int,double> GetMax_for_boost(long long int &oracle_times,double eps,vector<Ai_for_boost> &A,S_gamma &Si,const vector<int> &available,const int &max_movie,const double &ratio)
{
    double eta=1.0;
    double eps_apostrophe=eta*eps;
    for(auto e=A.begin();e!=A.end();)//remove unfeasible element
    {
        /*****if weight<=0, then delete; faster than not delete ??????****/
        if((available[(*e).node]==0)||(!Si.is_feasible((*e).node,max_movie))||(*e).weight<=0)
        {
            e=A.erase(e);
        }
        else
        {
            e++;
        }
    }
    //sort
    //*
    sort(A.begin(),A.end(),cmp_for_boost);
    double m=0.0;//always record weight of aij
    int ai=-1;
    for(auto &u:A)
    {
        double wt=u.weight;
        if(m>=wt) break;
        //update weight
        oracle_times++;

        u.weight=Si.marginal_distorted(u.node,ratio);
        if(u.weight>=(wt/(1+eps_apostrophe)))
        {
            ai=u.node;
            m=u.weight;
            break;
        }
        else
        {
            if(u.weight>m)
            {
                m=u.weight;
                ai=u.node;
            }
            u.l++;
        }
    }
    //if Aij is empty, then return nothing
    //if(aij==-1)
    //{
    //    return pair<int,double>(-1,0);
    //}
    for(auto e=A.begin();e!=A.end();)
    {
        if((*e).l>(log(node_num/eps_apostrophe)/log(1.0+eps_apostrophe)))
        {
            e=A.erase(e);
        }
        else
        {
            e++;
        }
    }
    return pair<int,double>(ai,m);
    //*/
}
pair<S_gamma,long long int> Boost(double eps,dataset* da,const vector<int> &N,const double &ratio) {
    long long int oracle_times = 0;
    vector<S_gamma> S;
    vector<vector<Ai_for_boost>> A;
    int ell = 1;

    for (int i = 0; i < ell; i++) {
        S.push_back(S_gamma());
        A.push_back(vector<Ai_for_boost>());
    }
    vector<int> available(node_num, 0);//mark S_i selected and discarded
    //initial Ai
    for (const auto &it:N) {
        //if(node_budget_feasible(it,B,data))
        //{
        oracle_times++;
        Ai_for_boost temp(it, ratio);
        for (int j = 0; j < ell; j++) {
            A[j].push_back(temp);
        }
        available[it] = 1;
        //}
    }
    while (1) {
        for (int j = 0; j < S.size(); j++) {
            S[j].max_marginal = -999999999.0;
            S[j].max_element = -1;
            //A=empty, then return element=-1 and marginal gain =0
            pair<int, double> temp = GetMax_for_boost(oracle_times, eps, A[j], S[j], available, max_movie, ratio);
            S[j].max_marginal = temp.second;
            S[j].max_element = temp.first;
        }
        double max_marginal = 0.0;
        int max_element = -1;
        int max_solution = -1;
        for (int j = 0; j < S.size(); j++) {
            if (S[j].max_marginal > max_marginal) {
                max_solution = j;
                max_marginal = S[j].max_marginal;
                max_element = S[j].max_element;
            }
        }
        if (max_solution == -1 || max_marginal <= 0)//non element or marginal gain<=0
            break;
        //S<-S\cup {u}
        double real_marginal = max_marginal + (ratio - 1.0) * modular_cost[max_element];
        S[max_solution].add_element(real_marginal, max_element, da);
        available[max_element] = 0;//discarded or selected
    }
    S_gamma *S_star;
    double max_revenue = -999999999.0;

    for (int i = 0; i < S.size(); i++) {
        if (S[i].S_revenue > max_revenue) {
            S_star = &S[i];
            max_revenue = S[i].S_revenue;
        }
    }
    return pair<S_gamma, long long int>(*S_star, oracle_times);
}
#endif //MOV_BOOST_H
