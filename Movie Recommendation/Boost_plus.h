#include "map"
#include "utility_functions.h"
#include "time.h"
#include "algorithm"
class Ai_for_boost
{
public:
    Ai_for_boost(int node_temp,double ratio)
    {
        node=node_temp;
        //weight=f_u_distorted(node,ratio);
        l=1;
    }
    int node;
    //double weight;
    int l;
};
//bool cmp_for_boost(const Ai_for_boost &a,const Ai_for_boost &b)
//{
//    return a.weight>b.weight;
//}
pair<int,double> GetMax_for_boost(long long int &oracle_times,double eps,multimap<double,Ai_for_boost> &A,S_gamma &Si,const vector<int> &available,const int &max_node,const double &ratio)
{
/********check this when pop a node*************/
    double m=0.0;//always record weight of aij
    int ai=-1;
    while(!A.empty())
    {
        //check the end node, i.e., the node with biggest weight
        auto it=A.end();
        it--;

        double old_value=it->first;
        //if the value of maximum element is less than 0, then all element has value less than 0 due to the diminish property of submodularity, so we break and return empty element.
        if(old_value<=0)
        {
            A.clear();
            ai=-1;
            m=0.0;
            break;
        }
        //if not satisfy the k-system constraint or has been selected by other solution, then delete it and pop the next element
        if(!Si.is_feasible(it->second.node,max_node)||available[it->second.node]==0)
        {
            A.erase(it);
            continue;
        }
        //if the element is useful, then we compute its new weight, and let its update numbers +1
        oracle_times++;
        double new_value=Si.marginal_distorted(it->second.node,ratio);
        //the value of element in multimap can be change directly, but the key can not be change, we need to erase it and re-insert it
        it->second.l++;
        //if the value of the element diminishes not much, we can return it
        if(new_value>=old_value/(1.0+eps))
        {
            ai=it->second.node;
            m=new_value;
            break;
        }
            //or we re-insert now element and then check the next element
        else
        {
            Ai_for_boost temp=it->second;
            //erase the element to update its weight
            A.erase(it);
            //if its update numbers is not greater than max, then we re-insert it into the queue
            if(temp.l<=(log(node_num*2.0/eps)/log(1.0+eps))) {
                A.insert(pair<double, Ai_for_boost>(new_value, temp));
            }
        }
    }
    return pair<int,double>(ai,m);
}
pair<S_gamma,long long int> Boost(double eps,dataset* da,const vector<int> &N,const double &ratio,pair<int,double> u_star)
{
    long long int oracle_times=0;
    vector<S_gamma> S;
    vector<multimap<double,Ai_for_boost>> A;
    int ell=1;

    for(int i=0;i<ell;i++) {
        S.push_back(S_gamma());
        A.push_back(multimap<double,Ai_for_boost>());
    }
    vector<int> available(node_num,0);//mark S_i selected and discarded
    //initial Ai
    for(const auto &it:N)
    {
        /*******not add the u_star**********/
        if(it==u_star.second)
            continue;
        //if(node_budget_feasible(it,B,data))
        //{
        oracle_times++;
        Ai_for_boost temp(it,ratio);
        double value=f_u_distorted(it,ratio);

        for (int j = 0; j < ell; j++) {
            A[j].insert(pair<double,Ai_for_boost>(value,temp));
        }
        available[it]=1;
        //}
    }

    //S<-S\cup {u^*}
    S[0].add_element(u_star.second,u_star.first,da);
    available[u_star.first] = 0;//discarded or selected

    while(1)
    {
        for (int j = 0; j < S.size(); j++)
        {
            S[j].max_marginal=-999999999.0;
            S[j].max_element=-1;
            //A=empty, then return element=-1 and marginal gain =0
            pair<int,double> temp=GetMax_for_boost(oracle_times,eps,A[j],S[j],available,max_movie,ratio);
            S[j].max_marginal=temp.second;
            S[j].max_element=temp.first;
        }
        double max_marginal = 0.0;
        int max_element=-1;
        int max_solution = -1;
        for(int j=0;j<S.size();j++)
        {
            if (S[j].max_marginal > max_marginal) {
                max_solution=j;
                max_marginal=S[j].max_marginal;
                max_element=S[j].max_element;
            }
        }
        if (max_solution == -1||max_marginal<=0)//non element or marginal gain<=0
            break;
        //S<-S\cup {u}
        double real_marginal=max_marginal+(ratio-1.0)*modular_cost[max_element];
        S[max_solution].add_element(real_marginal,max_element,da);
        available[max_element] = 0;//discarded or selected
    }
    S_gamma *S_star;
    double max_revenue=-999999999.0;

    for(int i=0;i<S.size();i++)
    {
        if(S[i].S_revenue>max_revenue)
        {
            S_star=&S[i];
            max_revenue=S[i].S_revenue;
        }
    }
    return pair<S_gamma,long long int>(*S_star,oracle_times);
}

