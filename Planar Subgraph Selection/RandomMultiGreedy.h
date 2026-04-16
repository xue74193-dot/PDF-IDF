
#ifndef PLA_MULTIRANDOM_PLUS_H
#define PLA_MULTIRANDOM_PLUS_H
#include "ParallelDistortedFilter.h"
using namespace std;
class Ai_edge
{
public:
    Ai_edge(int edge_temp)
    {
        edge_intex=edge_temp;

        l=1;
    }
    int edge_intex;

    int l;
};

pair<int,double> GetMax(long long int &oracle_times,double eps,multimap<double,Ai_edge> &A,S_class &Si,const vector<int> &available,const int &max_edge)
{
/********check this when pop a node*************/
    double m=0.0; //always record weight of aij
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
        if(!Si.is_feasible(Groundset[it->second.edge_intex],max_edge)||available[it->second.edge_intex]==0)
        {
            A.erase(it);
            continue;
        }
        //if the element is useful, then we compute its new weight, and let its update numbers +1
        oracle_times++;
        double new_value=Si.marginal_gain_f(Groundset[it->second.edge_intex]);
        //the value of element in multimap can be change directly, but the key can not be change, we need to erase it and re-insert it
        it->second.l++;
        //if the value of the element diminishes not much, we can return it
        if(new_value>=old_value/(1.0+eps))
        {
            ai=it->second.edge_intex;
            m=new_value;
            break;
        }
            //or we re-insert now element and then check the next element
        else
        {
            Ai_edge temp=it->second;
            //erase the element to update its weight
            A.erase(it);
            //if its update numbers is not greater than max, then we re-insert it into the queue
            if(temp.l<=(log(edge_num*2.0/eps)/log(1.0+eps))) {
                A.insert(pair<double, Ai_edge>(new_value, temp));
            }
        }
    }
    return pair<int,double>(ai,m);
}
pair<S_class,long long int> multi_random_acc(double eps,const vector<int> &N,const int &max_node)
{
    double p=2.0/(1.0+sqrt((float)k));
    bernoulli_distribution u(p);
    default_random_engine e(100);

    long long int oracle_times=0;
    vector<S_class> S;

    vector<multimap<double,Ai_edge>> A;
    int ell=2;

    for(int i=0;i<ell;i++) {
        S.push_back(S_class());
        A.push_back(multimap<double,Ai_edge>());
    }
    vector<int> available(edge_num,0);//mark S_i selected and discarded
    //initial Ai
    for(const auto &it:N)
    {

        oracle_times++;

        Ai_edge temp(it);
        double value=f_u(Groundset[it]);

        for (int j = 0; j < ell; j++) {
            A[j].insert(pair<double,Ai_edge>(value,temp));
        }
        available[it]=1;

    }

    while(1)
    {
        for (int j = 0; j < S.size(); j++)
        {
            S[j].max_marginal=-999999999.0;
            S[j].max_element=-1;
            //A=empty, then return element=-1 and marginal gain =0
            pair<int,double> temp=GetMax(oracle_times,eps,A[j],S[j],available,max_node);
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
        if(u(e)==1)
        {
            S[max_solution].add_element(max_marginal,Groundset[max_element]);
        }
        available[max_element] = 0;//discarded or selected
    }

    S_class *S_star;
    double max_revenue=-999999999.0;

//    cout<<"MultiRandomAcc & m: "<<max_movie<<endl;
    for(int i=0;i<S.size();i++)
    {

        if(S[i].S_revenue>max_revenue)
        {
            S_star=&S[i];
            max_revenue=S[i].S_revenue;
        }
    }

    return pair<S_class,long long int>(*S_star,oracle_times);
}
Result RandomMultiGreedy(double eps,const vector<int> &N,const int &max_edge)
{
    cout<<"RandomMultiGreedy & Max edge: "<<max_edge<<"---------start---------"<<endl;

    //random
    double p=2.0/(1.0+sqrt((float)k));
    int ell=2;
    //determined
    // double p=1.0;
   // int ell=ceil(sqrt(k))+1;

    bernoulli_distribution u(p);
    default_random_engine e(100);

    long long int oracle_times=0;
    vector<S_class> S;
    vector<multimap<double,Ai_edge>> A;


    for(int i=0;i<ell;i++) {
        S.push_back(S_class());
        A.push_back(multimap<double,Ai_edge>());
    }
    int counter=0;

    vector<int> available(edge_num,0);//mark S_i selected and discarded
    //initial Ai
    for(const auto &it:N)
    {

        oracle_times++;

        Ai_edge temp(it);
        double value=f_u(Groundset[it]);

        for (int j = 0; j < ell; j++) {
            A[j].insert(pair<double,Ai_edge>(value,temp));
        }
        available[it]=1;

    }
    while(1)
    {
        if(counter%1==0)
        {
            for (int j = 0; j < S.size(); j++)
            {
                //cout<<"S_"<<j<<" size: "<<S[j].Set.size()<<endl;
            }
        }
        counter++;

        for (int j = 0; j < S.size(); j++)
        {
            S[j].max_marginal=-999999999.0;
            S[j].max_element=-1;
            //A=empty, then return element=-1 and marginal gain =0
            pair<int,double> temp=GetMax(oracle_times,eps,A[j],S[j],available,max_edge);
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
        if(u(e)==1)
        {
            S[max_solution].add_element(max_marginal,Groundset[max_element]);

        }
        available[max_element] = 0;//discarded or selected
    }

    S_class *S_star;
    double max_revenue=-999999999.0;


    for(int i=0;i<S.size();i++)
    {

        if(S[i].S_revenue>max_revenue)
        {
            S_star=&S[i];
            max_revenue=S[i].S_revenue;
        }
    }
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<(*S_star).S_revenue<<" size: "<<(*S_star).edge_set.size()<<" cost1: "<<(*S_star).S_cost[0]<<" cost2: "<<(*S_star).S_cost[1]<<endl;
//    cout << "Graph Edge: " << endl;
//    for (const auto &p: (*S_star).edge_set)
//        cout << "(" << p.x << "," << p.y << ")"  << p.cost[0]<< " ";
//    cout << endl;
//    cout << "Graph Node: " << endl;
//    for (const auto &p: (*S_star).node_set)
//        cout << p << " ";

    cout << "Objective Values: " << (*S_star).S_revenue << endl;
    cout << "Oracle Queries: " << oracle_times << endl;

    cout<<"RandomMultiGreedy ---------end--------- "<<endl<<endl;
    return Result((*S_star).S_revenue,(*S_star).S_cost,(*S_star).edge_set.size(),oracle_times);

}

#endif //PLA_MULTIRANDOM_PLUS_H
