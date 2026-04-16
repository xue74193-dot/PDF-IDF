#ifndef PLA_STREAMINGHEURISTIC_H
#define PLA_STREAMINGHEURISTIC_H

#include "RandomMultiGreedy.h"

using namespace std;


Result RandomSelection(double eps, int max_edge)
{
    cout<<"RandomSelection & Max Edge: "<<max_edge<<"---------start---------"<<endl;
    S_class S;
    bernoulli_distribution ran(1.0/2.0);
    default_random_engine e(12345);
    //for(int u=0;u<node_num;u++)
    for(const auto &u:Groundset)
    {
        if(S.is_feasible(u,max_edge))
        {
            if(ran(e)==1)
                S.add_element(0.0,u);
        }
    }
    S.S_revenue=S.f_S();

//    cout<<"RandomSelection & Max Edge: "<<max_edge<<endl;
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<S.S_revenue<<" cost_1: "<<S.S_cost[0]<<" cost_1: "<<S.S_cost[1]<<" size: "<<S.edge_set.size()<<endl;
//    for(int i=0;i<S.edge_set.size();i++)
//        cout<<"("<<S.edge_set[i].x<<", "<<S.edge_set[i].y<<")"<<" ";
//    cout<<endl;
//    cout<<"  genre limit: "<<endl;
//    for(int j=0;j<S.genres_sum.size();j++)
//    {
//        cout<<S.genres_sum[j]<<" ";
//    }

    cout << "Objective Values: " << S.S_revenue << endl;
    cout << "RandomSelection ---------end--------- "<<endl<<endl;
    return Result(S.S_revenue,S.S_cost,S.edge_set.size(),0);
}


Result StreamingHeuristic(double eps, int max_edge)
{
    cout<<"StreamingHeuristic & Max Edge: "<<max_edge<<"---------start---------"<<endl;
    S_class S;
    bernoulli_distribution ran(1.0);
    default_random_engine e(12345);
   //for(int u=0;u<node_num;u++)
    for(const auto &u:Groundset)
    {
        if(S.is_feasible(u,max_edge))
        {
            if(ran(e)==1)
                S.add_element(0.0,u);
        }
    }
    S.S_revenue=S.f_S();

    cout<<"StreamingHeuristic & Max edge: "<<max_edge<<endl;
    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S.S_revenue<<" cost_0: "<<S.S_cost[0]<<" cost_1: "<<S.S_cost[1]<<" size: "<<S.edge_set.size()<<endl;
    for(int i=0; i<S.edge_set.size(); i++)
        cout<<"("<<S.edge_set[i].x<<", "<<S.edge_set[i].y<<")"<<" ";
    cout<<endl;

    cout<<endl;

    return Result(S.S_revenue,S.S_cost,S.edge_set.size(),0);
}

#endif //PLA_STREAMINGHEURISTIC_H
