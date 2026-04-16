
#ifndef COVERAGE_RANDOMSELECTION_H
#define COVERAGE_RANDOMSELECTION_H
#include "RandomMultiGreedy.h"
Result RandomSelection(double eps,const int &max_node)
{
    cout<<"RandomSelection & Max node: "<<max_node<<"---------start---------"<<endl;
    S_class S;
    bernoulli_distribution ran(1.0/2.0);
    default_random_engine e(12345);
    for(int u=0;u<node_num;u++)
    {
        if(S.is_feasible(u,max_node))
        {
            if(ran(e)==1)
                S.add_element(0.0,u);
        }
    }
    S.S_revenue=S.f_S();

//    cout<<"StreamingHeuristic & Max node: "<<max_node<<endl;
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<S.S_revenue<<" cost: "<<S.S_cost<<" size: "<<S.solution.size()<<endl;
//    for(const auto &p:S.solution)
//        cout<<p<<" ";
//    cout<<endl;
//    cout<<" group limit: "<<endl;
//    for(int j=0;j<S.group_exist.size();j++)
//    {
//        cout<<S.group_exist[j]<<" ";
//    }
//    cout<<endl;

    cout<<"Objective Values: "<<S.S_revenue<<endl;

    cout<<"RandomSelection ---------end--------- "<<endl<<endl;
    return Result(S.S_revenue,S.S_cost,S.solution.size(),0);
}


Result StreamingHeuristic(double eps,const int &max_node)
{
    cout<<"StreamingHeuristic & Max node: "<<max_node<<"---------start---------"<<endl;
    S_class S;
    bernoulli_distribution ran(1.0);
    default_random_engine e(12345);
    for(int u=0;u<node_num;u++)
    {
        if(S.is_feasible(u,max_node))
        {
            if(ran(e)==1)
                S.add_element(0.0,u);
        }
    }
    S.S_revenue=S.f_S();

    cout<<"StreamingHeuristic & Max node: "<<max_node<<endl;
    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S.S_revenue<<" cost: "<<S.S_cost<<" size: "<<S.solution.size()<<endl;
    for(const auto &p:S.solution)
        cout<<p<<" ";
    cout<<endl;
    cout<<" group limit: "<<endl;
    for(int j=0;j<S.group_exist.size();j++)
    {
        cout<<S.group_exist[j]<<" ";
    }
    cout<<endl;

    return Result(S.S_revenue,S.S_cost,S.solution.size(),0);
}

#endif //COVERAGE_RANDOMSELECTION_H
