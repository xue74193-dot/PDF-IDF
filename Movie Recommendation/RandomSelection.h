

#ifndef MOV_RANDOMSELECTION_H
#define MOV_RANDOMSELECTION_H
#include "RandomMultiGreedy.h"
Result RandomSelection(double eps,const int &max_node,dataset* da)
{
    cout<<"RandomSelection & Max node: "<<max_node<<"---------start---------"<<endl;
    S_gamma S;
    bernoulli_distribution ran(1.0/2.0);
    default_random_engine e(12345);
    for(int u=0;u<node_num;u++)
    {
        if(S.is_feasible(u,max_node))
        {
            if(ran(e)==1)
                S.add_element(0.0,u,da);
        }
    }
    S.S_revenue=S.f_S(da);

//    cout<<"RandomSelection & Max node: "<<max_node<<endl;
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<S.S_revenue<<" cost_1: "<<S.S_cost<<" size: "<<S.Set.size()<<endl;
//    for(int i=0;i<S.Set.size();i++)
//        cout<<S.Set[i]<<" ";
//    cout<<endl;
//    cout<<"  genre limit: "<<endl;
//    for(int j=0;j<S.genres_sum.size();j++)
//    {
//        cout<<S.genres_sum[j]<<" ";
//    }
//    cout<<endl;

    cout<<"Objective Values: " << S.S_revenue << endl;

    cout<<"RandomSelection ---------end--------- "<<endl<<endl;
    return Result(S.S_revenue,S.S_cost,S.Set.size(),0);
}

Result StreamingHeuristic(double eps,const int &max_node,dataset* da)
{
    cout<<"StreamingHeuristic & Max node: "<<max_node<<"---------start---------"<<endl;
    S_gamma S;
    bernoulli_distribution ran(1.0);
    default_random_engine e(12345);
    for(int u=0;u<node_num;u++)
    {
        if(S.is_feasible(u,max_node))
        {
            if(ran(e)==1)
                S.add_element(0.0,u,da);
        }
    }
    S.S_revenue=S.f_S(da);

    cout<<"StreamingHeuristic & Max node: "<<max_node<<endl;
    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S.S_revenue<<" cost_1: "<<S.S_cost<<" size: "<<S.Set.size()<<endl;
    for(int i=0;i<S.Set.size();i++)
        cout<<S.Set[i]<<" ";
    cout<<endl;
    cout<<"  genre limit: "<<endl;
    for(int j=0;j<S.genres_sum.size();j++)
    {
        cout<<S.genres_sum[j]<<" ";
    }
    cout<<endl;

    return Result(S.S_revenue,S.S_cost,S.Set.size(),0);
}
#endif //MOV_RANDOMSELECTION_H
