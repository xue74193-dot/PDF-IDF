#ifndef STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
#define STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sys/types.h>
#include <sstream>
#include <climits>
#include <string>
#include <vector>
#include <cmath>
#include <unordered_set>
#include "list"
#include "random"
#include "PlanarTest.h"
using namespace std;
const int d=2;
const int d_1=0;

vector<float> modular_cost;

string edge_text="../Graph/road-euroroad-undirected.txt";
string node_text="../Graph/road-euroroad-node_weight.txt";

int node_num=0;
int edge_num=0;
 int k_for_haba = 3;
const int k=3;
double min_cost=999999999.0;
double max_cost=-1.0;
vector<Edge> Groundset;
vector<double> node_weight;

void read_data()
{

    ifstream in2(edge_text);
    int v1,v2;
    double cost1,cost2,weight;
    edge_num=0;
    unordered_set<int> node;
    while(!in2.eof())
    {
        in2>>v1>>v2>>cost1>>cost2>>weight;
        if (in2.fail())
            break;
        edge_num++;
        node.insert(v1);
        node.insert(v2);
        vector<double> cost={cost1,cost2};
        Groundset.emplace_back(v1,v2,weight,cost,edge_num-1);
//        double min_temp=min(cost1,cost2);
//        double max_temp=max(cost1,cost2);
//        if(min_temp<min_cost)
//            min_cost=min_temp;
//        if(max_temp>max_cost)
//            max_cost=max_temp;


         if(cost1<min_cost) min_cost=cost1;
         if(cost1>max_cost) max_cost=cost1;

    }
    node_num=node.size();
    cout<<"map node size: "<<node_num<<endl;
    cout<<"edge size: "<<edge_num<<endl;
    in2.close();

    ifstream in1(node_text);
    int node_temp;
    double weight_temp;
    while(!in1.eof())
    {
        in1>>node_temp>>weight_temp;
        if (in1.fail())
            break;
        //cout<<weight_temp<<endl;
        node_weight.emplace_back(weight_temp);
    }
    in1.close();

    // edge_cost: U1(1,10000), edge_weight: U2(100.0, 5000.0)

//   default_random_engine e(1235);
//    uniform_real_distribution<double> u1(1,10000);
//    uniform_real_distribution<double> u2(100.0, 5000.0);
//    double cost_2=0.0;
//    ofstream out("edge_weight_cost_new.txt");
//    for(auto p:Groundset)
//    {
//        out<<p.x<<'\t'<<p.y<<'\t'<<u1(e)<<'\t'<<cost_2<<'\t'<<u2(e)<<endl;
//    }
//    out.close();

//    node_weight: u3(100,5000)

//    uniform_real_distribution<double> u3(100,5000);
//    ofstream node_out("node_weight_new.txt");
//    for(int i=0; i<node_num; i++){
//        node_out << i<<"\t"<<u3(e)<<endl;
//    }
//    node_out.close();
//
//
//    int max_id = 0;
//    for(auto p : Groundset) {
//        if(p.x > max_id) max_id = p.x;
//        if(p.y > max_id) max_id = p.y;
//    }
//    cout << "Max Node ID: " << max_id << endl;
//    cout << "Node Set Size (node_num): " << node_num << endl;
}

bool node_budget_feasible(const Edge &edge,double Budget)
{
    /**********check single node Budget constraint**********/
//    for(int i=0;i<d;i++)
//    {
//        if(edge.cost[i]>Budget)
//        {
//            return false;
//        }
//    }
        if(edge.cost[d_1]>Budget)
        {
            return false;
        }
    return true;
/*    if(edge.cost>Budget)
        return false;
    else
        return true;*/
}

double g_u(const Edge &edge)
{
    //return edge.weight;
    return node_weight[edge.x]+node_weight[edge.y];
}

double f_u(const Edge &edge)
{
    //return edge.weight;
//    double cost=0;
//    for(int i=0; i<d; i++){
//        cost+=edge.cost[i];
//    }
    double cost=edge.cost[d_1];
    return node_weight[edge.x]+node_weight[edge.y]-cost;
}

double f_u_distorted(const Edge &edge,const double &ratio)
{
    //return edge.weight;
    double cost=0;
//    for(int i=0; i<d; i++){
//        cost+=edge.cost[i];
//    }
     cost=edge.cost[d_1];
    return node_weight[edge.x]+node_weight[edge.y]-ratio*cost;
}
double cal_singlton_cost(const Edge &edge)
{
    double u_sum_cost=0.0;
//    for(int r=0;r<d;r++)
//    {
//        u_sum_cost+=edge.cost[r];
//    }
    u_sum_cost+=edge.cost[d_1];
    return u_sum_cost;
}


class S_class {
public:
    S_class() {
        for (int i = 0; i < d; i++)
            S_cost.push_back(0.0);
        S_revenue = 0.0;
        //for(int i=0)
        is_selected.resize(edge_num,false);
        record_feasible.resize(edge_num,true);
        max_marginal=-999999999.0;
        max_element=-1;
    }
    vector<Edge> edge_set;
    unordered_set<int> node_set;
    vector<double> S_cost;
    double S_revenue;
    vector<bool> is_selected;
    vector<bool> record_feasible;

    double max_marginal;
    int max_element;

    void clear()
    {
        for (int i = 0; i < d; i++)
            S_cost[i]=0.0;
        S_revenue = 0.0;
        is_selected.resize(edge_num,false);
        record_feasible.resize(edge_num,true);
        edge_set.clear();
        node_set.clear();
    }
    void replace_with_singleton(const double &marginal,const Edge &e)
    {
        clear();
        add_element(marginal,e);
    }
    //1. value==false, last time we find it is feasible, then it will never be feasible.
    //2. value==true, last time is feasible, but this time we still need to check
    //vector<int> edge_feasible;

    double g_S()
    {
        double value=0.0;
        for(const auto &p:node_set)
        {
            value+=node_weight[p];
        }
        return value;
        //return (double)edge_set.size();
    }
    double f_S()
    {
        double value=0.0;
        double total_cost=0.0;

        for(const auto &p:node_set)
        {
            value+=node_weight[p];
        }

        for(const auto &e:edge_set)
        {
//            for(int i=0;i<d;i++){
//                total_cost+=e.cost[i];
//            }
            total_cost+=e.cost[d_1];
        }
        return value-total_cost;
        //return (double)edge_set.size();
    }
    double f_S_distorted()
    {
        double value=0.0;
        double total_cost=0.0;
        double ratio=1.0;

        for(const auto &p:node_set)
        {
            value+=node_weight[p];
        }

        for(const auto &e:edge_set)
        {
           /* for(int i=0;i<d;i++){
                total_cost+=e.cost[i];
            }*/
            total_cost+=e.cost[d_1];
        }
        return value-ratio*total_cost;
        //return (double)edge_set.size();
    }
    double f_S_sub_u(const Edge&e)
    {
        unordered_set<int> temp_node_set;
        for(const auto &p:edge_set)
        {
            if(p.index!=e.index)
            {
                temp_node_set.emplace(p.x);
                temp_node_set.emplace(p.y);
            }
        }
        double value=0.0;
        for(const auto &p:temp_node_set)
        {
            value+=node_weight[p];
        }
        return value;
    }

    double marginal_gain_g(const Edge &edge)
    {
        double value=0.0;

        if(!node_set.count(edge.x)) {
            value+=node_weight[edge.x];
        }
        if(!node_set.count(edge.y)) {
            value+=node_weight[edge.y];
        }

        return value;
    }
//    double modular_cost(const Edge &edge)
//    {    double cost=0;
//  //      for(int i=0; i<d; i++){
//   //         cost+=edge.cost[i];
//   //     }
//        cost+=edge.cost[d_1];
//        return cost;
//    }

    double marginal_gain_f(const Edge &edge)
    {
        double cost=0;
        double value=0.0;
        if(!node_set.count(edge.x)) {
            value+=node_weight[edge.x];
        }
        if(!node_set.count(edge.y)) {
            value+=node_weight[edge.y];
        }

//        for(int i=0; i<d; i++){
//            cost+=edge.cost[i];
//        }
        cost+=edge.cost[d_1];
      //  value=value-cost;
        return value-cost;
    }

    double marginal_gain_f_distorted(const Edge &edge,const double &ratio )
    {
        double cost=0;
        double value=0.0;
        if(!node_set.count(edge.x)) {
            value+=node_weight[edge.x];
        }
        if(!node_set.count(edge.y)) {
            value+=node_weight[edge.y];
        }
//        for(int i=0; i<d; i++){
//            cost+=edge.cost[i];
//        }
        cost+=edge.cost[d_1];
        value=value-ratio*cost;
        return value;
    }

    bool budget_feasible(const Edge &edge,double Budget)
    {
        /**********check Budget constraint**********/
//        for(int i=0;i<d;i++)
//        {
//            if(S_cost[i]+edge.cost[i]>Budget)
//            {
//                return false;
//            }
//        }

            if(S_cost[d_1]+edge.cost[d_1]>Budget)
            {
                return false;
            }
        return true;
    }
    bool is_feasible(const Edge &edge,int max_edge)
    {
        /*************check set size*************/
        if(edge_set.size()>=max_edge)
            return false;
        /**********have been find feasible, never be feasible again**********/
        if(!record_feasible[edge.index])
            return false;

        /**********check whether is selected**********/
        if(is_selected[edge.index])
            return false;
        /**********check k-system constraint**********/
        //calculate total number of node in the graph
        int new_node_num=0;
        if(!node_set.count(edge.x)) {
            new_node_num++;
        }
        if(!node_set.count(edge.y)) {
            new_node_num++;
        }
        edge_set.push_back(edge);
        bool flag=test_planar(node_set.size()+new_node_num,edge_set);
        edge_set.pop_back();

        if(!flag)
            record_feasible[edge.index]=false;

        return flag;
    }
    void add_element(const double &marginal,const Edge &e)
    {
        S_revenue+=marginal;
//        for(int r=0;r<d;r++)
//            S_cost[r]+=e.cost[r];
        S_cost[d_1]+=e.cost[d_1];
        node_set.emplace(e.x);
        node_set.emplace(e.y);
        edge_set.emplace_back(e);
        is_selected[e.index]=true;
    }
    void delete_element(const double &revenue,const Edge &e)
    {
        S_revenue=revenue;
//        for(int r=0;r<d;r++)
//            S_cost[r]-=e.cost[r];
        S_cost[d_1]-=e.cost[d_1];
        for(vector<Edge>::iterator p=edge_set.begin();p!=edge_set.end();)
        {
            if((*p).index==e.index)
            {
                p=edge_set.erase(p);
                break;
            }
            else{
                p++;
            }
        }
        node_set.clear();
        for(const auto &e:edge_set) {
            node_set.emplace(e.x);
            node_set.emplace(e.y);
        }
        is_selected[e.index]=false;

    }

    bool haba_feasible(const Edge &edge,int max_edge)
    {

        if(is_feasible(edge,max_edge))
            return true;
        else
            return false;
    }
/********The version of renuming the graph***********/
    /*   bool is_feasible(const Edge &edge)
       {
           *//**********check k-system constraint**********//*
        //calculate total number of node in the graph
        int new_node_num=0;
        int new_node_array[2];
        if(!node_set.count(edge.x)) {
            new_node_array[new_node_num]=edge.x;
            new_node_num++;
        }
        if(!node_set.count(edge.y)) {
            new_node_array[new_node_num]=edge.y;
            new_node_num++;
        }

        *//********renum the graph***********//*
        int new_num=1;
        for(const auto &p:node_set)
        {
            renum[p]=new_num;
            new_num++;
        }
        for(int i=0;i<new_node_num;i++)
        {
            renum[new_node_array[i]]=new_num;
            new_num++;
        }
        edge_set.push_back(edge);
        bool flag=is_planar(node_set.size()+new_node_num,edge_set);
        edge_set.pop_back();
        return flag;
    }*/

};
class Result
{
public:
    Result(){}
    Result(double rev,const vector<double> &cos,int siz,long long int ora)
    {
        revenue=rev;
        cost=cos;
        size=siz;
        oracle=ora;
    }
    double revenue;
    long long int oracle;
    long long int round;
    vector<double> cost;
    int size;
    long long int time=0;
    long long int max_query=0;
};
class SolutionWithStore
{
public:
    SolutionWithStore()
    {
        for(int i=0;i<edge_num;i++)
            is_planar.push_back(-1);
    };
    S_class solution;
    //1. is_planar=-1, we don't know
    //2. is_planar=0, not planar;
    //3. is_planar=1, is planar.
    vector<int> is_planar;
};

#endif //STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
