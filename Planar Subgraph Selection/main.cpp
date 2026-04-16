#include"utilityfunction.h"
#include "time.h"
#include "IterativeDistortedFilter.h"

int main(int argc,char *argv[]) {
    read_data();
    cout<<"max_cost: "<<max_cost<<endl;
    cout<<"min_cost: "<<min_cost<<endl;
  //  k_for_haba=k+(max_cost/min_cost)*d;
    cout<<"k for haba "<<k_for_haba<<endl;

    time_t nowtime;
    struct tm *p;;
    time(&nowtime);
    p = localtime(&nowtime);
    string::size_type pos1, pos2, posend;
    pos1 = edge_text.find_last_of("/");
    pos2 = edge_text.rfind("/", pos1 - 1);
    posend = edge_text.find_last_not_of("/");
    string name1 = edge_text.substr(pos2 + 1, pos1 - pos2 - 1);
    string name2 = edge_text.substr(pos1 + 1, posend);
    string result_name = name1 + "_" + name2;
    //cout<<result_name<<endl;
    string outtext =
            "../result/result_" + result_name + "_" + to_string(p->tm_mon + 1) + "." + to_string(p->tm_mday) + "_" +
            to_string(p->tm_hour) + "_" + to_string(p->tm_min) + "_" + to_string(p->tm_sec) + ".txt";



    vector<Result> heuristicrandom_result;
    vector<Result> distorted_result;
    vector<Result> iterative_result;
    vector<Result> multiran_result;
    vector<Result> haba_result;


    vector<int> ground_set;
    for(int i=0;i<edge_num;i++)
        ground_set.push_back(i);


    double eps=0.1;

    int E_start=100;
    int E_end=300;
    int E_step=20;

    for(int E=E_start;E<=E_end;E+=E_step)
    {
         distorted_result.push_back(ParallelDistortedFilter(eps,E));
         iterative_result.push_back(IterativeDistortedFilter(eps, E));
         haba_result.push_back(Haba(eps,E));
         multiran_result.push_back(RandomMultiGreedy( eps, ground_set, E));
         heuristicrandom_result.push_back(RandomSelection(eps,E));
    }


    ofstream out(outtext);
    out<<"eps: "<<eps<<endl;
    //out<<"B: "<<endl;
    for(int E=E_start;E<=E_end;E+=E_step)
    {
        out<<E<<"\t";
    }
    out<<endl;

    out<<"ParallelDistortedFilter "<<endl;
    out<<"revenue: "<<endl;
    for(auto p: distorted_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p: distorted_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"cost: "<<endl;
    for(auto p: distorted_result)
    {
        out<<p.cost[0]<<" "<<p.cost[1]<<"\t";
    }

    out<<endl;

    out<<"IterativeDistortedFilter "<<endl;
    out<<"revenue: "<<endl;
    for(auto p:iterative_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:iterative_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"cost: "<<endl;
    for(auto p: iterative_result)
    {
        out<<p.cost[0]<<" "<<p.cost[1]<<"\t";
    }
    out<<endl;

    out<<"Haba "<<endl;
    out<<"revenue: "<<endl;
    for(auto p:haba_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:haba_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"cost: "<<endl;
    for(auto p:haba_result)
    {
        out<<p.cost[0]<<" "<<p.cost[1]<<"\t";
    }
    out<<endl;

    out<<"RandomMultiGreedy "<<endl;
    out<<"revenue: "<<endl;
    for(auto p: multiran_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p: multiran_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"cost: "<<endl;
    for(auto p:  multiran_result)
    {
        out<<p.cost[0]<<" "<<p.cost[1]<<"\t";
    }
    out<<endl;

    out<<"RandomSelection "<<endl;
    out<<"revenue: "<<endl;
    for(auto p: heuristicrandom_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p: heuristicrandom_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"cost: "<<endl;
    for(auto p:  heuristicrandom_result)
    {
        out<<p.cost[0]<<" "<<p.cost[1]<<"\t";
    }
    out<<endl;

    return 0;
}