#include "time.h"
#include "IterativeDistortedFilter.h"
int main(int argc,char *argv[]) {
    read_data();

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





    vector<Result> heuristic_random_result;
    vector<Result> distorted_result;
    vector<Result> iterative_result;
    vector<Result> multiran_result;
    vector<Result> haba_result;

    vector<int> ground_set;
    for(int i=0;i<node_num;i++)
        ground_set.push_back(i);

    double eps=0.1;
    int m_start=100;
    int m_end=1000;
    int m_step=100;


    for(int m=m_start;m<=m_end;m+=m_step)

    {
        cout<<"group limit: "<<group_limit<<endl;

        distorted_result.push_back(ParallelDistortedFilter(eps,m));
        iterative_result.push_back(IterativeDistortedFilter(eps, m));
        haba_result.push_back(Haba(eps,m));
        multiran_result.push_back(RandomMultiGreedy(eps,ground_set, m));
        heuristic_random_result.push_back(RandomSelection(eps,m));
    }

    ofstream out(outtext);
    out<<"eps: "<<eps<<endl;
    out<<"max node: "<<endl;
    for(int m=m_start;m<=m_end;m+=m_step)
    {
        out<<m<<"\t";
    }
    out<<endl;

    out<<"ParallelDistortedFilter "<<endl;
    out<<"revenue: "<<endl;
    for(auto p:distorted_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:distorted_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"memory: "<<endl;
    for(auto p:distorted_result)
    {
        out<<p.memory<<"\t";
    }
    out<<endl;


    out<<"IterativeDistortedFilter "<<endl;
    out<<"revenue: "<<endl;
    for(auto p: iterative_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p: iterative_result)
    {
        out<<p.oracle<<"\t";
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
    out<<"memory: "<<endl;
    for(auto p:haba_result)
    {
        out<<p.memory<<"\t";
    }
    out<<endl;

    out<<"RandomMultiGreedy "<<endl;
    out<<"revenue: "<<endl;
    for(auto p:multiran_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:multiran_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;

    out<<"RandomSelection "<<endl;
    out<<"revenue: "<<endl;
    for(auto p:heuristic_random_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:heuristic_random_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"memory: "<<endl;
    for(auto p:heuristic_random_result)
    {
        out<<p.size<<"\t";
    }
    out<<endl;


    return 0;
}