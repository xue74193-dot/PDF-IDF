#include "time.h"
#include "IterativeDistortedFilter.h"

int main(int argc,char *argv[]) {


   double eps=0.1;
    cout<<"eps: "<<eps<<endl;

    string filename = "../data/movies.txt";
    dataset *da = new dataset(filename);

    time_t nowtime;
    struct tm* p;;
    time(&nowtime);
    p = localtime(&nowtime);
    string outtext="../result/movie_result_"+to_string((int)ave_num)+"_"+to_string(p->tm_mon+1)+"."+to_string(p->tm_mday)+"_"+to_string(p->tm_hour)+"_"+to_string(p->tm_min)+"_"+to_string(p->tm_sec)+".txt";

    lambda_f=1.0;


    vector<Result> heuristic_random_result;
    vector<Result> distorted_result;
    vector<Result> iterative_result;
    vector<Result> multiran_result;
    vector<Result> haba_result;

    vector<int> ground_set;
    for(int i=0;i<node_num;i++)
        ground_set.push_back(i);

    int m_start=10;
    int m_end=30;
    int m_step=2;

    for(max_movie=m_start;max_movie<=m_end;max_movie+=m_step)
    {
        /**********calculate average for random algorithm************/
/*        double revenue_ave=0.0;
        long long int oracle_ave=0;
        int repeat_time=10;
        for(int i=0;i<repeat_time;i++) {
            Result temp_result = DistortedFiltering(eps,da);
            revenue_ave+=temp_result.revenue;
            oracle_ave+=temp_result.oracle;
        }
        Result pricing_result(revenue_ave/repeat_time,-1,-1,oracle_ave/repeat_time);
        distorted_result.emplace_back(pricing_result);*/
        /**********calculate average for random algorithm************/
      //  distorted_result.push_back(DistortedFiltering(eps,da));
      /*
      int R=50;
      vector<int> seeds;
      for (int i=0;i<R;i++)  seeds.push_back(1+i);

      Result best_result;
      double best_value=-1e18;

      for(int i=0;i<R;i++){
          int seed=seeds[i];
          Result res = MultiPassFilteringAlgorithm(eps,da,seed);
          if(res.revenue>best_value){
              best_value=res.revenue;
              best_result=res;
          }
      }
      cout<<"*********************************************************************************";
    */
        distorted_result.push_back(ParallelDistortedFilter(eps,da));
        iterative_result.push_back(IterativeDistortedFilter(eps,da));
        haba_result.push_back(Haba(eps,da));
        multiran_result.push_back(RandomMultiGreedy(eps,da,ground_set,max_movie));
        heuristic_random_result.push_back(RandomSelection(eps,max_movie,da));

    }

    ofstream out(outtext);
    out<<"eps: "<<eps<<" lambda_f:"<<lambda_f<<" lambda:"<<da->lambda<<" genres limit:"<<genres_limit<<endl;
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


    return 0;
}
