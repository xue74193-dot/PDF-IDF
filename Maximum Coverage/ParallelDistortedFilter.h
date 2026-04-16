
#ifndef COVERAGE_PARALLELDISTORTEDFILTER_H
#define COVERAGE_PARALLELDISTORTEDFILTER_H

#include "utilityfunction.h"
class ALG_known
{
public:
    double tau;
    vector<S_class> S;
    ALG_known(const double t)
    {
        tau=t;
    }
};
class ALG_Distorted
{
public:
    ALG_Distorted(){};
    vector<ALG_known> alg;
    double max_tau=0.0;//max tau in C, help to update tau
    int min_tau_index_in_C=-1;//-1 represent C is empty, else represent min tau index which is still available now
    int max_tau_index_in_C=-1;
};
Result ParallelDistortedFilter(double eps,const int &max_node)
{
    bernoulli_distribution ran(1.0);
    default_random_engine e(12345);

    cout<<"ParallelDistortedFilter & Max node: "<<max_node<<"---------start---------"<<endl;

    long long int oracle_times=0;

    double alpha_1=1.0+1.0/(k+1.0);
    //double alpha_1=1.0;
    double delta=3.0;
    S_class S_best;
    ALG_Distorted Gamma;
    S_class G;
    /*****only used for generate C, its first value seems useless?*****/
    pair<int,double> M=make_pair(-1,0.0);
    int counter=0;

    //for(const auto &u:Groundset)
    //for(int u=node_num-1;u>=0;u--)
    for(int u=0;u<node_num;u++)
    {
        if(counter%1000==0)
            //cout<<"Counter: "<<counter<<endl;
        counter++;
        /*******k-system constraint***********/
        if(G.is_feasible(u,max_node)){
            //marginal is useless for G, so we use a arbitrary value
            G.add_element(0.0,u);
        }
        /*****singleton must satisfy planar, so only check budget******/
        double fu_value = f_u(u);
        if (fu_value > M.second) {
            M.first = u;
            M.second = fu_value;
        }
        oracle_times++;
        if (fu_value > S_best.S_revenue) {
            S_best.replace_with_singleton(fu_value, u);
        }
        if(fu_value<=0) continue;

        //generate ALG2
        double left_temp=M.second;
        double right_temp=M.second*pow(1.0+delta,1.0+2.0*log(k)/log(1.0+delta)+2.0*log(G.solution.size())/log(1.0+delta) );

        if(left_temp>right_temp) continue;

        double left=ceil(log(left_temp)/log(1.0+delta));
        double right=floor(log(right_temp)/log(1.0+delta));

        if(Gamma.min_tau_index_in_C==-1)//first visit C
        {
            for(int t=left;t<=right;t++)
            {
                double tau_temp=pow(1.0+delta,t);
                Gamma.alg.push_back(ALG_known(tau_temp));
            }
            Gamma.min_tau_index_in_C=0;//the index of min tau which is still available now
            Gamma.max_tau=pow(1.0+delta,right);//now max tau
            Gamma.max_tau_index_in_C=Gamma.alg.size()-1;//max tau index in C now
        }
        else//not first visit
        {
            double now_min_tau_in_C=pow(1.0+delta,left);//min tau in C now
            double now_max_tau_in_C=pow(1.0+delta,right);//max tau in C now
            if(now_min_tau_in_C>Gamma.max_tau)//all old S_array should be removed
            {
                Gamma.min_tau_index_in_C=Gamma.alg.size();//the index of min tau which is still available now
                for(int t=left;t<=right;t++)//then add new S gamma pair
                {
                    double tau_temp=pow(1.0+delta,t);
                    Gamma.alg.push_back(ALG_known(tau_temp));
                }
                Gamma.max_tau_index_in_C=Gamma.alg.size()-1;//max gamma index in C now
                Gamma.max_tau = Gamma.alg.back().tau;//the last element always is the max_gamma anyway
            }
            else//else find where is the min gamma index now, which is equivalent to remove all S whose gamma < left
            {
                bool need_update= true;//judge S_array should be updated or not
                for (int z = Gamma.min_tau_index_in_C; z < Gamma.alg.size(); z++)
                {
                    if (Gamma.alg[z].tau < now_min_tau_in_C)
                        Gamma.min_tau_index_in_C++;
                    if (Gamma.alg[z].tau >= now_max_tau_in_C)
                    {
                        Gamma.max_tau_index_in_C=z;
                        need_update=false;
                        break;
                    }
                }
                if(need_update)
                {
                    //finally, go through all gamma now, put new pair (S,gamma) if needed
                    for (int t = left; t <= right; t++)
                    {
                        double tau_temp = pow(1.0+delta, t);
                        if (tau_temp > Gamma.max_tau)
                            Gamma.alg.push_back(ALG_known(tau_temp));
                    }
                    Gamma.max_tau_index_in_C=Gamma.alg.size()-1;//max gamma index in C now
                    Gamma.max_tau = Gamma.alg.back().tau;//the last element always is the max_gamma anyway
                }
            }
        }

        for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
        {
            //if F_max<=0, then continue;
            if(Gamma.alg[b].tau<=0) continue;

            int m=ceil(1.0+2.0*log(k)/log(1.0+delta)+2.0*log(G.solution.size())/log(1.0+delta));
            if(m>Gamma.alg[b].S.size())
            {
                for(int i=Gamma.alg[b].S.size();i<m;i++)
                {
                    Gamma.alg[b].S.push_back(S_class());
                }
            }
            for(int i=1;i<=m;i++)
            {
                double alpha_i=alpha_1;
                //double alpha_i=i*alpha_1;
                double threshold_i= Gamma.alg[b].tau/pow(1.0+delta,i-1.0);
                /*****calculate Cup distorted marginal profit******/
                double g_mariginal=0.0;
                //bool circle=false;
                /****check every end_vertex v in N(u) whether has been covered by any solution******/
                for(const auto &v:Edge[u])
                {
                    bool is_covered=false;
                    for(int x=0;x<i;x++)
                    {
                        if(Gamma.alg[b].S[x].node_covered.count(v))
                        {
                            is_covered=true;
                            break;
                        }
                    }
                    //if not covered, then we can sum weight as marginal gain
                    if(!is_covered)
                        g_mariginal+=Groundset[v].weight;

/*                    //mark whether has circle
                    if(v==u)
                        circle=true;*/
                }
/*                *//****check u whether has been covered by any solution******//*
                bool u_is_covered=false;
                for(int x=0;x<i;x++)
                {
                    if(Gamma.alg[b].S[x].node_covered.count(u))
                    {
                        u_is_covered=true;
                        break;
                    }
                }
                if(!u_is_covered&&!circle)
                    g_mariginal+=Groundset[u].weight;*/

                double cost=Groundset[u].cost;
                double distorted=g_mariginal-alpha_i*cost;
                /*****calculate Cup distorted marginal profit******/
                oracle_times++;

                bool flag=false;
                if(distorted>=threshold_i)
                {
                    flag=true;
                    if(Gamma.alg[b].S[i-1].is_feasible(u,max_node))
                    {
                        //the marginal gain input g()-c();
                        if(ran(e)==1)
                        {
                            double single_set_marginal=Gamma.alg[b].S[i - 1].marginal_gain(u);
                            Gamma.alg[b].S[i - 1].add_element(single_set_marginal, u);
                            //update S*
                            if(Gamma.alg[b].S[i - 1].S_revenue>S_best.S_revenue)
                            {
                                S_best=Gamma.alg[b].S[i - 1];
                            }
                        }

                        break;
                    }
                }

//                if(flag)
//                    break;
            }
        }
    }

//    cout<<"max_tau_index: "<<Gamma.max_tau_index_in_C<<endl;
//    cout<<"min_tau_index: "<<Gamma.min_tau_index_in_C<<endl;
    int m=ceil(1.0+2.0*log(k)/log(1.0+delta)+2.0*log(G.solution.size())/log(1.0+delta));
//    cout<<"m: "<<m<<endl;

/*    for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
    {
        cout<<"now tau: "<<Gamma.alg[b].tau<<endl;
        for(int i=0;i<m;i++)
        {
            cout<<"  revenue: "<<Gamma.alg[b].S[i].S_revenue<<"  Real revenue: "<<Gamma.alg[b].S[i].f_S()<<" cost: "<<Gamma.alg[b].S[i].S_cost<<" size: "<<Gamma.alg[b].S[i].solution.size()<<endl;
            for(const auto &p:Gamma.alg[b].S[i].solution)
                cout<<p<<" ";
            cout<<endl;
            cout<<" grounp limit: "<<endl;
            for(int j=0;j<Gamma.alg[b].S[i].group_exist.size();j++)
            {
                cout<<Gamma.alg[b].S[i].group_exist[j]<<" ";
            }
            cout<<endl;
        }
    }*/
    long long int saved_element=0;
    for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
    {
        for(const auto &p:Gamma.alg[b].S)
        {
            S_class temp_best=p;

            //sum element/memory
            saved_element+=temp_best.solution.size();

            for(const auto &it1:Gamma.alg[b].S) {
                for(const auto &it2:it1.solution)
                {
                    if(temp_best.is_feasible(it2,max_node))
                    {
                        temp_best.add_element(0.0,it2);
                    }
                }
            }
            temp_best.S_revenue=temp_best.f_S();

            if (temp_best.S_revenue > S_best.S_revenue)
                S_best = temp_best;
        }
    }

    /*****call offline for every gamma *******/
/*    for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
    {
        vector<int> U;
        for(const auto &p:Gamma.alg[b].S) {
            U.insert(U.end(), p.solution.begin(), p.solution.end());
        }
        pair<S_class,long long int> temp=Boost(eps,U,max_node,alpha_1,M);
        oracle_times+=temp.second;
        if (temp.first.S_revenue > S_best.S_revenue)
            S_best = temp.first;

//        cout<<"offline result: "<<endl;
//        cout<<"  revenue: "<<temp.first.S_revenue<<"  Real revenue: "<<temp.first.f_S()<<" cost: "<<temp.first.S_cost<<" size: "<<temp.first.solution.size()<<endl;
//        for(const auto &p:temp.first.solution)
//            cout<<p<<" ";
//        cout<<endl;
//        cout<<" grounp limit: "<<endl;
//        for(int j=0;j<temp.first.group_exist.size();j++)
//        {
//            cout<<temp.first.group_exist[j]<<" ";
//        }
//        cout<<endl;
    }*/

    /*****call offline for all existing sets *******/
/*    vector<int> U;
    for(int b=Gamma.min_tau_index_in_C;b<=Gamma.max_tau_index_in_C;b++)
    {
        for(const auto &p:Gamma.alg[b].S) {
            U.insert(U.end(), p.solution.begin(), p.solution.end());
        }
    }
    pair<S_class,long long int> temp=Boost(eps,U,max_node,alpha_1);
    oracle_times+=temp.second;
    if (temp.first.S_revenue > S_best.S_revenue)
        S_best = temp.first;*/


//    cout<<"DistortedFiltering & Max node: "<<max_node<<endl;
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<S_best.S_revenue<<" cost: "<<S_best.S_cost<<" size: "<<S_best.solution.size()<<endl;
//    for(const auto &p:S_best.solution)
//        cout<<p<<" ";
//    cout<<endl;
//    cout<<" group limit: "<<endl;
//    for(int j=0;j<S_best.group_exist.size();j++)
//    {
//        cout<<S_best.group_exist[j]<<" ";
//    }
//    cout<<endl;
//    cout<<"delta: "<<delta<<endl;
//    cout<<"memory: "<<saved_element<<endl;

    cout<<"Objective Values: " << S_best.S_revenue << endl;
    cout<<"Oracle Queries: " << oracle_times << endl;

    cout<<"ParallelDistortedFilter ----------------end---------------"<< endl << endl;
    return Result(S_best.S_revenue,S_best.S_cost,S_best.solution.size(),oracle_times,saved_element);
}
#endif //COVERAGE_PARALLELDISTORTEDFILTER_H
