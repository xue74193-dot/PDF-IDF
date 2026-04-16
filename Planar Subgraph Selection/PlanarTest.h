#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <vector>
using namespace std;
class Edge{
public:
    Edge(){};
    Edge(int a,int b,double w,const vector<double> &c,int i)
    {
        x=a;
        y=b;
        weight=w;
        cost=c;
        index=i;
    };
    int x;
    int y;
    double weight;
    vector<double> cost;
    //bool is_selected=false;
    int index;
};
bool test_planar(const int &n,const vector<Edge> &edge_set)
{
    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
    boost::property< boost::vertex_index_t, int > >
            graph;

    graph G(n);

    for(const auto &e:edge_set)
    {
        if (e.x == e.y) continue;
        add_edge(e.x, e.y, G);
    }
    return boyer_myrvold_planarity_test(G);
}
