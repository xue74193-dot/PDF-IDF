#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <map>
#include <list>
using namespace std;

inline long long read(){
    long long ans = 0, f = 1;
    char ch = getchar();
    while(!isdigit(ch))
        f *= (ch == '-') ? -1 : 1, ch = getchar();
    do ans = (ans << 1) + (ans << 3) + (ch ^ 48), ch = getchar();
    while(isdigit(ch));
    return ans * f;
}
//limit the max number of node to be 205, while this should be changed
const int MAXN = 205;
int sta[MAXN], dfn[MAXN], low[MAXN], vis[MAXN], isEmbed[MAXN];
vector<int> Plane[MAXN<<1], book[MAXN<<1];
int PlaneNum = 1;

struct Graph{
    map<int, int> head;
    vector<int> next, last, val, att;
    int atp, atpPos;
    void clear(){
        head.clear(), next.clear(), last.clear(), val.clear(), att.clear(), atp = atpPos = 0;
        next.push_back(0), last.push_back(0), val.push_back(0);
        next.push_back(0), last.push_back(0), val.push_back(0);
    }
    Graph(){clear();}
    void add(int x,int y){
        next.push_back(head[x]), last.push_back(y), val.push_back(1), head[x] = next.size() - 1;
    }
    const bool operator < (const Graph &temp) const{
        return atp < temp.atp;
    }
}Tot;

void getAtp(Graph &G){
    sort(G.att.begin(), G.att.end()), G.atp = 0;
    for(int i=1; i<=PlaneNum; i++){
        if(book[i].size() < G.att.size()) continue;
        int now = 0;
        for(int j=0; j<G.att.size(); j++){
            while(now < book[i].size() - 1 && book[i][now] < G.att[j]) now++;
            if(book[i][now] != G.att[j]) break;
            else if(j == G.att.size() - 1)
                G.atp++, G.atpPos = i;
        }
    }
}

void embed(int pos){
    for(int i=1; i<=sta[0]; i++) isEmbed[sta[i]] = true;
    int l = 0, r = Plane[pos].size() - 1;
    while(Plane[pos][l] != sta[1] && Plane[pos][l] != sta[sta[0]]) l++;
    while(Plane[pos][r] != sta[1] && Plane[pos][r] != sta[sta[0]]) r--;
    vector<int> temp1, temp2;
    for(int i=0; i<l; i++) temp1.push_back(Plane[pos][i]);
    if(Plane[pos][l] == sta[1]) for(int i=1; i<=sta[0]; i++) temp1.push_back(sta[i]);
    else for(int i=sta[0]; i>=1; i--) temp1.push_back(sta[i]);
    for(int i=r+1; i<Plane[pos].size(); i++) temp1.push_back(Plane[pos][i]);
    for(int i=r-1; i>l; i--) temp2.push_back(Plane[pos][i]);
    if(Plane[pos][l] == sta[1]) for(int i=1; i<=sta[0]; i++) temp2.push_back(sta[i]);
    else for(int i=sta[0]; i>=1; i--) temp2.push_back(sta[i]);
    Plane[pos] = book[pos] = temp1, ++PlaneNum;
    Plane[PlaneNum] = book[PlaneNum] = temp2;
    sort(book[pos].begin(), book[pos].end()), sort(book[PlaneNum].begin(), book[PlaneNum].end());
}

bool match(int x,int goal,Graph &G){
    vis[x] = true;
    for(int l=G.head[x]; l; l=G.next[l]){
        int y = G.last[l];
        if(vis[y]) continue;
        if(y == goal || (!isEmbed[y] && match(y, goal, G))){
            G.val[l] = G.val[l^1] = 0;
            if(y == goal) sta[++sta[0]] = y;
            sta[++sta[0]] = x;
            return true;
        }
    }
    return false;
}

void findGraph(Graph &G,int l,Graph &ret){
    int x = G.last[l], fa = G.last[l^1];
    ret.add(x, fa), ret.add(fa, x), G.val[l] = G.val[l^1] = 0;
    if(!isEmbed[x]) for(int lk=G.head[x]; lk; lk=G.next[lk]){
            if(G.val[lk]) findGraph(G, lk, ret);
        }else if(!vis[x])
        ret.att.push_back(x), vis[x] = true;
}

bool Solve(list<Graph> &Lis){
    if(!Lis.size()) return true;
    list<Graph>::iterator it = Lis.begin();
    int cnt = Lis.size() - 1;
    while(!Lis.empty()){
        Graph &Now = *it;
        getAtp(Now), cnt++;
        if(!Now.atp) return false;
        if(cnt == Lis.size() || Now.atp == 1){
            memset(vis, 0, sizeof(vis));
            sta[0] = 0, match(Now.att[0], Now.att[1], Now);
            embed(Now.atpPos), memset(vis, 0, sizeof(vis));
            for(int j=2; j<sta[0]; j++) for(int l=Now.head[sta[j]]; l; l=Now.next[l]) if(Now.val[l]){
                        Graph temp;
                        findGraph(Now, l, temp);
                        if(!vis[sta[j]]) temp.att.push_back(sta[j]);
                        for(int k=0; k<temp.att.size(); k++) vis[temp.att[k]] = 0;
                        Lis.push_back(temp);
                    }
            list<Graph>::iterator temp = it++;
            Lis.erase(temp), cnt = 0, it--;
        }
        it++;
        if(it == Lis.end()) it = Lis.begin();
    }
    return true;
}

void Tarjan(int x,int fa,vector<Graph> &ret){
    dfn[x] = low[x] = ++dfn[0];
    for(int l=Tot.head[x]; l; l=Tot.next[l]){
        int y = Tot.last[l];
        if(y == fa) continue;
        if(!dfn[y]) Tarjan(y, x, ret), low[x] = min(low[x], low[y]);
        else low[x] = min(low[x], dfn[y]);
    }
    if(dfn[x] <= low[x]){
        Graph temp;
        for(int l=Tot.head[x]; l; l=Tot.next[l]) if(Tot.val[l] && dfn[Tot.last[l]] > dfn[x])
                findGraph(Tot, l, temp);
        ret.push_back(temp);
    }
}

void findCircle(Graph &G){
    int x = G.last[2];
    while(!vis[x]){
        vis[x] = true;
        for(int l=G.head[x]; l; l=G.next[l]) if((l ^ 1) != sta[sta[0]]){
                x = G.last[l], sta[++sta[0]] = l;
                break;
            }
    }
    int l = 1, r = sta[0];
    while(G.last[sta[l] ^ 1] != x) l++;
    sta[0] = 0;
    for(int i=l; i<=r; i++)
        G.val[sta[i]] = G.val[sta[i] ^ 1] = 0, sta[++sta[0]] = G.last[sta[i] ^ 1];
}

int main(){
    int T = read();
    while(T--)
    {
        //input the number of node and edge
        int n = read(), m = read();
        cout<<"node num: "<<n<<endl;
        cout<<"edge num: "<<m<<endl;
        vector<Graph> Div;
        //Empty the the graph structure, and then add a new graph
        Tot.clear();
        //read edge
        for(int i=1; i<=m; i++){
            int x = read(), y = read();
            if(x == y) continue;
            Tot.add(x, y), Tot.add(y, x);
            cout<<"edge "<<i<<" :"<<x<<" "<<y<<endl;
            //take as undirected graph, but there is no difference for planar graph?
        }
        //read Hamiltonian Circuit but not used
       // for(int i=1; i<=n; i++)
        //    read();
        if(m > 3 * n - 6 && m > 1){
            printf("NO\n");
            continue;
            //return 0;
        }
        //clean up
        memset(dfn, 0, sizeof(dfn));
        memset(low, 0, sizeof(low));
        memset(isEmbed, 0, sizeof(isEmbed));
        memset(vis, 0, sizeof(vis));
        for(int i=1; i<=n; i++) if(!dfn[i])
                Tarjan(i, -1, Div);
        bool flag = true;
        for(int i=0; i<Div.size(); i++){
            if(!Div[i].head.size()) continue;
            sta[0] = 0, findCircle(Div[i]);
            Plane[1].push_back(sta[1]), Plane[1].push_back(sta[sta[0]]);
            embed(1);
            list<Graph> ret;
            memset(vis, 0, sizeof(vis));
            for(int j=1; j<=sta[0]; j++) for(int l=Div[i].head[sta[j]]; l; l=Div[i].next[l]) if(Div[i].val[l]){
                        Graph temp;
                        findGraph(Div[i], l, temp);
                        if(!vis[sta[j]]) temp.att.push_back(sta[j]);
                        for(int k=0; k<temp.att.size(); k++) vis[temp.att[k]] = 0;
                        ret.push_back(temp);
                    }
            flag &= Solve(ret);
            for(int j=1; j<=PlaneNum; j++) Plane[j].clear(), book[j].clear();
            PlaneNum = 1;
            if(!flag) break;
        }
        if(flag) printf("YES\n");
        else printf("NO\n");
    }
    return 0;
}