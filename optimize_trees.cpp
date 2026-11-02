// Exact optimization over rooted-tree polynomial states, truncated at degree 3.
// Build: g++ -O3 -std=c++17 optimize_trees.cpp -o optimize_trees
// Run:   ./optimize_trees 151 > tree_optima.csv
// All integers fit uint64_t for the enforced range 2 <= N <= 301.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using U = std::uint64_t;
using Key = std::array<U,4>;
struct Hash {
    std::size_t operator()(const Key& q) const {
        std::size_t h=0;
        for (U x:q) h ^= std::hash<U>{}(x)+0x9e3779b9+(h<<6)+(h>>2);
        return h;
    }
};
struct Back { int j, f, t; };
struct State { Key q; Back back; };
std::vector<std::vector<State>> forest;
const std::vector<State> leaf={{{1,0,1,0},{0,0,0}}};

std::array<U,2> product(U a,U b,int pa,U c,U d,int pc) {
    if (pa+pc>=2) return {0,a*c};
    return {a*c,a*d+b*c};
}
Key attach(const Key& f,int r,const Key& t,int j) {
    auto a=product(f[0],f[1],r%2,t[2],t[3],j%2);
    auto b=product(f[2],f[3],1-r%2,t[2],t[3],j%2);
    auto c=product(f[0],f[1],r%2,t[0],t[1],1-j%2);
    return {a[0],a[1],b[0]+c[0],b[1]+c[1]};
}
std::string forest_word(int m,int i);
std::string tree_word(int n,int i) {
    return n==1 ? "()" : "("+forest_word(n-1,i)+")";
}
std::string forest_word(int m,int i) {
    if (m==0) return "";
    auto b=forest[m][i].back;
    return forest_word(m-b.j,b.f)+tree_word(b.j,b.t);
}

int main(int argc,char** argv) {
    int N=argc>1 ? std::stoi(argv[1]) : 151;
    if (N<2 || N>301) throw std::runtime_error("N must be in [2,301]");
    forest.push_back({{{1,0,0,0},{0,0,0}}});
    std::cout << "n,max_s3,states,witness_parentheses\n";
    std::cout << "2,0,1,(())\n";
    for (int m=1;m<=N-2;++m) {
        std::unordered_map<Key,Back,Hash> candidates;
        for (int j=1;j<=m;++j) {
            const auto& trees=j==1 ? leaf : forest[j-1];
            const auto& rest=forest[m-j];
            for (int fi=0;fi<(int)rest.size();++fi)
                for (int ti=0;ti<(int)trees.size();++ti)
                    candidates.emplace(attach(rest[fi].q,m-j,trees[ti].q,j),Back{j,fi,ti});
        }
        std::vector<State> sorted;
        sorted.reserve(candidates.size());
        for (const auto& kv:candidates) sorted.push_back({kv.first,kv.second});
        candidates.clear(); candidates.rehash(0);
        std::sort(sorted.begin(),sorted.end(),[](auto& a,auto& b){return a.q>b.q;});
        std::vector<State> keep;
        for (const auto& s:sorted) {
            bool dominated=false;
            for (const auto& k:keep) {
                if (k.q[0]>=s.q[0] && k.q[1]>=s.q[1] &&
                    k.q[2]>=s.q[2] && k.q[3]>=s.q[3]) {dominated=true;break;}
            }
            if (!dominated) keep.push_back(s);
        }
        forest.push_back(std::move(keep));
        int n=m+2,wi=0;
        U best=0;
        if(n%2) {
            for(int i=0;i<(int)forest[m].size();++i) {
                U v=forest[m][i].q[1]+forest[m][i].q[3];
                if(v>best){best=v;wi=i;}
            }
        }
        std::cout << n << ',' << best << ',' << forest[m].size() << ','
                  << '(' << tree_word(n-1,wi) << ")\n" << std::flush;
        std::cerr << n << " states=" << forest[m].size()
                  << " candidates=" << sorted.size() << '\n' << std::flush;
    }
}
