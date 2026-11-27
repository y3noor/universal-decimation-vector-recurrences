#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>
// Exhaustive, early-pruned RREF enumeration for m=2, W containing 1.
// Each partial W must already have no primitive ratio. Every rejected partial
// space has a permanent witness, so pruning cannot remove a valid final W.
static int n,N,Q,tau;
static std::vector<int> ex,lg;
static std::vector<unsigned char> bad;
static uint64_t nodes[24],trials[24];
static int mul(int a,int b){return a&&b?ex[(lg[a]+lg[b])%tau]:0;}
static bool primitive_ratio(int a,int b){int d=lg[a]-lg[b];if(d<0)d+=tau;return bad[d];}
static bool extend_ok(const std::vector<int>& W,int x){
    // old/new pairs plus new/new pairs. A new/new pair has difference a in W*;
    // its primitive-ratio condition is t/(t+a) primitive, with t=x+s.
    for(int s:W){int t=x^s;
        for(size_t k=1;k<W.size();++k){int a=W[k];
            int d=lg[t]-lg[a];if(d<0)d+=tau;
            if(bad[d])return false;
            int u=ex[d]^1;
            if(u){int e=d-lg[u];if(e<0)e+=tau;if(bad[e])return false;}
        }
    }
    return true;
}
static void search(std::vector<int>& W,std::vector<int>& B,int last,int pivots){
    int r=B.size();nodes[r]++;
    if(r==n){
        bool isE=true;
        for(int x:B){int y=x;for(int j=0;j<n;j++)y=mul(y,y);if(y!=x)isE=false;}
        printf("GOOD is_E=%d basis=",isE);for(int x:B)printf(" %x",x);puts("");fflush(stdout);return;
    }
    for(int p=last+1;p<=N-(n-r);++p){
        std::vector<int> fp;for(int j=1;j<p;j++)if(!(pivots&(1<<j)))fp.push_back(j);
        for(int mask=0;mask<(1<<fp.size());++mask){
            int x=1<<p;for(size_t j=0;j<fp.size();j++)if(mask&(1<<j))x|=1<<fp[j];
            trials[r]++;if(!extend_ok(W,x))continue;
            size_t z=W.size();for(size_t j=0;j<z;j++)W.push_back(W[j]^x);
            B.push_back(x);search(W,B,p,pivots|(1<<p));B.pop_back();W.resize(z);
        }
    }
}
int main(int argc,char**argv){
    if(argc!=3){puts("usage: search_quadratic n primitive_polynomial_hex");return 1;}
    n=atoi(argv[1]);N=2*n;if(n<2||N>24)return 2;Q=1<<N;tau=Q-1;
    int poly=strtol(argv[2],nullptr,16);ex.resize(tau);lg.assign(Q,-1);bad.resize(tau);
    int x=1;for(int i=0;i<tau;i++){
        if(x==0||x>=Q||lg[x]!=-1){puts("invalid primitive polynomial");return 3;}
        ex[i]=x;lg[x]=i;bad[i]=std::gcd(i,tau)==1;
        x<<=1;if(x&Q)x^=poly;
    }if(x!=1)return 4;
    auto t=std::chrono::steady_clock::now();
    std::vector<int> W{0,1},B{1};search(W,B,0,1);
    printf("n=%d m=2 polynomial=%x\n",n,poly);
    for(int r=1;r<=n;r++)printf("dimension=%d surviving_partial_spaces=%llu extension_trials=%llu\n",r,(unsigned long long)nodes[r],(unsigned long long)trials[r]);
    printf("seconds=%.3f\n",std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count());
}
