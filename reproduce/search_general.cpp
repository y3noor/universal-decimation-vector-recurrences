#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>
static int n,m,N,Q,tau;
static std::vector<int> ex,lg,ds;
static unsigned long long nodes[25],trials[25];
static int mul(int a,int b){return a&&b?ex[(lg[a]+lg[b])%tau]:0;}
static bool independent(const std::vector<int>& B,int d){
    int basis[32]={0};
    for(int j=0;j<m;j++)for(int a:B){
        int x=ex[(lg[a]+(long long)j*d)%tau];
        while(x){int k=31-__builtin_clz((unsigned)x);if(!basis[k]){basis[k]=x;break;}x^=basis[k];}
        if(!x)return false;
    }
    return true;
}
static bool admissible(const std::vector<int>& B){for(int d:ds)if(!independent(B,d))return false;return true;}
static void search(std::vector<int>& B,int last,int pivots){
    int r=B.size();nodes[r]++;
    if(r==n){
        bool isE=true;for(int x:B){int y=x;for(int j=0;j<n;j++)y=mul(y,y);if(y!=x)isE=false;}
        printf("GOOD is_E=%d basis=",isE);for(int x:B)printf(" %x",x);puts("");fflush(stdout);return;
    }
    for(int p=last+1;p<=N-(n-r);p++){
        std::vector<int> fp;for(int j=1;j<p;j++)if(!(pivots&(1<<j)))fp.push_back(j);
        for(int mask=0;mask<(1<<fp.size());mask++){
            int x=1<<p;for(size_t j=0;j<fp.size();j++)if(mask&(1<<j))x|=1<<fp[j];
            trials[r]++;B.push_back(x);
            if(admissible(B))search(B,p,pivots|(1<<p));
            B.pop_back();
        }
    }
}
int main(int argc,char**argv){
    if(argc!=4){puts("usage: search_general n m primitive_polynomial_hex");return 1;}
    n=atoi(argv[1]);m=atoi(argv[2]);N=n*m;if(n<2||m<2||N>24)return 2;
    Q=1<<N;tau=Q-1;int poly=strtol(argv[3],nullptr,16);
    ex.resize(tau);lg.assign(Q,-1);int x=1;
    for(int i=0;i<tau;i++){
        if(x==0||x>=Q||lg[x]!=-1){puts("invalid primitive polynomial");return 3;}
        ex[i]=x;lg[x]=i;if(std::gcd(i,tau)==1)ds.push_back(i);
        x<<=1;if(x&Q)x^=poly;
    }if(x!=1)return 4;
    auto t=std::chrono::steady_clock::now();std::vector<int>B{1};search(B,0,1);
    printf("n=%d m=%d polynomial=%x primitive_elements=%zu\n",n,m,poly,ds.size());
    for(int r=1;r<=n;r++)printf("dimension=%d surviving_partial_spaces=%llu extension_trials=%llu\n",r,nodes[r],trials[r]);
    printf("seconds=%.3f\n",std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count());
}
