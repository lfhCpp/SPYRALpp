#include "sim/Stopping.h"
#include "catima/storage.h"
#include <thread>
#include <vector>
#include <array>
#include <atomic>
#include <cmath>
#include <stdexcept>
int main(){
    GasTarget gas({{1,1,2}},600.);
    const int a[4]={1,2,4,12},z[4]={1,1,2,4};const double energy[4]={20,40,80,150};
    std::array<std::array<double,20>,4> expected;
    for(int p=0;p<4;++p){
        StoppingPower power;if(!power.Init(a[p],z[p],gas))return 1;
        for(int i=0;i<20;++i)expected[p][i]=power.EnergyAfter(energy[p],.05*i);
    }
    std::atomic<int> started{0};
    std::atomic<int> failures{0};std::vector<std::thread> pool;
    std::array<const void*,8> cache_addresses{};
    for(int worker=0;worker<8;++worker)pool.emplace_back([&,worker]{
        cache_addresses[worker]=&catima::_storage;
        ++started;while(started.load()<8)std::this_thread::yield();
        std::array<StoppingPower,4> powers;
        for(int p=0;p<4;++p)if(!powers[p].Init(a[p],z[p],gas))++failures;
        for(int round=0;round<100;++round)for(int p=0;p<4;++p)for(int i=0;i<20;++i){
            const double result=powers[p].EnergyAfter(energy[p],.05*i);
            if(!std::isfinite(result)||std::abs(result-expected[p][i])>1e-12)++failures;
        }
    });
    for(auto& t:pool)t.join();
    for(int i=0;i<8;++i)for(int j=i+1;j<8;++j)if(cache_addresses[i]==cache_addresses[j])++failures;
    if(failures)throw std::runtime_error("Threaded energy loss differs from sequential reference");
}
