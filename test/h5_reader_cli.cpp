#include "reco/H5TraceReader.h"
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=2)return 1;
    try {
        rk4::H5TraceReader reader;
        if(!reader.Open(argv[1]))return 1;
        rk4::TraceEvent event;
        while(reader.ReadNext(event))
            std::cout << event.number << ' ' << event.traces.size() << ' ' << event.ic_downscale << '\n';
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
}
