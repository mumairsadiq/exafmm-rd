#include "type.h"
#include "tree.h"
#include "argument.h"
#include "body.h"
#include "tree.h"

int main(int argc, char* argv[])
{
    rtfmm::title("rtfmm/3dirpp/test/test_tree");
    rtfmm::Argument args(argc, argv);
    args.show();

    /* prepare bodies */
    rtfmm::Bodies bs = rtfmm::generate_random_bodies(args.n, args.r0, args.x0, args.seed, args.zero_netcharge);

    rtfmm::Tree tree(bs, args);
    std::cout<<"tree.cs.size() = "<<tree.cs.size()<<std::endl;
    std::cout<<"tree.depth_range = "<<tree.depth_range<<std::endl;

    for(int i = 0; i < tree.cs.size(); i++)
        std::cout<<tree.cs[i]<<std::endl;
    std::cout<<tree.find_cidx(rtfmm::LogicCoord(-2,1,1,1))<<std::endl;
   
    int bcnt = 0;
    for(int i = 0; i < tree.cs.size(); i++)
    {
        rtfmm::Cell& c = tree.cs[i];
        if(c.is_leaf)
        {
            bcnt += c.brange[1];
        }
        for(int j = 0; j < c.brange[1]; j++)
        {
            int bidx = c.brange[0] + j;
            rtfmm::Body& b = bs[bidx];
            rtfmm::vec3r dx = (b.x - c.xphys).abs();
            if(dx[0] > c.r || dx[1] > c.r || dx[2] > c.r) 
            {
                std::cout<<"error! "<<std::endl;
                std::cout<<c<<std::endl;
                std::cout<<b<<std::endl;
                exit(1);
            }
        }
    }
    std::cout<<"bcnt = "<<bcnt<<std::endl;
    rtfmm::assert_exit(bcnt == bs.size(), "bcnt error!");

    std::cout<<"test pass!"<<std::endl;

    return 0;
}