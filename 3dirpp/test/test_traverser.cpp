#include "type.h"
#include "tree.h"
#include "argument.h"
#include "body.h"
#include "tree.h"
#include "traverser.h"
#include <omp.h>

int main(int argc, char* argv[])
{
    rtfmm::title("rtfmm/3dirpp/test/test_traverser");
    rtfmm::Argument args(argc, argv);
    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    printf("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies bs = rtfmm::generate_random_bodies(args.n, args.r0, args.x0, args.seed, args.zero_netcharge);

    rtfmm::Tree tree(bs, args);
    std::cout<<"tree.cs.size() = "<<tree.cs.size()<<std::endl;
    std::cout<<"tree.depth_range = "<<tree.depth_range<<std::endl;

    rtfmm::Traverser traverser;
    traverser.traverse(tree);

    tbegin(fmm);
    // p2m
    #pragma omp parallel for
    for(int j = 0; j < traverser.leaf_cell_xmem.size(); j++)
    {
        rtfmm::Cell& c = tree.cs[traverser.leaf_cell_xmem[j]];
        int M = 0;
        for(int i = 0; i < c.brange[1]; i++)
        {
            M += 1;
        }
        c.M = M;
    }

    // m2m
    for(int depth = tree.depth_range[1]; depth >= tree.depth_range[0]; depth--)
    {
        std::vector<int>& depth_xmems = traverser.depth_xmem_map[depth];
        #pragma omp parallel for
        for(int j = 0; j < depth_xmems.size(); j++)
        {
            rtfmm::Cell& c = tree.cs[depth_xmems[j]];
            if(!c.is_leaf)
            {
                for(int octant = 0; octant < 8; octant++)
                {
                    int idx = tree.find_cidx(c.xlogic.child(octant));
                    rtfmm::Cell& child = tree.cs[idx];
                    c.M += child.M;
                }
            }
        }
    }

    // m2l
    #pragma omp parallel for
    for(int j = 0; j < traverser.m2l_map.size(); j++)
    {
        auto m2l = traverser.m2l_map.begin();
        std::advance(m2l, j);
        rtfmm::Cell& cj = tree.cs[m2l->first];
        std::vector<int>& srcs = m2l->second;
        for(int kj = 0; kj < 8; kj++)
        {
            rtfmm::LogicCoord lgcjj = cj.xlogic.child(kj);
            rtfmm::Cell& cjj = tree.cs[tree.find_cidx(lgcjj)];
            cjj.L = 0;
        }
        for(int i = 0; i < srcs.size(); i++)
        {
            rtfmm::Cell& ci = tree.cs[srcs[i]];
            // 8x8
            for(int kj = 0; kj < 8; kj++)
            {
                rtfmm::LogicCoord lgcjj = cj.xlogic.child(kj);
                rtfmm::Cell& cjj = tree.cs[tree.find_cidx(lgcjj)];
                for(int ki = 0; ki < 8; ki++)
                {
                    rtfmm::LogicCoord lgcii = ci.xlogic.child(ki);
                    rtfmm::Cell& cii = tree.cs[tree.find_cidx(lgcii)];
                    if(!rtfmm::LogicCoord::adjacent(lgcjj,lgcii))
                    {
                        cjj.L += cii.M;
                    }
                }
            }
        }
    }

    // p2p
    #pragma omp parallel for
    for(int j = 0; j < traverser.p2p_map.size(); j++)
    {
        auto p2p = traverser.p2p_map.begin();
        std::advance(p2p, j);
        rtfmm::Cell& cj = tree.cs[p2p->first];
        std::vector<int>& srcs = p2p->second;
        for(int i = 0; i < srcs.size(); i++)
        {
            rtfmm::Cell& ci = tree.cs[srcs[i]];
            for(int k = 0; k < cj.brange[1]; k++)
            {
                bs[cj.brange[0] + k].p += ci.brange[1];
            }
        }
    }

    // m2p
    #pragma omp parallel for
    for(int j = 0; j < traverser.m2p_map.size(); j++)
    {
        auto m2p = traverser.m2p_map.begin();
        std::advance(m2p, j);
        rtfmm::Cell& cj = tree.cs[m2p->first];
        std::vector<int>& srcs = m2p->second;
        for(int i = 0; i < srcs.size(); i++)
        {
            rtfmm::Cell& ci = tree.cs[srcs[i]];
            for(int k = 0; k < cj.brange[1]; k++)
            {
                bs[cj.brange[0] + k].p += ci.M;
            }
        }
    }    

    // p2l
    #pragma omp parallel for
    for(int j = 0; j < traverser.p2l_map.size(); j++)
    {
        auto p2l = traverser.p2l_map.begin();
        std::advance(p2l, j);
        rtfmm::Cell& cj = tree.cs[p2l->first];
        std::vector<int>& srcs = p2l->second;
        for(int i = 0; i < srcs.size(); i++)
        {
            rtfmm::Cell& ci = tree.cs[srcs[i]];
            cj.L += ci.brange[1];
        }
    }

    // l2l
    for(int depth = tree.depth_range[0]; depth <= tree.depth_range[1]; depth++)
    {
        std::vector<int>& depth_xmems = traverser.depth_xmem_map[depth];
        #pragma omp parallel for
        for(int j = 0; j < depth_xmems.size(); j++)
        {
            rtfmm::Cell& c = tree.cs[depth_xmems[j]];
            if(!c.is_leaf)
            {
                for(int octant = 0; octant < 8; octant++)
                {
                    int idx = tree.find_cidx(c.xlogic.child(octant));
                    rtfmm::Cell& child = tree.cs[idx];
                    child.L += c.L;
                }
            }
        }
    }

    // l2p
    #pragma omp parallel for
    for(int j = 0; j < traverser.leaf_cell_xmem.size(); j++)
    {
        rtfmm::Cell& c = tree.cs[traverser.leaf_cell_xmem[j]];
        for(int i = 0; i < c.brange[1]; i++)
        {
            int idx = c.brange[0] + i;
            bs[idx].p += c.L;
        }
    }

    tend(fmm);

    bs = rtfmm::sort_bodies_by_idx(bs);
    for(int j = 0; j < bs.size(); j++)
    {
        if(bs[j].p != bs.size())
        {
            std::cout<<"error : "<<bs[j]<<std::endl;
            exit(0);
        }
    }
    
    std::cout<<"test pass!"<<std::endl;

    return 0;
}