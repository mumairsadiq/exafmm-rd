#include "type.h"
#include "tree.h"
#include "argument.h"
#include "body.h"
#include "tree.h"
#include "traverser.h"
#include <omp.h>
#include "kernel.h"
#include "surface.h"

int main(int argc, char* argv[])
{
    rtfmm::title("rtfmm/3dirpp/test/test_kernel");
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
    tbegin(traverse);
    traverser.traverse(tree);
    tend(traverse);

    // init cell
    #pragma omp parallel for
    for(int j = 0; j < traverser.leaf_cell_xmem.size(); j++)
    {
        rtfmm::Cell& c = tree.cs[traverser.leaf_cell_xmem[j]];
        for(int i = 0; i < c.brange[1]; i++)
        {
            c.bs.push_back(bs[c.brange[0] + i]);
        }   
    }
    #pragma omp parallel for
    for(int i = 0; i < tree.cs.size(); i++)
    {
        rtfmm::Cell& c = tree.cs[i];
        c.q_equiv = rtfmm::Matrix(rtfmm::get_surface_point_num(args.P), 1);
        c.p_check = rtfmm::Matrix(rtfmm::get_surface_point_num(args.P), 1);
    }

    rtfmm::LaplaceKernel kernel;
    tbegin(precompute);
    kernel.precompute(args.P, args.r0, args.images);
    kernel.precompute_m2l(args.P, args.r0, args.images);
    tend(precompute);

    tbegin(fmm);
    tbegin(p2m);
    // p2m
    #pragma omp parallel for
    for(int j = 0; j < traverser.leaf_cell_xmem.size(); j++)
    {
        rtfmm::Cell& c = tree.cs[traverser.leaf_cell_xmem[j]];
        kernel.p2m(args.P, c);
    }
    tend(p2m);

    // m2m
    tbegin(m2m);
    for(int depth = tree.depth_range[1]; depth >= tree.depth_range[0]; depth--)
    {
        std::vector<int>& depth_xmems = traverser.depth_xmem_map[depth];
        #pragma omp parallel for
        for(int j = 0; j < depth_xmems.size(); j++)
        {
            rtfmm::Cell& c = tree.cs[depth_xmems[j]];
            if(!c.is_leaf)
            {
                kernel.m2m(args.P, c, tree);
            }
        }
    }
    tend(m2m);

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
    tbegin(p2p);
    #pragma omp parallel for
    for(int j = 0; j < traverser.p2p_map.size(); j++)
    {
        auto p2p = traverser.p2p_map.begin();
        std::advance(p2p, j);
        rtfmm::Cell& ctar = tree.cs[p2p->first];
        kernel.p2p(tree.cs, p2p->second, ctar);
    }
    tend(p2p);

    // m2p
    tbegin(m2p);
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
            kernel.m2p(args.P, ci, cj);
        }
    }    
    tend(m2p);

    // p2l
    tbegin(p2l);
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
            kernel.p2l(args.P, ci, cj);
        }
    }
    tend(p2l);

    // l2l
    tbegin(l2l);
    for(int depth = tree.depth_range[0]; depth <= tree.depth_range[1]; depth++)
    {
        std::vector<int>& depth_xmems = traverser.depth_xmem_map[depth];
        #pragma omp parallel for
        for(int j = 0; j < depth_xmems.size(); j++)
        {
            rtfmm::Cell& c = tree.cs[depth_xmems[j]];
            if(!c.is_leaf)
            {
                kernel.l2l(args.P, c, tree);
            }
        }
    }
    tend(l2l);

    // l2p
    tbegin(l2p);
    #pragma omp parallel for
    for(int j = 0; j < traverser.leaf_cell_xmem.size(); j++)
    {
        rtfmm::Cell& c = tree.cs[traverser.leaf_cell_xmem[j]];
        kernel.l2p(args.P, c);
    }
    tend(l2p);

    tend(fmm);
    
    std::cout<<"test pass!"<<std::endl;

    return 0;
}