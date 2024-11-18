#include "fmm.h"
#include "tree.h"
#include "mathfunc.h"
#include <omp.h>
#include "surface.h"
#include "direct/fmm.h"

rtfmm::LaplaceFMM::LaplaceFMM(const Bodies3& bs_, const Argument& args_) 
    : bs(bs_), args(args_)
{
    assert_exit(bs.size() == args.n, "LaplaceFMM init body size error");
}

rtfmm::Bodies3 rtfmm::LaplaceFMM::solve()
{
    /* build tree */
    tbegin(build_and_traverse);
    tbegin(build_tree);
    Tree tree;
    tree.build(bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    // tree.build(bs, args.x, args.r, 3, Tree::TreeType::uniform);
    cs = tree.get_cells();
    if(verbose && args.check_tree) check_tree(cs);
    tend(build_tree);

    /* traverse to get interaction list */
    tbegin(traverse);
    traverser.traverse(tree, args.cycle, args.images, args.P);
    cs = traverser.get_cells();
    if(verbose && args.check_tree) check_traverser(traverser);
    if(verbose && args.check_tree) check_cells(cs);
    tree_depth_range = get_min_max_depth(cs);
    if(verbose) std::cout<<"# of cells = "<<cs.size()<<std::endl;
    if(verbose) std::cout<<"tree_depth_range = "<<tree_depth_range<<std::endl;
    tend(traverse);
    tend(build_and_traverse);

    tbegin(init_cell_matrix_func);
    init_cell_matrix(cs);
    tend(init_cell_matrix_func);
    tbegin(init_reg_body_func);
    init_reg_body(cs);
    tend(init_reg_body_func);
    

    if(args.use_precompute)
    {
        TIME_BEGIN(precompute);
        TIME_BEGIN(precompute_others);
        kernel.precompute(args.P, args.r, args.images);
        TIME_END(precompute_others);
        TIME_BEGIN(precompute_m2l);
        kernel.precompute_m2l(args.P, args.r, cs, traverser.get_m2l_map(), args.images);
        TIME_END(precompute_m2l);
        if(args.timing) {TIME_END(precompute);}
    }

    tbegin(FMM_kernels);
    P2M();
    M2M();
    M2L();
    P2L();
    L2L();
    L2P();
    M2P();
    // P2P();


    TIME_BEGIN(fmm_p2p_standalone);
    std::vector<gmx::RVec> coordinates;
    std::vector<gmx::real> charges;

    for (int i = 0; i < bs.size(); i++)
    {
        coordinates.push_back(bs[i].x);
        charges.push_back(bs[i].q);
    }

    gmx::fmm::FMMDirectInteractions fmm_direct_inters(coordinates, charges, args.x, args.r, args.ncrit, args.rega);
    auto forces_and_potentials = fmm_direct_inters.execute_direct_kernel();

    for (int i = 0; i < bs.size(); i++)
    {
        bs[i].f += forces_and_potentials[i].first;
        bs[i].p += forces_and_potentials[i].second;
    }
    TIME_END(fmm_p2p_standalone);

    bs = sort_bodies_by_idx(bs);
    for(int i = 0; i < traverser.leaf_cell_idx.size(); i++)
    {
        int idx = traverser.leaf_cell_idx[i];
        Cell3& cell = cs[idx];
        for(int j = 0; j < cell.bodies.size(); j++)
        {
            Body3& b = cell.bodies[j];
            bs[b.idx].p += b.p;
            bs[b.idx].f += b.f;
        }
    }
    tend(FMM_kernels);

    if(args.dipole_correction)
        dipole_correction(bs, args.cycle);
        
    if(args.divide_4pi)
        scale_bodies(bs);
    //bs = sort_bodies_by_idx(bs);
    return bs;
}

void rtfmm::LaplaceFMM::P2M()
{
    Indices lcidx = get_leaf_cell_indices(cs);
    int leaf_number = lcidx.size();
    if(verbose) std::cout<<"P2M leaf_number = "<<leaf_number<<std::endl;
    TIME_BEGIN(P2M);
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        if(args.use_precompute)
            kernel.p2m_precompute(args.P, c);
        else
            kernel.p2m(args.P, bs, c);
    }
    if(args.timing)
    {
        TIME_END(P2M);
    }
}

void rtfmm::LaplaceFMM::M2M()
{
    TIME_BEGIN(M2M);
    for(int depth = tree_depth_range[1]; depth > tree_depth_range[0]; depth--)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            if(depth >= 0)
            {
                if(args.use_precompute)
                    kernel.m2m_precompute(args.P, ci, cs);
                else
                    kernel.m2m(args.P, ci, cs);
            }
            else
            {
                if(args.use_precompute)
                    kernel.m2m_img_precompute(args.P, ci, cs, args.cycle * std::pow(3, -depth - 1));
                else
                    kernel.m2m_img(args.P, ci, cs, args.cycle * std::pow(3, -depth - 1));
            }
        }
    }
    if(args.timing)
    {
        TIME_END(M2M);
    }
}

void rtfmm::LaplaceFMM::M2L()
{
    TIME_BEGIN(M2L);
    //PeriodicInteractionPairs m2l_pairs = traverser.get_pairs(OperatorType::M2L);
    //PeriodicInteractionMapM2L m2l_map = traverser.get_map(OperatorType::M2L);
    PeriodicM2LMap m2l_parent_map = traverser.get_M2L_parent_map();
    //PeriodicInteractionMapM2L m2l_map2 = traverser.get_m2l_map_from_m2l_parent_map();
    if(args.use_precompute)
    {
        //kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map2);
        //kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map);
        //kernel.m2l_fft_precompute_advanced3(args.P, cs, m2l_map, m2l_pairs);
        kernel.m2l_fft_precompute_t(args.P, cs, m2l_parent_map);
    }
    else
    {
        PeriodicInteractionMapM2L m2l_map = traverser.get_m2l_map();
        for(int i = 0; i < m2l_map.size(); i++)
        {
            auto m2l_list = m2l_map.begin();
            std::advance(m2l_list, i);
            Cell3& ctar = cs[m2l_list->first];
            for(int j = 0; j < m2l_list->second.size(); j++)
            {
                Cell3& csrc = cs[m2l_list->second[j].first];
                vec3r offset_src = m2l_list->second[j].second;
                if(args.use_fft) 
                {
                    kernel.m2l_fft(args.P, csrc, ctar, offset_src);
                }
                else
                {
                    kernel.m2l(args.P, csrc, ctar, offset_src);
                }
            }
        }
    }
    if(args.timing)
    {
        TIME_END(M2L);
    }
}

void rtfmm::LaplaceFMM::L2L()
{
    TIME_BEGIN(L2L);
    for(int depth = tree_depth_range[0]; depth < tree_depth_range[1]; depth++)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            if(args.use_precompute)
            {
                if(depth >= 0)
                    kernel.l2l_precompute(args.P, ci, cs);
                else
                    kernel.l2l_img_precompute(args.P, ci, cs);
            }
            else
                kernel.l2l(args.P, ci, cs);
        }
    }
    if(args.timing)
    {
        TIME_END(L2L);
    }
}

void rtfmm::LaplaceFMM::L2P()
{
    TIME_BEGIN(L2P);
    Indices lcidx = get_leaf_cell_indices(cs);
    int leaf_number = lcidx.size();
    if(verbose) std::cout<<"L2P leaf_number = "<<leaf_number<<std::endl;
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        if(args.use_precompute)
            kernel.l2p_precompute(args.P, c);
        else
            kernel.l2p(args.P, bs, c);
    }
    if(args.timing)
    {
        TIME_END(L2P);
    }
}

void rtfmm::LaplaceFMM::P2P()
{
    TIME_BEGIN(P2P);
    PeriodicInteractionMapP2P p2p_map = traverser.get_p2p_map();
    std::vector<int> leaves = traverser.leaf_cell_idx;
    #pragma omp parallel for
    for(int i = 0; i < p2p_map.size(); i++)
    {
        Cell3& ctar = cs[leaves[i]];
        if(args.use_simd)
        {
            /*std::vector<std::pair<int, vec3r>> temp;
            for(int j = 0; j < leaves.size(); j++)
            {
                temp.push_back(std::make_pair(leaves[j],vec3r(0,0,0)));
            }
            kernel.p2p_1toN_256(cs, temp, ctar);*/
            kernel.p2p_1toN_256(cs, p2p_map[i], ctar);
        }
        else
        {
            for(int j = 0; j < p2p_map[i].size(); j++)
            {
                Cell3& csrc = cs[p2p_map[i][j].first];
                kernel.p2p(bs,bs,csrc,ctar,p2p_map[i][j].second, 0);
            }
        }
    }
    if(args.timing)
    {
        TIME_END(P2P);
    }
}

void rtfmm::LaplaceFMM::M2P()
{
    TIME_BEGIN(M2P);
    PeriodicInteractionPairs m2p_pairs = traverser.get_pairs(OperatorType::M2P);
    for(auto m2p : m2p_pairs)
    {
        Cell3& ctar = cs[m2p.first];
        Cell3& csrc = cs[m2p.second.first];
        kernel.m2p(args.P, bs, csrc, ctar, m2p.second.second);
    }
    if(args.timing)
    {
        TIME_END(M2P);
    }
}

void rtfmm::LaplaceFMM::P2L()
{
    TIME_BEGIN(P2L);
    PeriodicInteractionPairs p2l_pairs = traverser.get_pairs(OperatorType::P2L);
    for(auto p2l : p2l_pairs)
    {
        Cell3& ctar = cs[p2l.first];
        Cell3& csrc = cs[p2l.second.first];
        kernel.p2l(args.P, bs, csrc, ctar, p2l.second.second);
    }
    if(args.timing)
    {
        TIME_END(P2L);
    }
}

void rtfmm::LaplaceFMM::init_cell_matrix(Cells3& cells)
{
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3& cell = cells[i];
        cell.q_equiv = Matrix(get_surface_point_num(args.P), 1);
        cell.p_check = Matrix(get_surface_point_num(args.P), 1);
    }
}

void rtfmm::LaplaceFMM::init_reg_body(Cells3& cells)
{
    std::vector<int> leaves = traverser.leaf_cell_idx;

    if(args.rega > 0)
    {
        printf("search 0,0,0\n");
        // search reg body in 0,0,0
        printf("~~~~~~~~~~ rega!\n");

        std::vector<Indices> boundary_bodies_idxs(leaves.size());
        // find weights for existing bodies within the same cell
        #pragma omp parallel for
        for(int leaf_idx = 0; leaf_idx < leaves.size(); leaf_idx++)
        {
            Cell3& cell = cells[leaves[leaf_idx]];
            cell.weights.resize(cell.brange.number);
            cell.bodies.resize(cell.brange.number);
            int rx = 0;
            for(int body_idx = cell.brange.offset; body_idx < cell.brange.offset + cell.brange.number; body_idx++)
            {
                const Body3 &body = bs[body_idx];
                const vec3r dx = body.x - cell.x;
                const vec3r dx_simcenter =
                    (body.x - args.x).abs() + args.rega;
                vec3r ws = get_w_xyz(dx, cell.r, args.rega);

                if(args.images == 0)
                {
                    for (int d = 0; d <= 2; d++)
                    {
                        if (dx_simcenter[d] >= args.r)
                        {
                            ws[d] = 1;
                        }
                    }
                }
                const real w = ws.mul();

                if (w < 1)
                {
                    boundary_bodies_idxs[leaf_idx].push_back(body_idx);
                }

                cell.bodies[rx] = body;
                cell.weights[rx++] = w;
            }
        }

        PeriodicInteractionMapP2P p2p_map = traverser.get_p2p_map();

        #pragma omp parallel for
        for(int leaf_idx = 0; leaf_idx < leaves.size(); leaf_idx++)
        {
            Cell3& cell = cells[leaves[leaf_idx]];
            const auto &aCells = p2p_map[cell.leaf_idx];
            

            // find out bodies from adjacent cells that are influencing this cell
            for (size_t j = 0; j < aCells.size(); j++)
            {
                auto& adj_cell_info = aCells[j];
                if (cell.idx != cells[adj_cell_info.first].idx)
                {
                    // check for regularization bodies from adjacent cells
                    for (const int &body_idx : boundary_bodies_idxs[cells[adj_cell_info.first].leaf_idx])
                    {

                        const Body3 &body = bs[body_idx];
                        const vec3r dx = body.x + adj_cell_info.second - cell.x;
                        const real w = get_w_xyz(dx, cell.r, args.rega).mul();

                        if (w > 0)
                        {
                            cell.reg_body_idx.push_back(std::make_pair(body_idx, adj_cell_info.second));
                        }
                    }
                }
            }
        }
        
        // compute weight and cache
        #pragma omp parallel for
        for(int leaf_idx = 0; leaf_idx < leaves.size(); leaf_idx++)
        {
            Cell3& cell = cells[leaves[leaf_idx]];
            for(int i = 0; i < cell.reg_body_idx.size(); i++)
            {
                int bidx = cell.reg_body_idx[i].first;
                Body3 body = bs[bidx];
                body.x += cell.reg_body_idx[i].second;
                vec3r dx = body.x - cell.x;
                vec3r dx_simcenter = (body.x - args.x).abs() + args.rega;
                vec3r ws = get_w_xyz(dx, cell.r, args.rega);
                if(args.images == 0)
                {
                    for(int d = 0; d <= 2; d++)
                    {
                        if(dx_simcenter[d] >= args.r)
                        {
                            ws[d] = 1;
                        }
                    }
                }
                real w = ws.mul();
                if(w > 0)
                {
                    cell.bodies.push_back(body);
                    cell.weights.push_back(w);
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for
        for(int leaf_idx = 0; leaf_idx < leaves.size(); leaf_idx++)
        {
            Cell3& cell = cells[leaves[leaf_idx]];
            for(int i = cell.brange.offset; i < cell.brange.offset + cell.brange.number; i++)
            {
                Body3 body = bs[i];
                cell.bodies.push_back(body);
                cell.weights.push_back(1);
                if (verbose)
                {
                    if(body.idx == 0)
                    {
                        std::cout<<body.x<<", "<<body.idx<<", "<<cell.idx<<std::endl;
                    }
                }
            }
        }
    }
}

rtfmm::vec2i rtfmm::LaplaceFMM::get_min_max_depth(const Cells3& cells)
{
    vec2i res = {cells[0].depth, cells[0].depth};
    for(int i = 0; i < cells.size(); i++)
    {
        res[0] = std::min(res[0], cells[i].depth);
        res[1] = std::max(res[1], cells[i].depth);
    }
    return res;
}

rtfmm::Indices rtfmm::LaplaceFMM::get_leaf_cell_indices(const Cells3& cells)
{
    Indices res;
    for(int i = 0; i < cells.size(); i++)
    {
        if(cells[i].crange.number == 0)
        {
            res.push_back(i);
        }
    }
    return res;
}

rtfmm::Indices rtfmm::LaplaceFMM::get_nonleaf_cell_indices(const Cells3& cells, int depth)
{
    Indices res;
    for(int i = 0; i < cells.size(); i++)
    {
        if(cells[i].depth == depth && cells[i].crange.number > 0)
        {
            res.push_back(i);
        }
    }
    return res;
}

void rtfmm::LaplaceFMM::check_tree(const Cells3& cells)
{
    printf("check tree\n");
    printf("cells.size() = %ld\n", cells.size());
    int bcnt = 0;
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3 c = cells[i];
        for(int j = 0; j < c.brange.number; j++)
        {
            Body3 b = bs[c.brange.offset + j];
            vec3r dx = (b.x - c.x).abs();
            assert_exit(dx[0] <= c.r && dx[1] <= c.r && dx[2] <= c.r, "tree body range error");
        }
        if(c.crange.number == 0)
        {
            bcnt += c.brange.number;
        }
    }
    assert_exit(bcnt == bs.size(), "tree body sum error");
}

void rtfmm::LaplaceFMM::check_traverser(Traverser& traverser)
{
    PeriodicInteractionPairs M2L_pairs = traverser.get_pairs(OperatorType::M2L);
    PeriodicInteractionPairs M2P_pairs = traverser.get_pairs(OperatorType::M2P);
    PeriodicInteractionPairs P2L_pairs = traverser.get_pairs(OperatorType::P2L);

    int m2l_num = M2L_pairs.size();
    int m2p_num = M2P_pairs.size();
    int p2l_num = P2L_pairs.size();

    if(verbose) std::cout<<"M2L_pair# : "<<m2l_num<<std::endl;
    if(verbose) std::cout<<"M2P_pair# : "<<m2p_num<<std::endl;
    if(verbose) std::cout<<"P2L_pair# : "<<p2l_num<<std::endl;
    if(verbose) std::cout<<"sum# : "<<m2l_num+m2p_num+p2l_num<<std::endl;

    // check pairs
    if(args.images == 0)
    {
        for(auto m2l : M2L_pairs)
        {
            int found = 0;
            for(auto m2l2 : M2L_pairs)
            {
                if(m2l.first == m2l2.second.first && m2l.second.first == m2l2.first)
                {
                    found++;
                }
            }
            assert_exit(found == 1, "treep m2l error");
        }
        std::cout<<"m2l check ok"<<std::endl;
    }

    for(auto m2p : M2P_pairs)
    {
        int found = 0;
        for(auto p2l : P2L_pairs)
        {
            if(m2p.first == p2l.second.first && m2p.second.first == p2l.first)
            {
                found++;
            }
        }
        assert_exit(found == 1, "treep p2l error");
    }
    std::cout<<"m2p check ok"<<std::endl;

    for(auto p2l : P2L_pairs)
    {
        int found = 0;
        for(auto m2p : M2P_pairs)
        {
            if(m2p.first == p2l.second.first && m2p.second.first == p2l.first)
            {
                found++;
            }
        }
        assert_exit(found == 1, "treep m2p error");
    }
    std::cout<<"p2l check ok"<<std::endl;
}

void rtfmm::LaplaceFMM::check_cells(const Cells3& cells)
{
    std::cout<<"check cells"<<std::endl;
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3 c = cells[i];
        assert_exit(c.idx == i, "cell idx error");
    }
}

rtfmm::real rtfmm::LaplaceFMM::reg_w(real x)
{
    return 0.25 * (2 + 3 * x - x * x * x);
}

rtfmm::real rtfmm::LaplaceFMM::get_w_single(real dx, real R, real rega)
{
    real r = std::abs(dx) - R + rega;
    real x;
    if(r <= 0) x = 1;
    else if(r > 2 * rega) x = -1;
    else x = 1 - r / rega;
    return reg_w(x);
}

rtfmm::real rtfmm::LaplaceFMM::get_w(vec3r dx, real R, real rega)
{
    real w = 1;
    for(int d = 0; d < 3; d++)
    {
        w *= get_w_single(dx[d], R, rega);
    }
    return w;
}

rtfmm::vec3r rtfmm::LaplaceFMM::get_w_xyz(vec3r dx, real R, real rega)
{
    rtfmm::vec3r res;
    for(int d = 0; d < 3; d++)
    {
        res[d] = get_w_single(dx[d], R, rega);
    }
    return res;
}