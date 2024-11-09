#include "epkifmm.h"
#include <iostream>
#include "ewald.h"
#include "kernel.h"

rtfmm::EpkiFMM::EpkiFMM(const Bodies3& bs_, const Argument& args_) : LaplaceFMM(bs_, args_)
{
    assert_exit(args.images == 1, "images should be 1");
    assert_exit(bs.size() == args.n, "EPKIFMM init body size error");
}

rtfmm::Bodies3 rtfmm::EpkiFMM::solve()
{
    /* build tree */
    tbegin(build_and_traverse);
    tbegin(build_tree);
    Tree tree;
    //tree.build(bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    tree.build(bs, args.x, args.r, 3, Tree::TreeType::uniform);
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

    tbegin(init_cell_matrix);
    init_cell_matrix(cs);
    init_reg_body(cs); // init body index in image-1 cells
    tend(init_cell_matrix);

    if(args.use_precompute)
    {
        TIME_BEGIN(precompute);
        TIME_BEGIN(precompute_others);
        kernel.precompute(args.P, args.r, args.images);
        TIME_END(precompute_others);
        TIME_BEGIN(precompute_m2l);
        kernel.precompute_m2l(args.P, args.r, cs, traverser.get_map(OperatorType::M2L), args.images);
        TIME_END(precompute_m2l);
        if(args.timing) {TIME_END(precompute);}
    }


    /* KIFMM upward */
    tbegin(EPKIFMM_kernels);
    P2M();
    M2M();
    M2L();

    /* calculate check surface of root cell */
    //calculate_root_check();
    //calculate_root_check2();
    //calculate_root_check3();
    //calculate_root_check4();
    //calculate_root_check5();
    calculate_root_check6();

    /* KIFMM downward */
    P2L();
    L2L();
    L2P();
    M2P();
    P2P();

    
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
    tend(EPKIFMM_kernels);

    if(args.dipole_correction)
        dipole_correction(bs, args.cycle);
        
    if(args.divide_4pi)
        scale_bodies(bs);
    //bs = sort_bodies_by_idx(bs);
    return bs;
}

void rtfmm::EpkiFMM::calculate_root_check()
{
    /* ewald[2,inf) part */
    Bodies3 ewald_bodies;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(args.P, cs[0].r * 1.05, cs[0].x);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 body;
        body.idx = i;
        body.p = 0;
        body.f = vec3r(0,0,0);
        body.x = x_equiv[i];
        body.q = cs[0].q_equiv[i];
        ewald_bodies.push_back(body);
    }
    Argument ewald_args = args;
    ewald_args.r = cs[0].r * 3;
    EwaldSolver ewald_solver(ewald_bodies, args);
    ewald_bodies = ewald_solver.solve();
    dipole_correction(ewald_bodies, args.cycle, 1);

    /* kifmm[0,1] part */
    LaplaceKernel kernel2;
    std::vector<real> root_p_check = kernel2.m2m_ewald_root0(args.P, cs[0], cs, args.cycle);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        cs[0].p_check[i] += ewald_bodies[i].p - root_p_check[i];
    }
}

void rtfmm::EpkiFMM::calculate_root_check2()
{
    printf("calculate_root_check2\n");
    /* direct[2,args.images] part */
    LaplaceKernel kernel1;
    Bodies3 direct_bodies;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(args.P, cs[0].r * 1.05, cs[0].x);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 b;
        b.idx = i;
        b.x = x_equiv[i];
        b.q = 0;
        b.p = 0;
        b.f = vec3r(0,0,0);
        direct_bodies.push_back(b);
    }
    Bodies3 direct_bodies2 = direct_bodies;

    kernel1.direct(bs, direct_bodies, 3, args.cycle, 0);
    

    /* kifmm[0,1] part */
    std::vector<real> root_p_check = kernel1.m2m_ewald_root0(args.P, cs[0], cs, args.cycle);
    //kernel1.direct(bs, direct_bodies2, 1, args.cycle, 0);

    /* merge */
    for(int i = 0; i < x_equiv.size(); i++)
    {
        //cs[0].p_check[i] += direct_bodies[i].p - direct_bodies2[i].p;
        cs[0].p_check[i] += direct_bodies[i].p - root_p_check[i];
    }
}

void rtfmm::EpkiFMM::calculate_root_check3()
{
    printf("calculate_root_check3\n");
    /* direct[2,args.images] part */
    LaplaceKernel kernel1;
    Bodies3 direct_bodies;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(args.P, cs[0].r * 1.05, cs[0].x);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 b;
        b.idx = i;
        b.x = x_equiv[i];
        b.q = 0;
        b.p = 0;
        b.f = vec3r(0,0,0);
        direct_bodies.push_back(b);
    }
    kernel1.direct(bs, direct_bodies, 3, args.cycle, 2);
    
    /* merge */
    for(int i = 0; i < x_equiv.size(); i++)
    {
        //cs[0].p_check[i] += direct_bodies[i].p - root_p_check[i];
        cs[0].p_check[i] += direct_bodies[i].p;
    }
}

void rtfmm::EpkiFMM::calculate_root_check4()
{
    LaplaceKernel kernel1;
    printf("calculate_root_check4\n");
    /* ewald[2,inf) part */
    Bodies3 ewald_bodies;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(args.P, cs[0].r * 1.05, cs[0].x);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 body;
        body.idx = i;
        body.p = 0;
        body.f = vec3r(0,0,0);
        body.x = x_equiv[i];
        body.q = cs[0].q_equiv[i];
        ewald_bodies.push_back(body);
    }
    Bodies3 direct_bodies2 = ewald_bodies;
    Argument args2 = args;
    args2.r *= 3;
    EwaldSolver ewald_solver(ewald_bodies, args2);
    ewald_bodies = ewald_solver.solve();
    dipole_correction(ewald_bodies, args2.cycle, 1);
    //kernel1.direct(bs, ewald_bodies, 3, args.cycle, 0);

    /* kifmm[0,1] part */
    LaplaceKernel kernel2;
    kernel2.direct(bs, direct_bodies2, 1, args2.cycle, 0);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        cs[0].p_check[i] += ewald_bodies[i].p - direct_bodies2[i].p;
    }
}

void rtfmm::EpkiFMM::calculate_root_check5()
{
    LaplaceKernel kernel1;
    printf("calculate_root_check5\n");
    /* ewald[2,inf) part */
    Bodies3 ewald_bodies;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(args.P, cs[0].r * 1.05, cs[0].x);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 body;
        body.idx = i;
        body.p = 0;
        body.f = vec3r(0,0,0);
        body.x = x_equiv[i];
        body.q = cs[0].q_equiv[i];
        ewald_bodies.push_back(body);
    }
    Bodies3 direct_bodies2 = ewald_bodies;
    Argument args2 = args;
    args2.r *= 1;
    EwaldSolver ewald_solver(ewald_bodies, bs, args2); 
    ewald_bodies = ewald_solver.solve();
    //dipole_correction(ewald_bodies, args.cycle, 1);

    {
        int num_body = bs.size();
        real coef = 4 * M_PI / (3 * std::pow(args2.cycle, 3));
        vec3r dipole(0,0,0);
        for (int i = 0; i < num_body; i++) 
        {
            dipole += bs[i].x * bs[i].q;
        }
        real dnorm = dipole.norm();
        for (int i = 0; i < num_body; i++) 
        {
            bs[i].p += 1 * coef * dnorm / num_body / bs[i].q; 
            bs[i].f += 1 * coef * dipole;
        }
    }

    /* kifmm[0,1] part */
    LaplaceKernel kernel2;
    kernel2.direct(bs, direct_bodies2, 1, args.cycle, 0);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        cs[0].p_check[i] += ewald_bodies[i].p - direct_bodies2[i].p;
    }
}


void rtfmm::EpkiFMM::calculate_root_check6()
{
    printf("calculate_root_check6()\n");

    // part 1
    LaplaceKernel kernel1;
    Bodies3 ewald_bodies1;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(args.P, cs[0].r * 1.05, cs[0].x);
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 body;
        body.idx = i;
        body.p = 0;
        body.f = vec3r(0,0,0);
        body.x = x_equiv[i];
        body.q = cs[0].q_equiv[i];
        ewald_bodies1.push_back(body);
    }
    Bodies3 direct_bodies1 = ewald_bodies1;
    Bodies3 direct_bodies2 = ewald_bodies1;

    EwaldSolver ewald_solver1(ewald_bodies1,ewald_bodies1, args);
    ewald_bodies1 = ewald_solver1.solve();
    kernel1.direct(direct_bodies2, direct_bodies1, 1, args.cycle, 0);

    // part 2
    LaplaceKernel kernel2;
    Bodies3 ewald_bodies2;
    for(int i = 0; i < x_equiv.size(); i++)
    {
        Body3 body;
        body.idx = i;
        body.p = 0;
        body.f = vec3r(0,0,0);
        body.x = vec3r(0.,0.,0.);
        body.q = cs[0].q_equiv[i];
        ewald_bodies2.push_back(body);
    }
    Bodies3 direct_source2 = ewald_bodies2;
    EwaldSolver ewald_solver2(ewald_bodies1, ewald_bodies2, args);
    ewald_bodies2 = ewald_solver2.solve();
    kernel2.direct(direct_source2, direct_bodies2, 1, args.cycle, 0);


    // 
    for(int i = 0; i < x_equiv.size(); i++)
    {
        cs[0].p_check[i] = (ewald_bodies1[i].p - direct_bodies1[i].p) 
        - (ewald_bodies2[i].p - direct_bodies2[i].p)
        ;
    }
}