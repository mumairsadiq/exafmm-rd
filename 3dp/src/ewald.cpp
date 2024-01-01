#include "ewald.h"

rtfmm::EwaldSolver::EwaldSolver(const Bodies3& bs_, const Argument& args_) 
: bs(bs_), args(args_)
{
    ksize = args.ewald_ksize;
    printf("ksize = %d\n", ksize);
    scale = 2 * M_PI / args.cycle;
    alpha = ksize / args.cycle;
    cutoff = args.cycle / 2;
    for (int k = -ksize; k <= ksize; k++)
    {                               
        for (int m = -ksize; m <= ksize; m++)
        {                          
            for (int p = -ksize; p <= ksize; p++)
            {                        
                if(k == 0 && m == 0 && p == 0) continue;
                Wave wave;
                wave.K = rtfmm::vec3r(k,m,p);
                wave.val = std::complex<double>(0,0);
                waves.push_back(wave);
            }
        }
    }
}

rtfmm::Bodies3 rtfmm::EwaldSolver::solve()
{
    rtfmm::Tree tree;
    tree.build(bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    //tree.build(bs, args.x, args.r, 3, Tree::TreeType::uniform);
    cells = tree.get_cells();
    TIME_BEGIN(real_part);
    real_part(0, 0);
    if(args.timing) {TIME_END(real_part);}
    TIME_BEGIN(fourier_part);
    fourier_part();
    if(args.timing) {TIME_END(fourier_part);}
    TIME_BEGIN(self_correction);
    self_correction();
    if(args.timing) {TIME_END(self_correction);}

    return sort_bodies_by_idx(bs);
}

void rtfmm::EwaldSolver::real_part(int this_cell_idx, int that_cell_idx)
{
    rtfmm::Cell3& this_cell = cells[this_cell_idx];
    if (this_cell.crange.number == 0) 
    {
        child(this_cell_idx, that_cell_idx);
    }
    else
    {
        for(int i = 0; i < this_cell.crange.number; i++)
        {
            real_part(this_cell.crange.offset + i, that_cell_idx);
        }
    }
}


void rtfmm::EwaldSolver::fourier_part()
{
    DFT(waves,bs);
    rtfmm::real coef = 4 * M_PI / std::pow(args.cycle,3);
    for (int w = 0; w < waves.size(); w++) 
    {
        rtfmm::real k2 = waves[w].K.norm();
        waves[w].val *= coef * std::exp(-k2 / (4 * alpha * alpha)) / k2;
    }
    IDFT(waves,bs);
}


void rtfmm::EwaldSolver::self_correction()
{
    for (int i = 0; i < bs.size(); i++) 
    {
        bs[i].p -= M_2_SQRTPI * bs[i].q * alpha;
    }
}


void rtfmm::EwaldSolver::child(int this_cell_idx, int that_cell_idx)
{
    rtfmm::Cell3& this_cell = cells[this_cell_idx];
    rtfmm::Cell3& that_cell = cells[that_cell_idx];
    rtfmm::vec3r offset(0,0,0);
    rtfmm::vec3r dx = this_cell.x - that_cell.x;
    for (int d = 0; d < 3; d++)
    {
        if(dx[d] < -args.cycle / 2)
        {
            offset[d]-=args.cycle;
        }
        else if(dx[d] >  args.cycle / 2)
        {
            offset[d]+=args.cycle;
        }
    }
    if ((dx - offset).r() - this_cell.r - that_cell.r < sqrtf(3) * cutoff)
    {
        if(that_cell.crange.number == 0)
        {
            real_p2p(this_cell_idx, that_cell_idx, offset);
        }
        else
        {
            for (int j = 0; j < that_cell.crange.number; j++)
            {
                child(this_cell_idx, that_cell.crange.offset + j);
            }
        }
    }
}


void rtfmm::EwaldSolver::real_p2p(int this_cell_idx, int that_cell_idx, rtfmm::vec3r offset)
{
    rtfmm::Cell3& this_cell = cells[this_cell_idx];
    rtfmm::Cell3& that_cell = cells[that_cell_idx];
    for (int i = 0; i < this_cell.brange.number; i++)
    {
        rtfmm::Body3& bi = bs[this_cell.brange.offset + i];
        for (int j = 0; j < that_cell.brange.number; j++)
        {  
            rtfmm::Body3& bj = bs[that_cell.brange.offset + j];
            rtfmm::vec3r dx = bi.x - bj.x - offset;
            rtfmm::real r = dx.r();
            if (r > 0 && r < cutoff)
            {
                bi.p += bj.q * std::erfc(alpha * r) / r;
                bi.f += -bj.q * dx / std::pow(r,3) * (std::erfc(alpha * r) + (2 * alpha * r * std::pow(M_E, -alpha * alpha * r * r)) / std::sqrt(M_PI));
            }
        }
    }
}


void rtfmm::EwaldSolver::DFT(std::vector<Wave>& ws, std::vector<rtfmm::Body3>& bs)
{
    for (int w = 0; w < ws.size(); w++) 
    {                     
        for (int i = 0; i < bs.size(); i++)  
        {                 
            rtfmm::real ph = (ws[w].K * bs[i].x).sum() * scale;                                
            ws[w].val += bs[i].q * std::pow(rtfmm::real(M_E), complexr(0, -ph));
        }                                                       
    } 
}


void rtfmm::EwaldSolver::IDFT(std::vector<Wave>& ws, std::vector<rtfmm::Body3>& bs)
{
    for (int i = 0; i < bs.size(); i++)
    {
        for (int w = 0; w < ws.size(); w++)
        {
            rtfmm::real ph = (ws[w].K * bs[i].x).sum() * scale;
            bs[i].p += std::real(ws[w].val * std::pow(rtfmm::real(M_E), complexr(0, ph)));
            bs[i].f += std::real(ws[w].val * std::pow(rtfmm::real(M_E), complexr(0, ph)) * complexr(0, 1) * scale) * ws[w].K; // actually this is dp/dx, while real f is -dp/dx*q
        }              
    }   
}
