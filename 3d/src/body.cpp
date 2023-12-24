#include "body.h"
#include <random>
#include <algorithm>

void rtfmm::print_body(const Body3& b)
{
    printf("[%d],(%.4f,%.4f,%.4f),%.4f,%.12f,(%.8f,%.8f,%.8f)\n", 
        b.idx,
        b.x[0],b.x[1],b.x[2],
        b.q,
        b.p,
        b.f[0],b.f[1],b.f[2]
    );
}

void rtfmm::print_bodies(const rtfmm::Bodies3& bs, int num, int offset, std::string name)
{
    std::cout<<"\n"<<name<<":"<<std::endl;
    if(num == -1) num = bs.size();
    int s = std::min(num, (int)bs.size());
    for(int i = offset; i < offset + s; i++)
    {
        print_body(bs[i]);
    }
}

std::vector<rtfmm::vec3r> rtfmm::get_bodies_x(rtfmm::Bodies3& bs, Range range)
{
    std::vector<rtfmm::vec3r> res(range.number);
    for(int i = 0; i < range.number; i++)
    {
        res[i] = bs[range.offset + i].x;
    }
    return res;
}

rtfmm::Matrix rtfmm::get_bodies_q(rtfmm::Bodies3& bs, Range range)
{
    Matrix res(range.number, 1);
    for(int i = 0; i < range.number; i++)
    {
        res[i] = bs[range.offset + i].q;
    }
    return res;
}

void rtfmm::set_boides_p(Bodies3& bs, Matrix& ps, Range range)
{
    assert(ps.m * ps.n == range.number);
    int num = bs.size();
    for(int i = 0; i < num; i++)
    {
        bs[range.offset + i].p = ps[i];
    }
}

void rtfmm::set_boides_f(Bodies3& bs, Matriv& fs, Range range)
{
    assert(fs.m * fs.n == range.number);
    int num = bs.size();
    for(int i = 0; i < num; i++)
    {
        bs[range.offset + i].f = fs[i];
    }
}

rtfmm::Bodies3 rtfmm::generate_random_bodies(int num, rtfmm::real r, vec3r offset)
{
    Bodies3 bodies;
    double q_avg = 0;
    srand48(0);
	for(int i = 0; i < num; i++)
	{
        Body3 body;
        for(int d = 0; d < 3; d++)
        {
            body.x[d] = drand48() * r * 2 - r;
        }
        double q = drand48() - 0.5;
        q_avg += q;
        body.idx = i;
        body.q = q;
        body.p = 0;
        body.f = vec3r(0,0,0);
        body.x += offset;
        bodies.push_back(body);
	}
    q_avg /= num;
    for(int i = 0; i < num; i++)
	{
        bodies[i].q -= q_avg;
    }

    return bodies;
}

rtfmm::BodyCompareResult rtfmm::compare(const Bodies3& bs1, const Bodies3& bs2, std::string name1, std::string name2)
{
    assert(bs1.size() == bs2.size());

    BodyCompareResult res;
    res.name1 = name1;
    res.name2 = name2;

    int num = bs1.size();

    real psum1 = 0, psum2 = 0, rmsp = 0;
    real fdif = 0, fnrm = 0;
    for(int i = 0; i < num; i++)
    {
        Body3 b1 = bs1[i];
        Body3 b2 = bs2[i];
        psum1 += b1.p * b1.q;                   
        psum2 += b2.p * b2.q;
        rmsp += std::pow(b1.p - b2.p, 2);
        fdif += (b1.f - b2.f).norm();
        fnrm += b1.f.norm();
    }

    res.rmsp = std::sqrt(rmsp / num);
    res.rmsf = std::sqrt(fdif / num);
    res.l2p = sqrt(std::pow(psum1 - psum2, 2) / psum1 / psum1);  
    res.l2f = std::sqrt(fdif / fnrm);
    res.epot1 = psum1;
    res.epot2 = psum2;

    return res;
}

rtfmm::Bodies3 rtfmm::sort_bodies_by_idx(const Bodies3& bs)
{
    rtfmm::Bodies3 bodies = bs;
    std::sort(
        bodies.begin(), 
        bodies.end(),
        [](const Body3& a, const Body3& b)
        {
            return a.idx < b.idx; 
        }
    );

    return bodies;
}
