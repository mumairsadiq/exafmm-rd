#include "body.h"
#include "argument.h"

rtfmm::Bodies rtfmm::generate_random_bodies(int num, rtfmm::real r, vec3r offset, int seed, int zero_netcharge)
{
    Bodies bodies;
    double q_avg = 0;
    srand48(seed);
	for(int i = 0; i < num; i++)
	{
        Body body;
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
    if(zero_netcharge)
    {
        for(int i = 0; i < num; i++)
        {
            bodies[i].q -= q_avg;
        }
    }

    return bodies;
}

rtfmm::BodyCompareResult rtfmm::compare(const Bodies& bs1, const Bodies& bs2, std::string name1, std::string name2, int num_compare)
{
    assert_exit(bs1.size() == bs2.size(), "inconsistent size in comparison");

    BodyCompareResult res;
    res.name1 = name1;
    res.name2 = name2;

    int num = num_compare == -1 ? bs1.size() : std::min(num_compare, (int)bs1.size());
    res.num_compared = num;

    real pdif = 0, pnrm = 0;
    real fdif = 0, fnrm = 0;
    real esum1 = 0, esum2 = 0;
    for(int i = 0; i < num; i++)
    {
        Body b1 = bs1[i];
        Body b2 = bs2[i];
        esum1 += b1.p * b1.q;                   
        esum2 += b2.p * b2.q;
        pdif += std::pow(b1.p - b2.p, 2);
        pnrm += std::pow(b2.p, 2);
        vec3r diff = b1.f - b2.f;
        fdif += diff.norm();
        fnrm += b2.f.norm();
        int flag = diff.r() > 1 ? 1 : 0;
        /*if(flag)
        printf("[%d]  %d  %.4f(%.4f,%.4f,%.4f)      %.8f (%.8f,%.8f,%.8f)   %.8f (%.8f,%.8f,%.8f)   %.8f (%.8f,%.8f,%.8f)\n", 
            i, flag,
            b1.q, b1.x[0], b1.x[1], b1.x[2],
            b1.p, b1.f[0], b1.f[1], b1.f[2],
            b2.p, b2.f[0], b2.f[1], b2.f[2],
            std::abs(b1.p - b2.p), diff[0], diff[1], diff[2]);*/
    }
    //printf("pnrm = %.8f\n", pnrm);
    res.rmsp = std::sqrt(pdif / num);
    res.rmsf = std::sqrt(fdif / num);
    res.l2p = std::sqrt(pdif / pnrm);
    res.l2f = std::sqrt(fdif / fnrm);
    res.epot1 = esum1;
    res.epot2 = esum2;
    res.l2e = std::sqrt(std::pow(esum1 - esum2, 2) / esum2 / esum2);

    return res;
}

void rtfmm::BodyCompareResult::show()
{
    int bar_num = 72 - name1.size() - name2.size() - 4;
    int bar_num_left = bar_num / 2;
    for(int i = 0; i < bar_num_left; i++) printf("-");
    std::cout<<name1<<" vs "<<name2;
    for(int i = 0; i < bar_num - bar_num_left; i++) printf("-");
    printf("[%d]\n", num_compared);
    printf("%-8s : %8.5e   %-8s : %8.5e   %-8s : %8.5e\n", "L2  (p)", l2p , "L2  (f)", l2f, "L2  (e)", l2e); 
    printf("%-8s : %8.5e   %-8s : %8.5e\n", "Rms (p)", rmsp, "Rms (f)", rmsf);
    printf("p-energy1 : %8.12e\n", epot1);
    printf("p-energy2 : %8.12e\n", epot2);
    printf("\n");
}

void rtfmm::dipole_correction(Bodies& bs, real cycle)
{
    if(verbose) std::cout<<"dipole correction"<<std::endl;
    int num_body = bs.size();
    real coef = 4 * M_PI / (3 * cycle * cycle * cycle);
    vec3r dipole(0,0,0);
    for (int i = 0; i < num_body; i++) 
    {
        dipole += bs[i].x * bs[i].q;
    }
    real dnorm = dipole.norm();
    for (int i = 0; i < num_body; i++) 
    { 
        bs[i].p -= coef * dnorm / num_body / bs[i].q; 
        bs[i].f -= coef * dipole;
    }
}

std::vector<rtfmm::vec3r> rtfmm::get_bodies_x(Bodies& bs, vec2i range, vec3r offset)
{
    std::vector<rtfmm::vec3r> res(range[1]);
    for(int i = 0; i < range[1]; i++)
    {
        res[i] = bs[range[0] + i].x + offset;
    }
    return res;
}

rtfmm::Matrix rtfmm::get_bodies_q(rtfmm::Bodies& bs, vec2i range)
{
    Matrix res(range[1], 1);
    for(int i = 0; i < range[1]; i++)
    {
        res[i] = bs[range[0] + i].q;
    }
    return res;
}

rtfmm::Matriv rtfmm::get_force_naive(
    std::vector<vec3r>& x_src, 
    std::vector<vec3r>& x_tar, 
    Matrix& q_src
)
{
    Matriv res(x_tar.size(), 1);
    for(int j = 0; j < x_tar.size(); j++)
    {   
        vec3r force(0,0,0);
        for(int i = 0; i < x_src.size(); i++)
        {
            vec3r dx = x_tar[j] - x_src[i];
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            force += q_src[i] * invr * invr * invr * (-dx);
        }
        res[j] = force;
    }
    return res;
}