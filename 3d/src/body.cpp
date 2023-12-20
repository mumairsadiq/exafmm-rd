#include "body.h"
#include <random>

void rtfmm::print_body(const Body3& b)
{
    printf("[%d],(%.4f,%.4f,%.4f),%.4f,%.8f,(%.8f,%.8f,%.8f)\n", 
        b.idx,
        b.x[0],b.x[1],b.x[2],
        b.q,
        b.p,
        b.f[0],b.f[1],b.f[2]
    );
}

void rtfmm::print_bodies(const rtfmm::Bodies3& bs, int num, int offset)
{
    printf("bodies : \n");
    if(num == -1) num = bs.size();
    int s = std::min(num, (int)bs.size());
    for(int i = offset; i < offset + s; i++)
    {
        print_body(bs[i]);
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