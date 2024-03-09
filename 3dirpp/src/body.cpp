#include "body.h"

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