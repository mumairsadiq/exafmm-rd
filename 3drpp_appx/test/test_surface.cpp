#include "body.h"
#include "surface.h"
#include "type.h"

int main(int argc, char *argv[])
{
    std::cout << "rtfmm_3d_test_surface" << std::endl;

    std::vector<rtfmm::vec3r> points = rtfmm::get_surface_points(16);
    std::cout << points.size() << std::endl;

    return 0;
}