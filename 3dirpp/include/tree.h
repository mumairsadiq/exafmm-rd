#pragma once
#include "type.h"
#include <vector>
#include <map>
#include "argument.h"
#include "body.h"

namespace rtfmm
{

class LogicCoord : public vec4i
{
public:
    LogicCoord(){}
    LogicCoord(int v1, int v2, int v3, int v4) : vec4i(v1, v2, v3, v4) {}
    LogicCoord(const vec4i &v) 
    {
        data[0] = v[0];
        data[1] = v[1];
        data[2] = v[2];
        data[3] = v[3];
    }
    vec3i idx3()
    {
        return vec3i(data[1], data[2], data[3]);
    }
    LogicCoord parent()
    {
        return LogicCoord(data[0] - 1, data[1] / 2, data[2] / 2, data[3] / 2);
    }
    LogicCoord child(const int& octant)
    {
        // TODO : 27 octants
        return LogicCoord(
            data[0] + 1, 
            data[1] * 2 + ((octant >> 0) & 1),  // note : this should match rtfmm::Tree::get_body_octant
            data[2] * 2 + ((octant >> 1) & 1), 
            data[3] * 2 + ((octant >> 2) & 1)
        );
    }
    int depth_wise_idx() const
    {
        int N = data[0] >= 0 ? std::pow(2, data[0]) : 3;
        return data[1] + data[2] * N + data[3] * N * N;
    }
    int hash(int images) const
    {
        int num_imgcell = std::min(data[0] + images, images) * 27;
        int num_nonimgcell = (std::pow(8, data[0]) - 1) / 7;
        int depth_wise_offset = depth_wise_idx();
        return num_imgcell + num_nonimgcell + depth_wise_offset;
    }
    int legal() const
    {
        if(data[0] >= 0)
        {
            int N = std::pow(2, data[0]);
            if(data[1] >= 0 && data[1] < N && data[2] >= 0 && data[2] < N && data[3] >= 0 && data[3] < N)
            {
                return 1;
            }
        }
        else
        {
            if(data[1] >= -1 && data[1] <= 1 && data[2] >= -1 && data[2] <= 1 && data[3] >= -1 && data[3] < 1)
            {
                return 1;
            }
        }
        return 0;
    }
    static int adjacent(const LogicCoord& xj, const LogicCoord& xi)
    {
        LogicCoord a, b; // a.depth >= b.depth
        if(xj[0] >= xi[0])
        {
            a = xj; 
            b = xi;
        }
        else
        {
            a = xi; 
            b = xj;
        }
        int N = std::pow(2, a[0] - b[0]);
        int r = 1 + N;
        vec3i xa = 2 * a.idx3() + 1;
        vec3i xb = (2 * b.idx3() + 1 ) * N; 
        vec3i dx = (xa - xb).abs();
        return (dx[0] <= r && dx[1] <= r && dx[2] <= r);
    }
};

struct Cell
{
    // logic position of cell, e.g., (2,1,3,0) stands for cell in (x=1,y=3,z=0) posiiton in depth=2 among 2^2^3=64 cells 
    LogicCoord xlogic;

    // physic location of cell
    vec3r xphys;

    // index of layer's memory storing the cell
    int xmem;

    // half of edge length    
    real r;

    // hash value of xlogic
    int key;

    // if is a leaf cell, i.e., has neither nonregbody nor regbody
    int is_leaf;

    // range of nonregbody
    vec2i brange;

    /**
     * @brief equivalent charge matrix index, -1 means unassigned
     * @note different Cells could have same qidx, 
     * this often means they are image-copies of a certain cell.
     * (image-copies : cells have same M but different positions)
    */
    int qidx;

    // check potential matrix index, -1 means unassigned, this often means its L is unnecessary
    int pidx;    

    friend std::ostream &operator<<(std::ostream & os, const Cell & cell) 
    {
        os << "<Cell> "
        <<"xlogic=" << cell.xlogic
        << ", xphys=" << cell.xphys
        << ", xmem=" << cell.xmem
        << ", r=" << cell.r
        << ", key=" << cell.key
        << ", is_leaf=" << cell.is_leaf
        << ", brange=" << cell.brange
        << ", (qidx,pidx)=(" << cell.qidx << "," << cell.pidx << ")"
        << ", (M,L)=("<<cell.M<<","<<cell.L<<")";
        return os;
    }

    Cell()
    {
        M = 0;
        L = 0;
    }

    int M;
    int L;
};

using Cells = std::vector<Cell>;

class Tree
{
public:
    Tree(Bodies& bs, const Argument& args);

    int find_cidx(const LogicCoord& xlogic)
    {
        if(!xlogic.legal())
            return -1;
        int cidx = -1;
        const int key = xlogic.hash(images);
        if(cidx_map.count(key))
            cidx = cidx_map[key];
        return cidx;
    }

    static vec3r logic2phys(const LogicCoord& xlogic, const real& r0, const vec3r& x0)
    {
        vec3r res;

        if(xlogic[0] >= 0)
        {
            int N = std::pow(2, xlogic[0]);
            real r = r0 / N;
            res = x0 + (vec3r(xlogic[1],xlogic[2],xlogic[3]) * 2 - N + 1) * r;
        }
        else
        {
            int N = std::pow(3, -xlogic[0] - 1);
            real r = r0 * N;
            res = x0 + vec3r(xlogic[1]-1,xlogic[2]-1,xlogic[3]-1) * 2 * r;
        }

        return res;
    }

    static int get_body_octant(const vec3r& bx, const vec3r& cx, const real& cr)
    {
        // note : this should match rtfmm::LogicCoord::child
        int idx = ((bx[0] > cx[0]) << 0)  //x
                + ((bx[1] > cx[1]) << 1)  //y
                + ((bx[2] > cx[2]) << 2); //z
        return idx;
    }

public:
    real r0;
    vec3r x0;
    int ncrit;
    Cells cs;
    int images;

    std::map<int, int> cidx_map; // hash(xlogic) -> xmem, -1 mean non-existent

    vec2i depth_range;

private:
    void build_octree(Bodies& bs);
    void build_base27tree(int images);
    void precompute_attributes();
};

}