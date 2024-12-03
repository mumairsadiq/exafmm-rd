
#ifndef GMX_FMM_TYPES_H
#define GMX_FMM_TYPES_H

#include "type.h"
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gmx
{
using real = rtfmm::real;
using RVec = rtfmm::vec3r;
using BVec = rtfmm::vec<3, bool>;
namespace fmm
{

using FPIndices = std::vector<int>;
struct FBody
{
    int idx;
    RVec x;
    real q;
};

struct OffsetAndNumber
{
    int offset;
    int number;

    OffsetAndNumber() {}
    OffsetAndNumber(int offset_, int number_) : offset(offset_), number(number_)
    {
    }

    friend std::ostream &operator<<(std::ostream &os, const OffsetAndNumber &on)
    {
        os << "[" << on.offset << "," << on.number << "]";
        return os;
    }
};

using Range = OffsetAndNumber;

struct FBodies : public std::vector<FBody>
{

    FBodies() = default;
    // Constructor that calls initialize from the initializer list
    FBodies(const std::vector<RVec> &coordinates,
            const std::vector<real> &charges)
    {
        initialize(coordinates, charges);
    }

  private:
    // Helper function to initialize FBodies from coordinates and charges
    void initialize(const std::vector<RVec> &coordinates,
                    const std::vector<real> &charges)
    {
        this->reserve(coordinates.size()); // Reserve memory for efficiency

        for (size_t i = 0; i < coordinates.size(); ++i)
        {
            FBody body;
            body.idx = i;
            body.x = coordinates[i];
            body.q = charges[i];
            this->push_back(body);
        }
    }
};

struct FMMCell
{
    size_t index;      // Unique cell index
    int octant;        // Octant relative to the parent cell (0-7 for relative
                       // position, 13 for central/hitorikko)
    int depth;         // Depth of the cell
    real radius;       // Radius of the cell
    RVec center;       // Center coordinates of the cell
    Range crange;      // Children information of the cell
    real radiusParent; // Radius of parent of cell
    RVec centerParent; // Center coordinates of the parent of the cell
    FPIndices bodiesIndices; // Indices of bodies in the cell
    FPIndices bodiesIndicesWithReg;
    std::vector<real> weights;

    // Constructors
    FMMCell()
    {
        this->index = 0;
        this->crange = Range(0, 0);
        this->octant = 0;
        this->depth = 0;
    }

    bool isLeaf() const { return this->crange.number == 0; }

    // Public member functions
    bool operator<(const FMMCell &other) const
    {
        if (center[0] != other.center[0])
            return center[0] < other.center[0];
        if (center[1] != other.center[1])
            return center[1] < other.center[1];
        if (center[2] != other.center[2])
            return center[2] < other.center[2];
        return index < other.index; // If centers are equal, compare by index
    }
};

using FMMCells = std::vector<FMMCell>;

#ifndef _PRECISION_FMM_Adj
#define _PRECISION_FMM_Adj 1e-3
#endif

inline long coord_to_long(const RVec x, const RVec x_, const real r)
{
    // Define the origin of the box (corner) based on center and radius
    real ox = x_[0] - r;
    real oy = x_[1] - r;
    real oz = x_[2] - r;

    // Quantize coordinates directly
    long ix = static_cast<long>((x[0] - ox) / _PRECISION_FMM_Adj);
    long iy = static_cast<long>((x[1] - oy) / _PRECISION_FMM_Adj);
    long iz = static_cast<long>((x[2] - oz) / _PRECISION_FMM_Adj);

    return (ix << 42) | (iy << 21) | iz;
}

} // namespace fmm
} // namespace gmx
#endif