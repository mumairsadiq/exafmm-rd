#pragma once
#include <cstdlib>
#include <memory>
#include "type.h"

namespace rtfmm
{

const int MEM_ALIGN = 64;
const int CACHE_SIZE = 512;

template <typename T, size_t NALIGN>
struct AlignedAllocator : public std::allocator<T> {
  template <typename U>
  struct rebind {
    typedef AlignedAllocator<U, NALIGN> other;
  };

  T * allocate(size_t n) {
    void *ptr = nullptr;
    int rc = posix_memalign(&ptr, NALIGN, n * sizeof(T));
    if (rc != 0) return nullptr;
    if (ptr == nullptr) throw std::bad_alloc();
    return reinterpret_cast<T*>(ptr);
  }

  void deallocate(T * p, size_t) {
    return free(p);
  }
};

typedef AlignedAllocator<real, MEM_ALIGN> AlignAllocator;
typedef std::vector<real, AlignAllocator> AlignedVec;

}

