#ifndef LSV_CUBE_H
#define LSV_CUBE_H

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "bdd/cudd/cudd.h"
#include "bdd/cudd/cuddInt.h"

// 3-valued cube
// _v: value
// _d: don't care
struct Cube3
{
    Vec_Bit_t* _v;
    Vec_Bit_t* _d;
};

static inline Cube3 * Cube3Start( int nSize )
{
    Cube3 *c = ABC_ALLOC(Cube3, 1);
    c->_v = Vec_BitStart(nSize);
    c->_d = Vec_BitStartFull(nSize);
    return c;
}

static inline void Cube3Free( Cube3 * c )
{
    Vec_BitFree(c->_v);
    Vec_BitFree(c->_d);
    ABC_FREE(c);
}

static inline int Cube3Size( Cube3 * c )
{
    return Vec_BitSize(c->_v);
}

static inline int Cube3Entry( Cube3 * c, int i )
{
    if(Vec_BitEntry(c->_d, i)) return 2;
    return Vec_BitEntry(c->_v, i);
}

static inline void Cube3WriteEntry( Cube3 * c, int i, int Entry )
{
    assert(Entry>=0 && Entry<=2);
    if(Entry == 2) {
        Vec_BitWriteEntry(c->_v, i, 0);
        Vec_BitWriteEntry(c->_d, i, 1);
    }
    else{
        Vec_BitWriteEntry(c->_v, i, Entry);
        Vec_BitWriteEntry(c->_d, i, 0);
    }
}

#endif