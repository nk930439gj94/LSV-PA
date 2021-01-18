#ifndef LSV_CUBE_H
#define LSV_CUBE_H

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "bdd/cudd/cudd.h"
#include "bdd/cudd/cuddInt.h"
#include <string>

// 3-valued cube
// _v: value
// _c: care
struct Cube3
{
    Vec_Bit_t* _v;
    Vec_Bit_t* _c;
};

static inline Cube3 * Cube3Start( int nSize )
{
    Cube3 *c = ABC_ALLOC(Cube3, 1);
    c->_v = Vec_BitStart(nSize);
    c->_c = Vec_BitStart(nSize);
    return c;
}

static inline Cube3 * Cube3Dup( Cube3 * c )
{
    Cube3 * p = ABC_ALLOC(Cube3, 1);
    p->_v = Vec_BitDup(c->_v);
    p->_c = Vec_BitDup(c->_c);
    return p;
}

static inline void Cube3Free( Cube3 * c )
{
    Vec_BitFree(c->_v);
    Vec_BitFree(c->_c);
    ABC_FREE(c);
}

static inline int Cube3Size( Cube3 * c )
{
    return Vec_BitSize(c->_v);
}

static inline int Cube3Entry( Cube3 * c, int i )
{
    if(!Vec_BitEntry(c->_c, i)) return 2;
    return Vec_BitEntry(c->_v, i);
}

static inline void Cube3WriteEntry( Cube3 * c, int i, int Entry )
{
    assert(Entry>=0 && Entry<=2);
    if(Entry == 2) {
        Vec_BitWriteEntry(c->_c, i, 0);
    }
    else{
        Vec_BitWriteEntry(c->_v, i, Entry);
        Vec_BitWriteEntry(c->_c, i, 1);
    }
}

static inline std::string Cube3ToString( Cube3 * c )
{
    std::string str = "";
    int Entry, i, x;
    Vec_BitForEachEntry( c->_c, Entry, i ){
        x = Cube3Entry(c, i);
        if(x != 2){
            if(!str.empty()) str += ' ';
            str += std::to_string(i);
            if(x == 0) str += '\'';
        }
    }
    return str;
}

static inline Vec_Bit_t * Vec_BitXor( Vec_Bit_t * p0, Vec_Bit_t * p1 )
{
    assert(Vec_BitSize(p0) == Vec_BitSize(p1));
}

static inline Vec_Bit_t * Vec_BitOr( Vec_Bit_t * p0, Vec_Bit_t * p1 )
{
    assert(Vec_BitSize(p0) == Vec_BitSize(p1));
}

static inline int Cube3Distance( Cube3 * c0, Cube3 * c1 )
{
    assert(Cube3Size(c0) == Cube3Size(c1));
    Vec_Bit_t * v_xor = Vec_BitXor(c0->_v, c1->_v);
    Vec_Bit_t * c_xor = Vec_BitXor(c0->_c, c1->_c);
    Vec_Bit_t * diff = Vec_BitOr(v_xor, c_xor);
    int count = Vec_BitCount(diff);
    Vec_BitFree(v_xor);
    Vec_BitFree(c_xor);
    Vec_BitFree(diff);
    return count;
}

#endif