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
    c->_v->nSize = nSize;
    c->_c = Vec_BitStart(nSize);
    return c;
}

static inline Cube3 * Cube3Dup( Cube3 * c )
{
    Cube3 * p = ABC_ALLOC(Cube3, 1);
    int nSize = c->_v->nSize;
    c->_v->nSize = c->_c->nSize;
    p->_v = Vec_BitDup(c->_v);
    c->_v->nSize = nSize;
    p->_c = Vec_BitDup(c->_c);
    p->_v->nSize = nSize;
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

static inline int Cube3CountLiteral( Cube3 * c )
{
    return Vec_BitCount(c->_c);
}

static inline void Cube3MergeDist1( Cube3 * c0, Cube3 * c1 )
{
    // merge cubes with distance equal to 1
    // merged cube replaces c0
    int size = Cube3Size(c0);
    assert(Cube3Size(c1) == size);
    int e0, e1;
    for(int i = 0; i < size; ++i) {
        e0 = Cube3Entry(c0, i);
        e1 = Cube3Entry(c1, i);
        if(e0 != e1) {
            if(e0 != 0 && e1 != 0) Cube3WriteEntry(c0, i, 0);
            else if(e0 != 1 && e1 != 1) Cube3WriteEntry(c0, i, 1);
            else Cube3WriteEntry(c0, i, 2);
            break;
        }
    }
}

static inline std::string Cube3ToString( Cube3 * c, Vec_Ptr_t* PiNames = 0 )
{
    std::string str = "";
    int x;
    for(int i = 0; i < Cube3Size(c); ++i){
        x = Cube3Entry(c, i);
        if(x != 2){
            if(!str.empty()) str += ' ';
            if(PiNames) str += std::string((char*)Vec_PtrEntry(PiNames, i));
            else str += std::to_string(i);
            if(x == 0) str += '\'';
        }
    }
    return str.empty() ? "\"T\"" : str;
}

static inline Vec_Bit_t * Vec_BitXor( Vec_Bit_t * p0, Vec_Bit_t * p1, int realSize )
{
    assert(Vec_BitSize(p0) == Vec_BitSize(p1));
    int nWords = (p0->nSize >> 5) + ((p0->nSize & 31) > 0);
    Vec_Bit_t * v_xor = Vec_BitStart(realSize);
    int i;
    for ( i = 0; i < nWords; i++ )
        v_xor->pArray[i] = p0->pArray[i] ^ p1->pArray[i];
    if ( realSize & 31 )
        v_xor->pArray[nWords-1] &= ~(~0 << (realSize & 31));
    return v_xor;
}

static inline Vec_Bit_t * Vec_BitOr( Vec_Bit_t * p0, Vec_Bit_t * p1, int realSize )
{
    assert(Vec_BitSize(p0) == Vec_BitSize(p1));
    int nWords = (p0->nSize >> 5) + ((p0->nSize & 31) > 0);
    Vec_Bit_t * v_or = Vec_BitStart(realSize);
    int i;
    for ( i = 0; i < nWords; i++ )
        v_or->pArray[i] = p0->pArray[i] | p1->pArray[i];
    if ( realSize & 31 )
        v_or->pArray[nWords-1] &= ~(~0 << (realSize & 31));
    return v_or;
}

static inline int Cube3Distance( Cube3 * c0, Cube3 * c1 )
{
    int size = Cube3Size(c0);
    assert(Cube3Size(c1) == size);
    Vec_Bit_t * v_xor = Vec_BitXor(c0->_v, c1->_v, size);
    Vec_Bit_t * c_xor = Vec_BitXor(c0->_c, c1->_c, size);
    Vec_Bit_t * diff = Vec_BitOr(v_xor, c_xor, size);
    int count = Vec_BitCount(diff);
    Vec_BitFree(v_xor);
    Vec_BitFree(c_xor);
    Vec_BitFree(diff);
    return count;
}

#endif
