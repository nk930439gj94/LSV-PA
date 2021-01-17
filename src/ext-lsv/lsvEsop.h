#ifndef LSV_ESOP_H
#define LSV_ESOP_H

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "bdd/cudd/cudd.h"
#include "bdd/cudd/cuddInt.h"

#include "lsvCube.h"

extern "C"
{
  Abc_Ntk_t* Abc_NtkCollapse(Abc_Ntk_t* pNtk, int fBddSizeMax, int fDualRail, int fReorder, int fReverse, int fDumpOrder, int fVerbose);
  void Abc_NodeShowBdd(Abc_Obj_t * pNode, int fCompl);
}

#define AigNodeThreshold 20
#define Cudd_Index(node) ((Cudd_Regular(node))->index)

class CofactorTree;
class CofactorNode;
class TDD;
class TDDNode;

static int BDD_nCube(DdNode* n);
static Abc_Ntk_t* Cofactor(Abc_Ntk_t* pNtk, bool fPos, int iVar);


class CofactorNode
{
friend class CofactorTree;
public:
  CofactorNode(Abc_Ntk_t* pNtk) {
    _pNtk = pNtk;
    _l = NULL; _r = NULL;
  }
  bool isLeaf() {return !_l;}
private:
  union {
    Abc_Ntk_t* _pNtk;
    int _dvar;
  };
  CofactorNode* _l, * _r; // positive cofactor, negative cofactor
};

class CofactorTree
{
friend class TDD;
public:
  CofactorTree(Abc_Ntk_t* pNtkCone);
  ~CofactorTree();
  Vec_Ptr_t* toEsop();
  static void setGlobalNtk(Abc_Ntk_t* pNtk_global) {
    _pNtk_global = pNtk_global;
  }
  static void setGlobalPiReference() {
    Abc_Obj_t* pPi; int i;
    Abc_NtkForEachPi(_pNtk_global, pPi, i) pPi->iTemp = i;
  }
private:
  CofactorNode* _root;
  static Abc_Ntk_t* _pNtk_global;
  void CofactorTree_rec(CofactorNode* n, bool root = false);
  void CofactorTree_Delete_rec(CofactorNode* n, bool root = false);
  void toEsop_rec(CofactorNode* n, Cube3* factor, Vec_Ptr_t* cubes);
};

class TDDNode
{
friend class TDD;
public:
  TDDNode(DdNode* n) {
    _n = n;
    _l = _r = _x = 0;
  }
private:
  DdNode* _n;
  TDDNode* _l, * _r, * _x; // positive cofactor, negative cofactor, boolean difference
};

class TDD
{
public:
  TDD(DdNode* n, Abc_Ntk_t* pNtk);
  ~TDD();
  void toEsop(Cube3* cube, Vec_Ptr_t* cubes);
  int nCubes() {return _nCubes;}
private:
  TDDNode* _root;
  int _nCubes;
  Abc_Ntk_t* _pNtk;
  int TDD_rec(TDDNode* t);
  void TDD_Delete_rec(TDDNode* tn);
  void toEsop_rec(TDDNode* tn, Cube3* cube, Vec_Ptr_t* cubes);
};

#endif