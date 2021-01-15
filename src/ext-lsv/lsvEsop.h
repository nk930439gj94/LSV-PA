#ifndef LSV_ESOP_H
#define LSV_ESOP_H

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "sat/cnf/cnf.h"
#include "bdd/cudd/cudd.h"
#include "bdd/cudd/cuddInt.h"

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
public:
  CofactorNode(Abc_Ntk_t* pNtk) {
    _pNtk = pNtk;
    _l = NULL; _r = NULL;
  }
  union {
    Abc_Ntk_t* _pNtk;
    Abc_Obj_t* _dvar;
  };
  CofactorNode* _l, * _r; // positive cofactor, negative cofactor
};

class CofactorTree
{
public:
  CofactorTree(Abc_Ntk_t* pNtkCone, Abc_Ntk_t* pNtk_origin);
  ~CofactorTree();
  int CofactorTree_rec(CofactorNode* n, Abc_Ntk_t* pNtk_origin, bool root = false);
  void CofactorTree_Delete_rec(CofactorNode* n, bool root = false);
  CofactorNode* _root;
  int level;
};

class TDDNode
{
public:
  TDDNode(DdNode* n) {
    _n = n;
    _l = _r = _x = 0;
  }
  DdNode* _n;
  TDDNode* _l, * _r, * _x; // positive cofactor, negative cofactor, boolean difference
};

class TDD
{
public:
  TDD(DdNode* n, DdManager* dd);
  int TDD_rec(TDDNode* t, DdManager* dd);
  TDDNode* _root;
};

#endif