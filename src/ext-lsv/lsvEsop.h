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
  Abc_Ntk_t * Abc_NtkFromGlobalBdds( Abc_Ntk_t * pNtk, int fReverse );
  int Abc_NtkRewrite( Abc_Ntk_t * pNtk, int fUpdateLevel, int fUseZeros, int fVerbose, int fVeryVerbose, int fPlaceEnable );
  void Abc_NodeShowBdd(Abc_Obj_t * pNode, int fCompl);
  void Abc_NtkShow(Abc_Ntk_t * pNtk0, int fGateNames, int fSeq, int fUseReverse);
}

#define SupportThreshold 20
#define Cudd_Index(node) ((Cudd_Regular(node))->index)

class CofactorTree;
class CofactorNode;
class TDD;
class TDDNode;

int Abc_NtkSupportNum(Abc_Ntk_t* pNtk);
void dfs_rec(Abc_Obj_t* pNode);

static Abc_Ntk_t* Collapse_reservePi(Abc_Ntk_t * pNtk, int fReorder, Vec_Int_t*& perm);
static Abc_Ntk_t* Cofactor(Abc_Ntk_t* pNtk, bool fPos, int iVar);

void esopStats(Vec_Ptr_t* cubes);
void esopSimplify(Vec_Ptr_t* cubes);
void esopPrint(Vec_Ptr_t* cubes, Vec_Ptr_t* PiNames = 0);
void esopFree(Vec_Ptr_t* cubes);

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
  CofactorNode* _l;
  union {
    CofactorNode* _r;
    Vec_Int_t* _PiMap;
  };
};

class CofactorTree
{
friend class TDD;
public:
  CofactorTree(Abc_Ntk_t* pNtkCone);
  ~CofactorTree();
  Vec_Ptr_t* toEsop();
  static void setGlobalNtk(Abc_Ntk_t* pNtk_global);
private:
  CofactorNode* _root;
  static Abc_Ntk_t* _pNtk_global;
  void CofactorTree_rec(CofactorNode* n, bool root = false);
  void CofactorTree_Delete_rec(CofactorNode* n, bool root = false);
  void toEsop_rec(CofactorNode* n, Cube3* factor, Vec_Ptr_t* cubes);
  void setPiMap(Abc_Ntk_t* pNtk, Vec_Int_t* PiMap);
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
  TDD(DdNode* n, Abc_Ntk_t* pNtk, Vec_Int_t* PiMap);
  ~TDD();
  void toEsop(Cube3* cube, Vec_Ptr_t* cubes);
  int nCubes() {return _nCubes;}
private:
  TDDNode* _root;
  int _nCubes;
  Abc_Ntk_t* _pNtk;
  Vec_Int_t* _PiMap;
  int TDD_rec(TDDNode* t);
  void TDD_Delete_rec(TDDNode* tn);
  void toEsop_rec(TDDNode* tn, Cube3* cube, Vec_Ptr_t* cubes);
};

#endif
