#include "lsvEsop.h"

Abc_Ntk_t* CofactorTree::_pNtk_global = 0;

CofactorTree::CofactorTree(Abc_Ntk_t* pNtkCone) {
  assert(_pNtk_global);
  _root = new CofactorNode(pNtkCone);
  setGlobalPiReference();
  CofactorTree_rec(_root, true);
}

void CofactorTree::CofactorTree_rec(CofactorNode* n, bool root) {
  Abc_Ntk_t* pNtk = n->_pNtk;
  if(Abc_NtkNodeNum(pNtk) <= AigNodeThreshold) {
    Abc_Ntk_t* pNtkRes = Abc_NtkCollapse(pNtk, ABC_INFINITY, 0, 1, 0, 0, 0);
    Abc_Obj_t* pPi; int i;
    Abc_NtkForEachPi(pNtkRes, pPi, i) pPi->iTemp = Abc_NtkFindCi(_pNtk_global, Abc_ObjName(pPi))->iTemp;
    Abc_NtkForEachPi(_pNtk_global, pPi, i) assert(pPi->iTemp == i);
    n->_pNtk = pNtkRes;
    if(!root) Abc_NtkDelete(pNtk);
    return;
  }
  Abc_Obj_t* pPi, *pPi_min = 0;
  int i, min = __INT_MAX__, tmp;
  Abc_Ntk_t* pLNtk, * pRNtk, * pLNtk_min = 0, * pRNtk_min = 0;
  Abc_NtkForEachPi(pNtk, pPi, i){
    pLNtk = Cofactor(pNtk, true, i);
    pRNtk = Cofactor(pNtk, false, i);
    tmp = Abc_NtkNodeNum(pLNtk) + Abc_NtkNodeNum(pRNtk);
    if(tmp < min){
      if(i) {
        Abc_NtkDelete(pLNtk_min);
        Abc_NtkDelete(pRNtk_min);
      }
      pLNtk_min = pLNtk;
      pRNtk_min = pRNtk;
      pPi_min = pPi;
      min = tmp;
    }
    else{
      Abc_NtkDelete(pLNtk);
      Abc_NtkDelete(pRNtk);
    }
  }

  n->_dvar = Abc_NtkFindCi(_pNtk_global, Abc_ObjName(pPi_min))->iTemp;
  if(!root) Abc_NtkDelete(pNtk);
  n->_l = new CofactorNode(pLNtk_min);
  n->_r = new CofactorNode(pRNtk_min);
  CofactorTree_rec(n->_l);
  CofactorTree_rec(n->_r);
}

CofactorTree::~CofactorTree() {
  CofactorTree_Delete_rec(_root, true);
}

void CofactorTree::CofactorTree_Delete_rec(CofactorNode* n, bool root) {
  if(n->_l) {
    CofactorTree_Delete_rec(n->_l);
    delete n->_l;
    CofactorTree_Delete_rec(n->_r);
    delete n->_r;
  }
  else {
    if(!root) Abc_NtkDelete(n->_pNtk);
  }
}

Vec_Ptr_t* CofactorTree::toEsop() {
  Cube3* factor = Cube3Start(Abc_NtkPiNum(_pNtk_global));
  Vec_Ptr_t* cubes = Vec_PtrStart(0);
  toEsop_rec(_root, factor, cubes);
  Cube3Free(factor);
  return cubes;
}

void CofactorTree::toEsop_rec(CofactorNode* n, Cube3* factor, Vec_Ptr_t* cubes) {
  if(n->isLeaf()) {
    Abc_Ntk_t* pNtk = n->_pNtk;
    Abc_Obj_t* pNode = Abc_ObjFanin0(Abc_NtkPo(pNtk, 0));
    DdNode* bdd = (DdNode *)pNode->pData;
    if(Abc_ObjFaninC0(Abc_NtkPo(pNtk, 0))) bdd = Cudd_Complement(bdd);
    TDD tdd(bdd, pNtk);
    Cube3* tmp = Cube3Dup(factor);
    tdd.toEsop(tmp, cubes);
    Cube3Free(tmp);
    static int counter = 0;
    ++counter;
    if(counter == 1){
      printf("%d\n", bdd->index);
      // printf("%s\n", Cube3ToString(factor).c_str());
      // Abc_Obj_t* pNode = Abc_ObjFanin0(Abc_NtkPo(pNtk, 0));
      // if(Abc_ObjFaninC0(Abc_NtkPo(pNtk, 0))) printf("!\n");
      // Abc_NodeShowBdd(pNode, 0);
      // Vec_PtrClear(cubes);
      // tdd.toEsop(factor, cubes);
      // Cube3* cube; int k;
      // Vec_PtrForEachEntry(Cube3*, cubes, cube, k) printf("%s\n", Cube3ToString(cube).c_str());
      // Abc_Obj_t* pPi;
      // Abc_NtkForEachPi(pNtk, pPi, k) printf("%s %d\n", Abc_ObjName(pPi), pPi->iTemp);
      exit(0);
    }
    return;
  }
  Cube3WriteEntry(factor, n->_dvar, 1);
  toEsop_rec(n->_l, factor, cubes);
  Cube3WriteEntry(factor, n->_dvar, 0);
  toEsop_rec(n->_r, factor, cubes);
}


TDD::TDD(DdNode* n, Abc_Ntk_t* pNtk) {
  _root = new TDDNode(n);
  _pNtk = pNtk;
  Cudd_Ref(n);
  _nCubes = TDD_rec(_root);
}

int TDD::TDD_rec(TDDNode* tn) {
  DdNode* bdd = tn->_n;
  DdManager * dd = (DdManager *)_pNtk->pManFunc;
  if(Cudd_IsConstant(bdd)) return !Cudd_IsComplement(bdd);
  tn->_l = new TDDNode(Cudd_T(bdd));
  Cudd_Ref(Cudd_T(bdd));
  tn->_r = new TDDNode(Cudd_E(bdd));
  Cudd_Ref(Cudd_E(bdd));
  DdNode* x = Cudd_bddXor(dd, Cudd_T(bdd), Cudd_E(bdd));
  tn->_x = new TDDNode(x);
  Cudd_Ref(x);

  int cost_l = TDD_rec(tn->_l);
  int cost_r = TDD_rec(tn->_r);
  int cost_x = TDD_rec(tn->_x);

  TDDNode* deleteNode;

  if(cost_l > cost_r && cost_l > cost_x) {
    // Postive Davio
    deleteNode = tn->_l;
    tn->_l = 0;
    cost_l = 0;
  }
  else if(cost_r > cost_l && cost_r > cost_x) {
    // negative Davio
    deleteNode = tn->_r;
    tn->_r = 0;
    cost_r = 0;
  }
  else {
    // Shannon
    deleteNode = tn->_x;
    tn->_x = 0;
    cost_x = 0;
  }

  TDD_Delete_rec(deleteNode);

  return cost_l + cost_r + cost_x;
}

TDD::~TDD() {
  TDD_Delete_rec(_root);
}

void TDD::TDD_Delete_rec(TDDNode* tn) {
  DdManager * dd = (DdManager *)_pNtk->pManFunc;
  if(tn->_l) TDD_Delete_rec(tn->_l);
  if(tn->_r) TDD_Delete_rec(tn->_r);
  if(tn->_x) TDD_Delete_rec(tn->_x);
  Cudd_RecursiveDeref(dd, tn->_n);
  delete tn;
}

void TDD::toEsop(Cube3* cube, Vec_Ptr_t* cubes) {
  toEsop_rec(_root, cube, cubes);
}

void TDD::toEsop_rec(TDDNode* tn, Cube3* cube, Vec_Ptr_t* cubes) {
  DdNode* bdd = tn->_n;
  if(Cudd_IsConstant(tn->_n)) {
    if(!Cudd_IsComplement(tn->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    return;
  }
  int i = Abc_NtkPi(_pNtk, int(Cudd_Index(bdd)))->iTemp;
  Abc_Obj_t* pPi; int k;
  Abc_NtkForEachPi(CofactorTree::_pNtk_global, pPi, k) assert(pPi->iTemp == k);
  Abc_NtkForEachPi(_pNtk, pPi, k) assert(pPi->iTemp == Abc_NtkFindCi(CofactorTree::_pNtk_global, Abc_ObjName(pPi))->iTemp);
  if(!tn->_l) {
    // positive Davio
    assert(tn->_r && tn->_x);
    if(Cudd_IsComplement(tn->_r->_n) && !Cudd_IsConstant(tn->_r->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_r, cube, cubes);
    Cube3WriteEntry(cube, i, 1);
    if(Cudd_IsComplement(tn->_x->_n) && !Cudd_IsConstant(tn->_x->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_x, cube, cubes);
    Cube3WriteEntry(cube, i, 2);
  }
  else if(!tn->_r) {
    // negative Davio
    assert(tn->_l && tn->_x);
    if(Cudd_IsComplement(tn->_r->_n) && !Cudd_IsConstant(tn->_r->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_l, cube, cubes);
    Cube3WriteEntry(cube, i, 0);
    if(Cudd_IsComplement(tn->_l->_n) && !Cudd_IsConstant(tn->_l->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_l, cube, cubes);
    Cube3WriteEntry(cube, i, 2);
  }
  else if(!tn->_x){
    // Shannon
    assert(tn->_l && tn->_r);
    Cube3WriteEntry(cube, i, 1);
    if(Cudd_IsComplement(tn->_l->_n) && !Cudd_IsConstant(tn->_l->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_l, cube, cubes);
    Cube3WriteEntry(cube, i, 0);
    if(Cudd_IsComplement(tn->_r->_n) && !Cudd_IsConstant(tn->_r->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_r, cube, cubes);
    Cube3WriteEntry(cube, i, 2);
  }
  else {
    assert(0);
  }
}


int BDD_nCube(DdNode* n) {
  if(Cudd_IsConstant(n)) return !Cudd_IsComplement(n);
  return BDD_nCube(Cudd_T(n)) + BDD_nCube(Cudd_E(n));
}

Abc_Ntk_t* Cofactor(Abc_Ntk_t* pNtk, bool fPos, int iVar) {
  Abc_Ntk_t* pCof = Abc_NtkDup(pNtk), *pNtkRes;
  Vec_Ptr_t* vNodes;
  Abc_Obj_t *pObj, *pNext, *pFanin;
  int i;
  assert(Abc_NtkIsStrash(pCof));
  assert(iVar < Abc_NtkCiNum(pCof));

  // collect the internal nodes
  pObj = Abc_NtkCi(pCof, iVar);
  vNodes = Abc_NtkDfsReverseNodes(pCof, &pObj, 1);

  // assign the cofactors of the CI node to be constants
  pObj->pCopy = Abc_ObjNotCond(Abc_AigConst1(pCof), !fPos);

  // quantify the nodes
  Vec_PtrForEachEntry(Abc_Obj_t*, vNodes, pObj, i) {
    for (pNext = pObj ? pObj->pCopy : pObj; pObj; pObj = pNext, pNext = pObj ? pObj->pCopy : pObj) {
      pFanin = Abc_ObjFanin0(pObj);
      if (!Abc_NodeIsTravIdCurrent(pFanin)) {
        pFanin->pCopy = pFanin;
      }
      pFanin = Abc_ObjFanin1(pObj);
      if (!Abc_NodeIsTravIdCurrent(pFanin)) {
        pFanin->pCopy = pFanin;
      }
      pObj->pCopy = Abc_AigAnd((Abc_Aig_t*)pCof->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj));
    }
  }
  Vec_PtrFree(vNodes);

  // update the affected COs
  Abc_NtkForEachCo(pCof, pObj, i) {
    if (!Abc_NodeIsTravIdCurrent(pObj)) continue;
    pFanin = Abc_ObjFanin0(pObj);
    // get the result of quantification
    pNext = Abc_ObjNotCond(Abc_ObjChild0Copy(pObj), Abc_ObjFaninC0(pObj));
    if (Abc_ObjRegular(pNext) == pFanin) continue;
    // update the fanins of the CO
    Abc_ObjPatchFanin(pObj, pFanin, pNext);
  }

  pNtkRes = Abc_NtkStrash(pCof, 0, 1, 0);
  Abc_NtkDelete(pCof);
  return pNtkRes;
}
