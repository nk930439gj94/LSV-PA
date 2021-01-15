#include "lsvEsop.h"

CofactorTree::CofactorTree(Abc_Ntk_t* pNtkCone, Abc_Ntk_t* pNtk_origin) {
  _root = new CofactorNode(pNtkCone);
  level = CofactorTree_rec(_root, pNtk_origin, true);
}

int CofactorTree::CofactorTree_rec(CofactorNode* n, Abc_Ntk_t* pNtk_origin, bool root) {
  Abc_Ntk_t* pNtk = n->_pNtk;
  if(Abc_NtkNodeNum(pNtk) <= AigNodeThreshold) {
    Abc_Ntk_t* pNtkRes = Abc_NtkCollapse(pNtk, ABC_INFINITY, 0, 1, 0, 0, 0);
    n->_pNtk = pNtkRes;
    if(!root) Abc_NtkDelete(pNtk);
    return 1;
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

  n->_dvar = Abc_NtkFindCi(pNtk_origin, Abc_ObjName(pPi_min));
  if(!root) Abc_NtkDelete(pNtk);
  n->_l = new CofactorNode(pLNtk_min);
  n->_r = new CofactorNode(pRNtk_min);
  int l_level = CofactorTree_rec(n->_l, pNtk_origin);
  int r_level = CofactorTree_rec(n->_r, pNtk_origin);
  if(l_level > r_level) return l_level + 1;
  else return r_level + 1;
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

TDD::TDD(DdNode* n, DdManager* dd) {
  _root = new TDDNode(n);
  Cudd_Ref(n);
  TDD_rec(_root, dd);
}

int TDD::TDD_rec(TDDNode* tn, DdManager* dd) {
  DdNode* bdd = tn->_n;
  if(Cudd_IsConstant(bdd)) return !Cudd_IsComplement(bdd);
  tn->_l = new TDDNode(Cudd_T(bdd));
  Cudd_Ref(Cudd_T(bdd));
  tn->_r = new TDDNode(Cudd_E(bdd));
  Cudd_Ref(Cudd_E(bdd));
  DdNode* x = Cudd_bddXor(dd, Cudd_T(bdd), Cudd_E(bdd));
  tn->_x = new TDDNode(x);
  Cudd_Ref(x);

  int cost_l = TDD_rec(tn->_l, dd);
  int cost_r = TDD_rec(tn->_r, dd);
  int cost_x = TDD_rec(tn->_x, dd);

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

  if(deleteNode->_l) {
    Cudd_RecursiveDeref(dd, deleteNode->_l->_n);
    delete deleteNode->_l;
  }
  if(deleteNode->_r) {
    Cudd_RecursiveDeref(dd, deleteNode->_r->_n);
    delete deleteNode->_r;
  }
  if(deleteNode->_x) {
    Cudd_RecursiveDeref(dd, deleteNode->_x->_n);
    delete deleteNode->_x;
  }
  Cudd_RecursiveDeref(dd, deleteNode->_n);
  delete deleteNode;

  return cost_l + cost_r + cost_x;
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
