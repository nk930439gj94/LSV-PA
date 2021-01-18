#include "lsvEsop.h"

Abc_Ntk_t* CofactorTree::_pNtk_global = 0;

void CofactorTree::setGlobalNtk(Abc_Ntk_t* pNtk_global) {
  _pNtk_global = pNtk_global;
  Abc_Obj_t* pPi; int i;
  Abc_NtkForEachPi(pNtk_global, pPi, i) pPi->iData = i;
}

CofactorTree::CofactorTree(Abc_Ntk_t* pNtkCone) {
  assert(_pNtk_global);
  _root = new CofactorNode(pNtkCone);
  CofactorTree_rec(_root, true);
}

void CofactorTree::CofactorTree_rec(CofactorNode* n, bool root) {
  Abc_Ntk_t* pNtk = n->_pNtk;
  if(Abc_NtkNodeNum(pNtk) <= AigNodeThreshold) {
    Abc_Ntk_t* pNtkRes = Collapse_reservePi(pNtk, 0);
#ifdef debug
    Abc_Obj_t* pPi; int i;
    Abc_NtkForEachPi(_pNtk_global, pPi, i) assert(pPi->iData == i);
#endif
    setGlobalPiReference(pNtkRes);
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

  n->_dvar = Abc_NtkFindCi(_pNtk_global, Abc_ObjName(pPi_min))->iData;
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
  esopSimplify(cubes);
  return cubes;
}

void CofactorTree::toEsop_rec(CofactorNode* n, Cube3* factor, Vec_Ptr_t* cubes) {
  if(n->isLeaf()) {
    Abc_Ntk_t* pNtk = n->_pNtk;
    Abc_Obj_t* pNode = Abc_ObjFanin0(Abc_NtkPo(pNtk, 0));
    DdNode* bdd = (DdNode *)pNode->pData;
    if(Abc_ObjFaninC0(Abc_NtkPo(pNtk, 0))) bdd = Cudd_Complement(bdd);
    TDD tdd(bdd, pNtk);
    tdd.toEsop(factor, cubes);
    return;
  }
  Cube3WriteEntry(factor, n->_dvar, 1);
  toEsop_rec(n->_l, factor, cubes);
  Cube3WriteEntry(factor, n->_dvar, 0);
  toEsop_rec(n->_r, factor, cubes);
  Cube3WriteEntry(factor, n->_dvar, 2);
}

void CofactorTree::setGlobalPiReference(Abc_Ntk_t* pNtk) {
  assert(_pNtk_global);
  Abc_Obj_t* pPi; int i;
  Abc_NtkForEachPi(pNtk, pPi, i) pPi->iData = Abc_NtkFindCi(_pNtk_global, Abc_ObjName(pPi))->iData;
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
  if(Cudd_IsComplement(_root->_n) && !Cudd_IsConstant(_root->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
  toEsop_rec(_root, cube, cubes);
}

void TDD::toEsop_rec(TDDNode* tn, Cube3* cube, Vec_Ptr_t* cubes) {
  DdNode* bdd = tn->_n;
  if(Cudd_IsConstant(tn->_n)) {
    if(!Cudd_IsComplement(tn->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    return;
  }
  int i = Abc_NtkPi(_pNtk, int(Cudd_Index(bdd)))->iData;
#ifdef debug
  Abc_Obj_t* pPi; int k;
  Abc_NtkForEachPi(CofactorTree::_pNtk_global, pPi, k) assert(pPi->iData == k);
  Abc_NtkForEachPi(_pNtk, pPi, k) assert(pPi->iData == Abc_NtkFindCi(CofactorTree::_pNtk_global, Abc_ObjName(pPi))->iData);
#endif
  if(!tn->_l) {
    // positive Davio
    if(Cudd_IsComplement(tn->_r->_n) && !Cudd_IsConstant(tn->_r->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_r, cube, cubes);
    Cube3WriteEntry(cube, i, 1);
    if(Cudd_IsComplement(tn->_x->_n) && !Cudd_IsConstant(tn->_x->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_x, cube, cubes);
    Cube3WriteEntry(cube, i, 2);
  }
  else if(!tn->_r) {
    // negative Davio
    if(Cudd_IsComplement(tn->_r->_n) && !Cudd_IsConstant(tn->_r->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_l, cube, cubes);
    Cube3WriteEntry(cube, i, 0);
    if(Cudd_IsComplement(tn->_l->_n) && !Cudd_IsConstant(tn->_l->_n)) Vec_PtrPush(cubes, Cube3Dup(cube));
    toEsop_rec(tn->_l, cube, cubes);
    Cube3WriteEntry(cube, i, 2);
  }
  else if(!tn->_x){
    // Shannon
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


Abc_Ntk_t * Collapse_reservePi( Abc_Ntk_t * pNtk, int fReorder )
{
    Abc_Ntk_t * pNtkNew;

    assert( Abc_NtkIsStrash(pNtk) );
    // compute the global BDDs
    if ( Abc_NtkBuildGlobalBdds(pNtk, ABC_INFINITY, 1, fReorder, 0, 0) == NULL )
        return NULL;

    // create the new network
    pNtkNew = Abc_NtkFromGlobalBdds( pNtk, 0 );
    Abc_NtkFreeGlobalBdds( pNtk, 1 );
    if ( pNtkNew == NULL )
        return NULL;

    if ( pNtk->pExdc )
        pNtkNew->pExdc = Abc_NtkDup( pNtk->pExdc );

    // make sure that everything is okay
    if ( !Abc_NtkCheck( pNtkNew ) )
    {
        printf( "Abc_NtkCollapse: The network check has failed.\n" );
        Abc_NtkDelete( pNtkNew );
        return NULL;
    }
    return pNtkNew;
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
  Abc_NtkRewrite(pNtkRes, 1, 0, 0, 0, 0);
  return pNtkRes;
}

void esopStats(Vec_Ptr_t* cubes) {
  int n_cubes = Vec_PtrSize(cubes);
  int n_literals = 0;
  Cube3* cube; int i;
  Vec_PtrForEachEntry(Cube3*, cubes, cube, i) n_literals += Cube3CountLiteral(cube);
  printf("Cube #   : %d\n", n_cubes);
  printf("Literal #: %d\n", n_literals);
}

void esopSimplify(Vec_Ptr_t* cubes) {
  int size = Vec_PtrSize(cubes);
  assert(size);
  Cube3* cube0, * cube1;
  int i, j, dist;
  for(i = 0; i < size - 1; ++i) {
    cube0 = (Cube3*)Vec_PtrEntry(cubes, i);
    if(!cube0) continue;
    for(j = i+1; j < size; ++j) {
      cube1 = (Cube3*)Vec_PtrEntry(cubes, j);
      if(!cube1) continue;
      dist = Cube3Distance(cube0, cube1);
      if(dist == 0) {
        Vec_PtrWriteEntry(cubes, i, 0);
        Vec_PtrWriteEntry(cubes, j, 0);
        Cube3Free(cube0);
        Cube3Free(cube1);
        break;
      }
      else if(dist == 1) {
        Cube3MergeDist1(cube0, cube1);
        Vec_PtrWriteEntry(cubes, j, 0);
        Cube3Free(cube1);
        break;
      }
    }
  }
  j = 0;
  for(i = 0; i < size; ++i) {
    cube0 = (Cube3*)Vec_PtrEntry(cubes, i);
    if(cube0) 
      Vec_PtrWriteEntry(cubes, j++, cube0);
  }
  Vec_PtrShrink(cubes, j);
}

void esopPrint(Vec_Ptr_t* cubes, Vec_Ptr_t* PiNames) {
  Cube3* cube; int i;
  Vec_PtrForEachEntry(Cube3*, cubes, cube, i) printf("%s\n", Cube3ToString(cube, PiNames).c_str());
}

void esopFree(Vec_Ptr_t* cubes) {
  Cube3* cube; int i;
  Vec_PtrForEachEntry(Cube3*, cubes, cube, i) Cube3Free(cube);
  Vec_PtrFree(cubes);
}
