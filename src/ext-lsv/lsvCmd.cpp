#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"

static int Lsv_CommandPrintNodes(Abc_Frame_t* pAbc, int argc, char** argv);
static int Lsv_CommandPrintSopUnate(Abc_Frame_t* pAbc, int argc, char** argv);

void init(Abc_Frame_t* pAbc) {
  Cmd_CommandAdd(pAbc, "LSV", "lsv_print_nodes", Lsv_CommandPrintNodes, 0);
  Cmd_CommandAdd(pAbc, "LSV", "lsv_print_sopunate", Lsv_CommandPrintSopUnate, 0);
}

void destroy(Abc_Frame_t* pAbc) {}

Abc_FrameInitializer_t frame_initializer = {init, destroy};

struct PackageRegistrationManager {
  PackageRegistrationManager() { Abc_FrameAddInitializer(&frame_initializer); }
} lsvPackageRegistrationManager;

void Lsv_NtkPrintNodes(Abc_Ntk_t* pNtk) {
  Abc_Obj_t* pObj;
  int i;
  Abc_NtkForEachNode(pNtk, pObj, i) {
    printf("Object Id = %d, name = %s\n", Abc_ObjId(pObj), Abc_ObjName(pObj));
    Abc_Obj_t* pFanin;
    int j;
    Abc_ObjForEachFanin(pObj, pFanin, j) {
      printf("  Fanin-%d: Id = %d, name = %s\n", j, Abc_ObjId(pFanin),
             Abc_ObjName(pFanin));
    }
    if (Abc_NtkHasSop(pNtk)) {
      printf("The SOP of this node:\n%s", (char*)pObj->pData);
    }
  }
}

int Lsv_CommandPrintNodes(Abc_Frame_t* pAbc, int argc, char** argv) {
  Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
  int c;
  Extra_UtilGetoptReset();
  while ((c = Extra_UtilGetopt(argc, argv, "h")) != EOF) {
    switch (c) {
      case 'h':
        goto usage;
      default:
        goto usage;
    }
  }
  if (!pNtk) {
    Abc_Print(-1, "Empty network.\n");
    return 1;
  }
  Lsv_NtkPrintNodes(pNtk);
  return 0;

usage:
  Abc_Print(-2, "usage: lsv_print_nodes [-h]\n");
  Abc_Print(-2, "\t        prints the nodes in the network\n");
  Abc_Print(-2, "\t-h    : print the command usage\n");
  return 1;
}

int sort_id_compare(Abc_Obj_t** a, Abc_Obj_t** b) {
  if(Abc_ObjId(*a) > Abc_ObjId(*b)) return 1;
  else return -1;
}

void Lsv_PrintSopUnate(Abc_Ntk_t* pNtk) {
  Abc_Obj_t* node;
  int i;
  Abc_NtkForEachNode(pNtk, node, i) {
    if(!Abc_NtkHasSop(pNtk)) continue;
    if(!Abc_ObjFaninNum(node)) continue;
    printf("node %s:\n", Abc_ObjName(node));
    char* sop = (char*)node->pData;
    int unate_info_n = 0, unate_info_p = -1, phase_info_n = 0, phase_info_p = 0;
    int j = 0;
    int b;
    while(sop[j] != '\0') {
      if(sop[j] == ' ') {
        ++j;
        assert(sop[j] == '0' || sop[j] == '1');
        b = (int(sop[j]) - int('0'));
        phase_info_n = phase_info_n << 1;
        phase_info_p = phase_info_p << 1;
        phase_info_n |= b;
        phase_info_p |= b;
        unate_info_p &= phase_info_p;
        unate_info_n |= phase_info_n;
        phase_info_n = phase_info_p = 0;
        assert(sop[++j] == '\n');
        ++j;
      }
      else {
        if(sop[j] == '-') {
          phase_info_n = phase_info_n << 1;
          phase_info_p = phase_info_p << 1;
          phase_info_p |= 0x1;
        }
        else {
          b = (int(sop[j]) - int('0'));
          phase_info_n = phase_info_n << 1;
          phase_info_p = phase_info_p << 1;
          phase_info_n |= b;
          phase_info_p |= b;
        }
        ++j;
      }
    }

    int nFanins = Abc_ObjFaninNum(node);
    Vec_Ptr_t* unate_vars_n = Vec_PtrAlloc(nFanins), *unate_vars_p = Vec_PtrAlloc(nFanins), *binate_vars = Vec_PtrAlloc(nFanins);
    Abc_Obj_t* pFanin;
    int k;
    bool nu, pu;
    assert( (unate_info_n & 0x1 ) == (unate_info_p & 0x1) );
    bool onset = unate_info_n & 0x1;
    Abc_ObjForEachFanin(node, pFanin, k){
      nu = ~( ( unate_info_n >> (nFanins - k) ) | (-1 << 1) );
      pu = ( unate_info_p >> (nFanins - k) & 0x1 );
      if(onset) {
        if(nu) Vec_PtrPush(unate_vars_n, pFanin);
        if(pu) Vec_PtrPush(unate_vars_p, pFanin);
      }
      else {
        if(nu) Vec_PtrPush(unate_vars_p, pFanin);
        if(pu) Vec_PtrPush(unate_vars_n, pFanin);
      }
      if(!nu && !pu) Vec_PtrPush(binate_vars, pFanin);
    }
    Vec_PtrSort(unate_vars_n, (int(*)())sort_id_compare);
    Vec_PtrSort(unate_vars_p, (int(*)())sort_id_compare);
    Vec_PtrSort(binate_vars, (int(*)())sort_id_compare);

    Abc_Obj_t* entry;
    if(Vec_PtrSize(unate_vars_p)) {
      printf("+unate inputs: ");
      Vec_PtrForEachEntry(Abc_Obj_t*, unate_vars_p, entry, k){
        if(k) printf(",");
        printf("%s", Abc_ObjName(entry));
      }
      printf("\n");
    }
    if(Vec_PtrSize(unate_vars_n)) {
      printf("-unate inputs: ");
      Vec_PtrForEachEntry(Abc_Obj_t*, unate_vars_n, entry, k){
        if(k) printf(",");
        printf("%s", Abc_ObjName(entry));
      }
      printf("\n");
    }
    if(Vec_PtrSize(binate_vars)) {
      printf("binate inputs: ");
      Vec_PtrForEachEntry(Abc_Obj_t*, binate_vars, entry, k){
        if(k) printf(",");
        printf("%s", Abc_ObjName(entry));
      }
      printf("\n");
    }

    Vec_PtrFree(unate_vars_n);
    Vec_PtrFree(unate_vars_p);
    Vec_PtrFree(binate_vars);
  }
}

int Lsv_CommandPrintSopUnate(Abc_Frame_t* pAbc, int argc, char** argv) {
  Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
  int c;
  Extra_UtilGetoptReset();
  while ((c = Extra_UtilGetopt(argc, argv, "h")) != EOF) {
    switch (c) {
      case 'h':
        goto usage;
      default:
        goto usage;
    }
  }
  if (!pNtk) {
    Abc_Print(-1, "Empty network.\n");
    return 1;
  }
  Lsv_PrintSopUnate(pNtk);
  return 0;

usage:
  Abc_Print(-2, "usage: lsv_print_sopunate [-h]\n");
  Abc_Print(-2, "\t        prints the unate information for each node whose function is expressed in the SOP form\n");
  Abc_Print(-2, "\t-h    : print the command usage\n");
  return 1;
}