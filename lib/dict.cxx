// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libdIdict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR_Dictionary();
   static void ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR_TClassManip(TClass*);
   static void *new_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p = nullptr);
   static void *newArray_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(Long_t size, void *p);
   static void delete_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p);
   static void deleteArray_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p);
   static void destruct_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ROOT::VecOps::RVec<TLorentzVector>*)
   {
      ROOT::VecOps::RVec<TLorentzVector> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(ROOT::VecOps::RVec<TLorentzVector>));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::VecOps::RVec<TLorentzVector>", -2, "ROOT/RVec.hxx", 1530,
                  typeid(ROOT::VecOps::RVec<TLorentzVector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR_Dictionary, isa_proxy, 4,
                  sizeof(ROOT::VecOps::RVec<TLorentzVector>) );
      instance.SetNew(&new_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR);
      instance.SetNewArray(&newArray_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR);
      instance.SetDelete(&delete_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR);
      instance.SetDestructor(&destruct_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< ROOT::VecOps::RVec<TLorentzVector> >()));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ROOT::VecOps::RVec<TLorentzVector>*)
   {
      return GenerateInitInstanceLocal(static_cast<ROOT::VecOps::RVec<TLorentzVector>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ROOT::VecOps::RVec<TLorentzVector>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ROOT::VecOps::RVec<TLorentzVector>*>(nullptr))->GetClass();
      ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ROOT::VecOps::RVec<TLorentzVector> : new ROOT::VecOps::RVec<TLorentzVector>;
   }
   static void *newArray_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ROOT::VecOps::RVec<TLorentzVector>[nElements] : new ROOT::VecOps::RVec<TLorentzVector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p) {
      delete (static_cast<ROOT::VecOps::RVec<TLorentzVector>*>(p));
   }
   static void deleteArray_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p) {
      delete [] (static_cast<ROOT::VecOps::RVec<TLorentzVector>*>(p));
   }
   static void destruct_ROOTcLcLVecOpscLcLRVeclETLorentzVectorgR(void *p) {
      typedef ROOT::VecOps::RVec<TLorentzVector> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ROOT::VecOps::RVec<TLorentzVector>

namespace ROOT {
   static TClass *ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR_Dictionary();
   static void ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR_TClassManip(TClass*);
   static void *new_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p = nullptr);
   static void *newArray_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(Long_t size, void *p);
   static void delete_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p);
   static void deleteArray_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p);
   static void destruct_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*)
   {
      ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >", -2, "ROOT/RVec.hxx", 1530,
                  typeid(ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >) );
      instance.SetNew(&new_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR);
      instance.SetNewArray(&newArray_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR);
      instance.SetDelete(&delete_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR);
      instance.SetDestructor(&destruct_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > >()));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*>(nullptr))->GetClass();
      ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > : new ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >;
   }
   static void *newArray_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >[nElements] : new ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p) {
      delete (static_cast<ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*>(p));
   }
   static void deleteArray_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p) {
      delete [] (static_cast<ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >*>(p));
   }
   static void destruct_ROOTcLcLVecOpscLcLRVeclEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPtEtaPhiM4DlEdoublegRsPgRsPgR(void *p) {
      typedef ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_dict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_dict) = RecordReadRules_dict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_dict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_dict_Impl() {
    static const char* headers[] = {
nullptr
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/usr/include/root",
"/eos/home-t/thkelly/code/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace ROOT{namespace Math{template <typename T> class __attribute__((annotate("$clingAutoload$Math/Vector4Dfwd.h")))  __attribute__((annotate("$clingAutoload$TLorentzVector.h")))  PtEtaPhiM4D;
}}
namespace ROOT{namespace Math{template <class CoordSystem> class __attribute__((annotate("$clingAutoload$Math/Vector4Dfwd.h")))  __attribute__((annotate("$clingAutoload$TLorentzVector.h")))  LorentzVector;
}}
namespace ROOT{namespace VecOps{template <typename T> class __attribute__((annotate(R"ATTRDUMP(__cling__ptrcheck(off))ATTRDUMP"))) __attribute__((annotate("$clingAutoload$ROOT/RVec.hxx")))  RVec;
}}
class __attribute__((annotate("$clingAutoload$TLorentzVector.h")))  TLorentzVector;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"

// rootcling -f lib/dict.cxx -I$(root-config --incdir) lib/linkdef.h
// g++ -shared -fPIC lib/dict.cxx `root-config --cflags --libs` -o lib/libdict.so

#ifdef __CLING__
#pragma link C++ class ROOT::VecOps::RVec<TLorentzVector>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector>+;
#endif


#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >", payloadCode, "@",
"ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector>", payloadCode, "@",
"ROOT::VecOps::RVec<TLorentzVector>", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict() {
  TriggerDictionaryInitialization_dict_Impl();
}
