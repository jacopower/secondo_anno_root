#ifndef MYCLASS_HPP
#define MYCLASS_HPP

#include "TH1F.h"
#include "TF1.h"

class myClass : public TObject
{
  public:
  myClass();
  myClass(TH1F *hIn);
  void Generate(TF1 *fIn, Int_t nGen);
  void ShowHisto() const;

  private:
  TH1F *h_;

  ClassDef(myClass, 1);
};

#endif