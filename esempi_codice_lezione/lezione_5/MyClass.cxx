#include "MyClass.h"

MyClass::MyClass() { h_ = new TH1F("h", "titolo", 100, -5, 5); }
MyClass::MyClass(TH1F *hIn) { h_ = new TH1F(*hIn); }

void MyClass::Generate(TF1 *f, int nGen)
{
  h_->FillRandom(f->GetName(), nGen);
}
void MyClass::ShowHisto() const { h_->Draw(); }
//for persistency
//ClassImp(MyClass)
