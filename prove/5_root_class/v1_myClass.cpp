#include "myClass.hpp"

myClass::myClass()
{
  h_ = new TH1F("h_", "title", 100, -5., 5.);
}

myClass::myClass(TH1F *hIn)
{
  h_ = new TH1F(*hIn);
}

void myClass::Generate(TF1 *fIn, Int_t nGen)
{
  h_->FillRandom(fIn->GetName(), nGen);
}

void myClass::ShowHisto() const
{
  h_->Draw();
}
