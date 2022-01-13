#include "TF1.h"
#include "TH1F.h"

//class MyClass: public TObject{
class MyClass
{
public:
  MyClass();
  MyClass(TH1F *hIn);
  //public methods
  void Generate(TF1 *f, int nGen);
  void ShowHisto() const;

private:
  TH1F *h_;
  //for persistency
  //ClassDef(MyClass,1)
};
