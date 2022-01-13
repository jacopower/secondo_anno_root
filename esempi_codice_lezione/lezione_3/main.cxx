#include <iostream>
#include "BaseClass.h"
#include "DerivedClass.h"
using namespace std;
int main()
{
  BaseClass myBase(3);
  DerivedClass myObject(1, 2);
  cout << " my object attribute 1 : " << myObject.GetBaseMember() << endl;
  cout << " my object attribute 2 : " << myObject.GetDerivedMember() << endl;
  cout << "counting objects:" << BaseClass::GetBaseCounter() << endl;
  return 0;
}
