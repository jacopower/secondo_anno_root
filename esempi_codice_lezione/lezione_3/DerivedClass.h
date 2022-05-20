#ifndef DERIVEDCLASS_H
#define DERIVEDCLASS_H
#include "BaseClass.h"
class DerivedClass : public BaseClass
{
public:
  DerivedClass(int baseMember, int derivedMember);
  int GetDerivedMember() const;

private:
  int derivedMember_;
};
#endif
