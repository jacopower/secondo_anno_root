#ifndef DERIVEDCLASS_HPP
#define DERIVEDCLASS_HPP

#include "BaseClass.hpp"

class DerivedClass : public BaseClass
{
public:
  DerivedClass(int baseMember, int derivedMember);
  int GetDerivedMember() const;

private:
  int derivedMember_;
};

#endif