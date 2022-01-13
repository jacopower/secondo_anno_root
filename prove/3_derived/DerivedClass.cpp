#include "DerivedClass.hpp"

DerivedClass::DerivedClass(int baseMember, int derivedMember) : BaseClass{baseMember}, derivedMember_{derivedMember} {}

int DerivedClass::GetDerivedMember() const { return derivedMember_; }
