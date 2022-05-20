#include "BaseClass.hpp"

int BaseClass::baseCounter_ = 0;

BaseClass::BaseClass(int baseMember) : baseMember_{baseMember} { baseCounter_++; }

int BaseClass::GetBaseMember() const { return baseMember_; }

int BaseClass::GetBaseCounter() { return baseCounter_; }

