#include <iostream>
#include "BaseClass.hpp"
#include "DerivedClass.hpp"

int main() {
  BaseClass myBase{3};
  DerivedClass myObject{1, 2};
  std::cout << "My object attributes (from base): " << myObject.GetBaseMember() << '\n';
  std::cout << "My object attributes (from derived): " << myObject.GetDerivedMember() << '\n';
  std::cout << "Counting my objects: " << BaseClass::GetBaseCounter() << '\n';
}