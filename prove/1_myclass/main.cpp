#include <iostream>
#include "myClass.hpp"

int main()
{
  MyClass myObject{1};
  std::cout << myObject.GetA() << '\n';
}