#include <iostream>
#include "Person.h"
#include "Student.h"
using namespace std;
int main()
{
  Person *test[2];
  test[0] = new Person("tizio");
  test[1] = new Student("caio", 223344);
  for (int i = 0; i < 2; i++)
    test[i]->PrintData();
  return 0;
}
