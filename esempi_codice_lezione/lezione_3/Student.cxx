#include <iostream>
#include "Student.h"
Student::Student(const char *name,
                 int code) : Person(name), code_(code) {}
void Student::PrintData() const
{
  Person::PrintData();
  std::cout << "Matricola:" << code_ << std::endl;
}
