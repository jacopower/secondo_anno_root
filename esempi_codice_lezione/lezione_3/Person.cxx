#include <iostream>
#include "Person.h"
Person::Person(const char *name) : name_(name) {}
void Person::PrintName() const { std::cout << "Nome:" << name_ << std::endl; }
void Person::PrintData() const { std::cout << "Nome:" << name_ << std::endl; }
