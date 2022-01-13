#include <iostream>

class Person
{
public:
  Person(char *const name);
  void PrintName() const;

private:
  char *const name_;
};

Person::Person(char *const name) : name_{name} {}
void Person::PrintName() const
{
  std::cout << "Name: " << name_ << '\n';
}

class Student : public Person
{
public:
  Student(char *const name, int code);

private:
  int code_;
};

Student::Student(char *const name, int code) : Person{name}, code_{code} {}

int main() {}