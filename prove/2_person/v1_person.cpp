#include <iostream>

// RAPPORTO  IS-A

class Person {
  public:
  Person(char* const name);
  void PrintName() const;

  private:
  char* const name_;
};

Person::Person(char* const name) : name_{name} {}
void Person::PrintName() const {
  std::cout << "Name: " << name_ << '\n';
}

class Student {
  public:
  Student(char* const name, int code);
  void PrintName() const;

  private:
  Person self_;
  int code_;
};

Student::Student(char* const name, int code) : self_{name}, code_{code} {}
void Student::PrintName() const {
  self_.PrintName();
}

int main() {}