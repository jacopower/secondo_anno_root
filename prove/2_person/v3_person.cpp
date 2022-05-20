#include <iostream>
#include <string>

class Person
{
public:
  Person(std::string const name);

  virtual void PrintData() const
  {
    std::cout << "Name: " << name_ << '\n';
  }

private:
  std::string const name_;
};

Person::Person(std::string const name) : name_{name} {}

class Student : public Person
{
public:
  Student(std::string const name, int code) : Person{name}, code_{code} {}
  void PrintData() const override; //è giusto override qui?
private:
  int code_;
};

void Student::PrintData() const
{
  Person::PrintData();
  std::cout << "Matricola: " << code_ << '\n';
}

int main()
{
  Person *test[2];
  test[0] = new Person("giovanni");
  test[1] = new Student("ciro", 34566);

  for (int i = 0; i < 2; ++i)
  {
    test[i]->PrintData();
  }
}