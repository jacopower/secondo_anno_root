#ifndef PERSON_H
#define PERSON_H
class Person
{
public:
  Person() = default;
  Person(const char *name);
  void PrintName() const;
  virtual void PrintData() const;

private:
  const char *name_;
};
#endif
