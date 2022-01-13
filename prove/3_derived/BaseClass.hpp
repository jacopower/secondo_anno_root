#ifndef BASECLASS_HPP
#define BASECLASS_HPP

class BaseClass
{
public:
  BaseClass(int baseMember);
  int GetBaseMember() const;
  static int GetBaseCounter();

private:
  int baseMember_;
  static int baseCounter_;
};

#endif