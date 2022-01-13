#ifndef BASECLASS_H
#define BASECLASS_H
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
