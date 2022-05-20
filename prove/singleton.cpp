class MySingleton
{
public:
  ~MySingleton() { ptr_ = 0; }
  static MySingleton *instance()
  {
    if (ptr_ == 0)
    {
      ptr_ = new MySingleton();
      return ptr_;
    }
  }

  // altri metodi

private:
  MySingleton() {}
  static MySingleton *ptr_;
};

MySingleton *MySingleton::ptr_ = 0; //initialize static member

int main() {
  MySingleton *p = MySingleton::instance();
  //p->method();
  delete p;
}