#include<iostream> 
using namespace std; 
  
// An abstract class with constructor 
class A 
{ 
protected: 
   int x; 
public: 
  virtual void fun() = 0; 
  A(int i) { x = i; } 
}; 
  
class B: public A 
{ 
    int y; 
public: 
    B(int i, int j):A(i), y(j) { } 
    void fun() 
    { 
      cout << "x = " << x << std::endl;
      cout << "y = " << y << std::endl;
    } 
}; 
  
class C: public B 
{ 
    int z; 
public: 
    C(int i, int j, int k):B(i, j), z(k){} 
    void fun() 
    { 
      B::fun();
      cout << "z= " << z << std::endl;
    } 
}; 

int main(void) 
{ 
    C d(4, 5, 6); 
    d.fun(); 
    return 0; 
} 
