#include <vector>

#include "TestMacro.h"

class test_obj
{
public:
  test_obj(): x(std::vector<int>({1, 2 ,3})){}

  void copy(std::vector<int> & y)
  {
    y = x; //y 是 x 的拷贝
  }

  std::vector<int> x;
};

void test_ref()
{
  test_obj x;
  std::vector<int> y;
  x.copy(y);
  ASSERT_EQUAL(x.x[0], y[0]);
  ASSERT_EQUAL(x.x[1], y[1]);
  ASSERT_EQUAL(x.x[2], y[2]);
  //x.x = std::vector<int>({4, 5, 6});
  ASSERT_EQUAL(x.x[0], y[0]);
  ASSERT_EQUAL(x.x[1], y[1]);
  ASSERT_EQUAL(x.x[2], y[2]);
}

int main()
{
  test_ref();
}
