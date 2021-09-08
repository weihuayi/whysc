#include <array>
#include <vector>

#include "TestMacro.h"

void build_array()
{
  auto a = std::array<int, 3>({0, 1, 2});
  ASSERT_EQUAL(0, a[0]);
  ASSERT_EQUAL(1, a[1]);
  ASSERT_EQUAL(2, a[2]);
}

void build_array_in_vector()
{
  auto a = std::vector<std::array<int, 3> >({{0, 1, 2}, {0, 1, 2}});
  ASSERT_EQUAL(2, a.size());
  ASSERT_EQUAL(0, a[0][0]);
  ASSERT_EQUAL(1, a[0][1]);
  ASSERT_EQUAL(2, a[0][2]);
  ASSERT_EQUAL(0, a[1][0]);
  ASSERT_EQUAL(1, a[1][1]);
  ASSERT_EQUAL(2, a[1][2]);
}

int main()
{
  build_array();
  build_array_in_vector();
}
