#include <algorithm>
#include <vector>
#include <set>
#include <iostream>


void test_minus()
{
  std::set<int> a({1, 2, 3, 4, 5});
  std::set<int> b({3, 4, 6, 10});
  std::set<int> c({1, 2, 5});
  std::vector<int> e(10);
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(), e.begin());
  std::cout<< e[0] << " " << e[1] << " " << e[2] <<std::endl;
}

int main()
{
  test_minus();
}
