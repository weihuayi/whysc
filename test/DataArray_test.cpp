#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "mesh/ComponentDataArray.h"

typedef WHYSC::Mesh::ComponentDataArray<double> CDataArray;


int main(int argc, char **argv)
{
  CDataArray array("test", 3);

  array.resize(10);
  std::cout << array.size() << std::endl;
  std::cout << array.number_of_components() << std::endl;

  return 0;
}
