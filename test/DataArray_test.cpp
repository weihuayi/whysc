#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "mesh/ComponentDataArray.h"
#include "mesh/OffsetDataArray.h"

typedef WHYSC::Mesh::ComponentDataArray<double> CDataArray;
typedef WHYSC::Mesh::OffsetDataArray<int> OffsetDataArray;


int main(int argc, char **argv)
{
  CDataArray a0("test", 3);
  a0.resize(10);
  std::cout << a0.size() << std::endl;
  std::cout << a0.number_of_components() << std::endl;

  OffsetDataArray a1("cell"); 
  auto& offset = a1.offset();
  offset.resize(3);
  offset[0] = 0;
  offset[1] = 3;
  offset[2] = 6;

  a1.resize(6);
  a1[0] = 0;
  a1[1] = 3;
  a1[2] = 4;
  a1[3] = 4;
  a1[4] = 1;
  a1[5] = 0;
  std::cout << a1.size() << std::endl;
  std::cout << a1.number_of_components() << std::endl;

  return 0;
}
