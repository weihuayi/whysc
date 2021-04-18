#include <fstream>
#include<stdlib.h>
#include <string>
#include <vector>
#include <iostream>

int main(int argc, char* argv[])
{
  const char* filename0 = (argc > 1) ? argv[1] : "data/pig.off";
  const char* filename1 = (argc > 2) ? argv[2] : "data/pig.off";
  std::ifstream input;
  input.open(filename0);

	if (!input.is_open())
	{
    std::cout<< "读取文件失败" <<std::endl;
		return 0;
	}
  std::vector<std::vector<double>> vdata;
  std::vector<std::vector<int>> cdata;
  std::string data;
  int vsize = 1022;
  int csize = 2040;

  int i = 0;
  while(getline(input, data))
  {
    if(data[0] == 'f')
    {
      continue;
    }

    if(i < vsize) 
    {
      std::vector<double> v;
      std::string x="";
      std::cout<< "data" << data <<std::endl;
      for(int j = 0; j < data.length(); j++)
      {
        std::cout<< x <<std::endl;
        if(data[j] !=' ')
        {
          x+=data[j];
        }

        if((data[j] == ' ' && x != "") || j == data.length()-1)
        {
          v.push_back(stod(x));
          std::cout<< x <<std::endl;
          x="";
        }
      }
      vdata.push_back(v);
    } 
    else 
    {
      std::vector<int> c;
      std::string x="";
      std::cout<< "data" << data <<std::endl;
      for(int j = 0; j < data.length(); j++)
      {
        std::cout<< x <<std::endl;
        if(data[j] !=' ')
        {
          x+=data[j];
        }

        if((data[j] == ' ' && x != "") || j == data.length()-1)
        {
          c.push_back(stod(x));
          std::cout<< x <<std::endl;
          x="";
        }
      }
      cdata.push_back(c);
    }

    i++;
  }
  std::cout<< csize <<std::endl;
  input.close();

  int vert_number = vdata.size();
  int face_number = cdata.size();

  std::ofstream out(filename1);
  out << "OFF"<< std::endl;
  out << vdata.size() <<" "<< cdata.size() <<" 0"<< std::endl;
  for (int j = 0; j <vdata.size() ; ++j) {
      out << (vdata[j][0]) <<" " << vdata[j][1]<< " "<< vdata[j][2]<<std::endl;
  }
  for (int j = 0; j <cdata.size() ; ++j) {
      out <<"3 " << (cdata[j][0]) <<" " << cdata[j][1] <<" " << cdata[j][2] <<std::endl;
  }
  out.close();
}
