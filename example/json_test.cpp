#include <iostream>
#include <sstream>
#include <vector>
#include "thirdparty/json.hpp"

using json = nlohmann::json;

int main(int argc, char * argv[])
{
  // create an empty structure (null)
  json j;

  // add a number that is stored as double (note the implicit conversion of j to an object)
  j["pi"] = 3.141;

  // add a Boolean that is stored as bool
  j["happy"] = true;

  // add a string that is stored as std::string
  j["name"] = "Niels";

  // add another null object by passing nullptr
  j["nothing"] = nullptr;

  // add an object inside the object
  j["answer"]["everything"] = 42;

  // add an array that is stored as std::vector (using an initializer list)
  j["list"] = { 1, 0, 2 };

  // add another object (using an initializer list of pairs)
  j["object"] = { {"currency", "USD"}, {"value", 42.99} };

  j["2"] = std::vector<int>(3);
  j["pds"]["2"] = std::vector<int>(3);

  auto & v = j["list"];
  v[0] = 3;

  for(auto i : v)
    std::cout << i << std::endl;


  for(auto i : j)
    std::cout << i << "ceshi" << std::endl;

  // instead, you could also write (which looks very similar to the JSON above)
  json j2 = {
    {"pi", 3.141},
    {"happy", true},
    {"name", "Niels"},
    {"nothing", nullptr},
    {"answer", {
      {"everything", 42}
    }},
    {"list", {1, 0, 2}},
    {"object", {
      {"currency", "USD"},
      {"value", 42.99}
    }}
  };
  return 0;
}
