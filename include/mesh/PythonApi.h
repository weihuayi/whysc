#include "Python.h"
#include <string>

class PythonApi
{
public:
  PythonApi(std::string fname, int n, char** argv)
  {
    for(auto key : argv)
    {
      if(key=='int_list')

    }
  }

private:
  PyObject* m_ptuple;

};
int main()
{
  Py_Initialize();
	if( !Py_IsInitialized()){
		cout << "python init fail" << endl;
		return 0;
	}
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('./test')");

	PyObject* pModule = PyImport_ImportModule("sayhallo");

	if( pModule == NULL ){
		cout <<"module not found" << endl;
		return 1;
	}
  std::cout<< pModule <<std::endl;

	PyObject* pH = PyObject_GetAttrString(pModule, "Histogram_plot");

	if(!pH || !PyCallable_Check(pH)){
		cout <<"not found function add_num" << endl;
		return 0;
	}

  PyObject* plist = PyList_New(1000);
  PyObject* ptuple = PyTuple_New(1);

  for(int i = 0; i < 1000; i++)
  {
    double ra = pow((double)rand()*1.0/RAND_MAX, 2);
    PyObject*  pra = Py_BuildValue("d", ra);
    PyList_SetItem(plist, i, pra);
  }
  PyTuple_SetItem(ptuple, 0, plist);
  PyObject *Item = PyList_GetItem(plist, 4);//获取List对象中的每一个元素
  double result;
  PyArg_Parse(Item, "d", &result);
  std::cout<< result <<std::endl;

	PyObject_CallObject(pH, ptuple);
	Py_Finalize();
  return 0;
};
