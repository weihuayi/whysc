#include <iostream>
#include <algorithm>
#include "Matrix.h"

namespace WHYSC {
namespace AlgebraObject {

template<typename F=double, typename I=int>
struct Eigenvalue
{
    float *data, *eig;
    
    template<typename Matrix> 
    Eigenvalue(Matrix & A)
    {
        eig = Eig(A);

    }
    void Eig(data)
    { 
        I S0 = data.shape[0];
        I S1 = data.shape[1];
        for(auto i=0; i<S0;i++)
        {
            for(auto j=0; j<S1; j++)
            {
                eig = data[i][j];
            }
        }
    }
};

}




































}
