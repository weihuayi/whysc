#ifndef HEXAHEDRON_H
#define HEXAHEDRON_H

#include "Point_3.h"
#include "Vector_3.h"
#include "constant.h"

#include <vector>
#include <cmath>

namespace iMath {
namespace GeometryObject {

class Hexahedron
{
    static int face[6][4];
    static int corner[8][3];

public:
    Hexahedron()
    {
        data[0] = Point_3(0.0, 0.0, 0.0);
        data[1] = Point_3(1.0, 0.0, 0.0);
        data[2] = Point_3(1.0, 1.0, 0.0);
        data[3] = Point_3(0.0, 1.0, 0.0);
        data[4] = Point_3(0.0, 0.0, 1.0);
        data[5] = Point_3(1.0, 0.0, 1.0);
        data[6] = Point_3(1.0, 1.0, 1.0);
        data[7] = Point_3(0.0, 1.0, 1.0);
    }

    Hexahedron(double h)
    {
        data[0] = Point_3(0.0, 0.0, 0.0);
        data[1] = Point_3(h, 0.0, 0.0);
        data[2] = Point_3(h,  h,  0.0);
        data[3] = Point_3(0.0, h, 0.0);
        data[4] = Point_3(0.0, 0.0, h);
        data[5] = Point_3(h, 0.0, h);
        data[6] = Point_3(h, h, h);
        data[7] = Point_3(0.0, h, h);
    }

    Point_3 & operator [] (const int i) 
    {
        return data[i];
    }

    const Point_3 & operator [] (const int i) const
    {
        return data[i];
    }

    Point_3 & vertex(int i)
    {
        return data[i];
    }

    double getCornerJacobi(int i)
    {
        int  j = corner[i][0];
        int  k = corner[i][1];
        int m = corner[i][2];
        Vector_3 vij = data[j] - data[i];
        Vector_3 vik = data[k] - data[i];
        Vector_3 vim = data[m] - data[i];

        double lij = vij.squaredLength();
        double lik = vik.squaredLength();
        double lim = vim.squaredLength();
        return 3*cross(vij, vik)*vim/(lij*std::sqrt(lij) + lik*std::sqrt(lik) + lim*std::sqrt(lim));
    }

    int getNumberOfMinusJacobi()
    {
        int num = 0;
        for(int i = 0; i < 8; i++)
        {
            double J = getCornerJacobi(i);
            if(J < EPS)
                num++;
        }
        return num;
    }
    
    bool isValid()
    {
        for(int i = 0; i < 8; i++)
        {
            double J = getCornerJacobi(i);
            if(J < EPS)
                return false;
        }
        return true;
    }

    double quality()
    {
        double sumq = 0.0;
        for(int i=0; i< 8; i++)
        {
            sumq += quality(i);
        }

        return sumq/8.0;
    }

    double quality(int i)
    {
        assert( i >= 0 && i < 8);
        int j = corner[i][0];
        int k = corner[i][1];
        int m = corner[i][2];
        Vector_3 vji = data[i] - data[j];
        Vector_3 vki = data[i] - data[k];
        Vector_3 vmi = data[i] - data[m];
        double lji = vji.squaredLength();
        double lki = vki.squaredLength();
        double lmi = vmi.squaredLength();
        double di = lji + lki + lmi;
        double Vi = cross(vki, vji)*vmi; 
        return di/std::pow(Vi, 2.0/3.0)/3.0;
    }

    
    double getUTQuality()
    {
        double sum_q = 0.0;
        for(int i = 0; i < 8; i++)
        {
            double q = getUTQuality(i);
            sum_q += q;
        }
        return sum_q/8.0;
    }

    double getMaxUTQuality()
    {
        double max_q = -1.0;
        for(int i = 0; i < 8; i++)
        {
            double q = getUTQuality(i);
            if(max_q < q)
                max_q = q;
        }
        return max_q;
    }

    double getUTQuality(int i)
    {
        int j = corner[i][0];
        int k = corner[i][1];
        int m = corner[i][2];
        Vector_3 vji = data[i] - data[j];
        Vector_3 vki = data[i] - data[k];
        Vector_3 vmi = data[i] - data[m];
        double lji = vji.squaredLength();
        double lki = vki.squaredLength();
        double lmi = vmi.squaredLength();
        double di = (lji*std::sqrt(lji) + lki*std::sqrt(lki) + lmi*std::sqrt(lmi));
        double Ji = dot(cross(vki, vji), vmi);
        double mui = 3*Ji/di;
        return (2 - mui)*(2 - mui);
    }


    void computeUTGradient(std::vector<Vector_3 > & grads, 
            std::vector<double> & ws, 
            double alpha,
            bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(8, Vector_3(0.0, 0.0, 0.0));
            ws.resize(8, 0.0);
        }

        for(int i = 0; i < 8; i++)
        {
            computeUTGradient(i, grads, ws, alpha);
        }
    }

    void computeUTGradient(int i, std::vector<Vector_3 > & grads,
            std::vector<double> & ws, double alpha)
    {
        int j = corner[i][0];
        int k = corner[i][1];
        int m = corner[i][2];
        Vector_3 vji = data[i] - data[j];
        Vector_3 vki = data[i] - data[k];
        Vector_3 vmi = data[i] - data[m];
        double lji = vji.squaredLength();
        double lki = vki.squaredLength();
        double lmi = vmi.squaredLength();
        double di = std::pow(lji, 1.5) + std::pow(lki, 1.5) + std::pow(lmi, 1.5);
        double Ji = dot(cross(vki, vji), vmi);
        double mui = 3*Ji/di;
        double q = (2-mui)*(2-mui)*(2-mui)*(2-mui);
        double ci = 6.0*alpha*std::pow(mui - 2.0, 2*alpha - 1); 

        lji = std::sqrt(lji);
        lki = std::sqrt(lki);
        lmi = std::sqrt(lmi);

        grads[i] += ci*((cross(vji, vmi) + cross(vmi, vki) + cross(vki, vji)) 
            - mui*(lji*vji + lki*vki + lmi*vmi))/di;
        ws[i] += std::abs(ci)*(1 + std::abs(mui))*(lji+lki+lmi)/di;

        grads[j] += ci*(cross(vki, vmi) + mui*lji*vji)/di;
        ws[j] += std::abs(ci)*(lki + std::abs(mui)*lji)/di;
        
        grads[k] += ci*(cross(vmi, vji) + mui*lki*vki)/di;
        ws[k] += std::abs(ci)*(lmi + std::abs(mui)*lki)/di;
        
        grads[m] += ci*(cross(vji, vki) + mui*lmi*vmi)/di;
        ws[m] += std::abs(ci)*(lji + std::abs(mui)*lmi)/di;
    }


    void computeTMOPGradient(std::vector<Vector_3 > & grads,  std::vector<double> & ws, bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(8, Vector_3(0.0, 0.0, 0.0));
            ws.resize(8, 0.0);
        }

        for(int i = 0; i < 8; i++)
        {
            double J = getCornerJacobi(i);
            if(J > EPS)
                computeTMOPGradient(i, grads, ws);
        }
    }

    void computeTMOPGradient(int i, std::vector<Vector_3 > & grads,  std::vector<double> & ws)
    {
        assert( i >= 0 && i < 8);
        int j = corner[i][0];
        int k = corner[i][1];
        int m = corner[i][2];
        Vector_3 vji = data[i] - data[j];
        Vector_3 vki = data[i] - data[k];
        Vector_3 vmi = data[i] - data[m];
        double lji = vji.squaredLength();
        double lki = vki.squaredLength();
        double lmi = vmi.squaredLength();
        double di = lji + lki + lmi;
        lji = std::sqrt(lji);
        lki = std::sqrt(lki);
        lmi = std::sqrt(lmi);
        double Vi = dot(cross(vki, vji), vmi); 
        double mui = di/std::pow(Vi, 2.0/3.0)/3.0;
        double w1 = 2.0*mui/di;
        double w2 = 2.0*mui/Vi/3.0;
        grads[i] += w1*(vji + vki + vmi) +w2*(cross(vmi, vji) + cross(vji, vki) + cross(vki, vmi));
        grads[j] += w2*cross(vmi, vki) - w1*vji;
        grads[k] += w2*cross(vji, vmi) - w1*vki;
        grads[m] += w2*cross(vki, vji) - w1*vmi;
        ws[i] += 3*w1 + w2*(lji + lki + lmi); 
        ws[j] += w1 + w2*lmi;
        ws[k] += w1 + w2*lji;
        ws[m] += w1 + w2*lki;
    }

    void computeMCGradient(std::vector<Vector_3 > & grads, std::vector<double> & ws, bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(8, Vector_3(0.0, 0.0, 0.0));
            ws.resize(8, 0.0);
        }

        for(int i = 0; i < 8; i++)
        {
            computeMCGradient(i, grads[i]);
            ws[i] = 1.0;
        }
        return;
    }

    void computeMCGradient(int i, Vector_3 & gradi)
    {
        int j = corner[i][0];
        int k = corner[i][1];
        int m = corner[i][2];
        computeMCGradient(i, j, m, gradi);
        computeMCGradient(i, m, k, gradi);
        computeMCGradient(i, k, j, gradi);
        gradi = gradi/3.0;
    }

    void computeMCGradient(int i, int j, int k,  Vector_3 & gradi)
    {
        double J = getCornerJacobi(i);
        Vector_3 vij = data[j] - data[i];
        Vector_3 vik = data[k] - data[i];
        double h = (std::sqrt(vij.squaredLength()) + std::sqrt(vik.squaredLength()))/2.0;
        Vector_3 va = cross(vij, vik);
        if(J < EPS) 
            gradi = gradi - va/h;    
        else
            gradi = gradi + va/h;
    }

    void computeLpGradient(std::vector<Vector_3 > & grads, std::vector<double> & ws, bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(8, Vector_3(0.0, 0.0, 0.0));
            ws.resize(8, 0.0);
        }

        for(int i = 0; i < 8; i++)
        {
            computeLpGradient(i, grads[i]);
            ws[i] = 1.0;
        }
        return;
    }

    void computeLpGradient(int i, Vector_3 & gradi)
    {
        int j = corner[i][0];
        int k = corner[i][1];
        int m = corner[i][2];
        Vector_3 vji = data[i] - data[j];
        Vector_3 vki = data[i] - data[k];
        Vector_3 vmi = data[i] - data[m];
        gradi += vji;
        gradi += vki;
        gradi += vmi;
        gradi = gradi/3.0;
    }
private:
    Point_3 data[8];
};


}

}

#endif
