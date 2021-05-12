#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "Point_3.h"
#include "Vector_3.h"
#include "constant.h"

#include <cmath>
#include <vector>

namespace  iMath {

namespace GeometryObject {

class Tetrahedron
{
    static int face[4][3];
    static int index[12][4];
    static int edge_idx[4][4];
    static int edge[6][2];

public:

    Tetrahedron()
    {
        data[0] = Point_3(0.0, 0.0, 0.0);
        data[1] = Point_3(1.0, 0.0, 0.0);
        data[2] = Point_3(0.0, 1.0, 0.0);
        data[3] = Point_3(0.0, 0.0, 1.0);
    }

    Point_3 & operator[] (const int i)
    {
        return data[i];
    }

    const Point_3 & operator[] (const int i) const
    {
        return data[i];
    }

    double quality()
    {
        return radiusRatio();
    }

    double radiusRatio()
    {
       double sa = surfaceArea();
       double vol = volume();
       Vector_3 d = direction(0);
       double di = std::sqrt(d.squaredLength());
       double R = di/vol/12.0;
       double r = 3.0*vol/sa;
       return R/r/3.0;
    }

    bool isValid()
    {
        double vol = volume();
        return vol > EPS;
    }

    double getUTQuality(double alpha = 2)
    {
        double sn = 0.0;
        double ln = 0.0;

        for(int i = 0; i < 4; i++)
        {
             sn +=faceArea(i);
        }

        for(int i = 0; i < 6; i++)
        {
            ln += sideLength(i);
        }

        double J = 6.0*volume();
        double d = 6.0*std::sqrt(6.0)/sn/ln; 
        return  std::pow(d*J - 2, 2*alpha); 
    }

    Vector_3 direction(int i)
    {
        Vector_3 v10 = data[index[3*i][0]] - data[index[3*i][1]];
        Vector_3 v20 = data[index[3*i][0]] - data[index[3*i][2]];
        Vector_3 v30 = data[index[3*i][0]] - data[index[3*i][3]];
        double L1 = v10.squaredLength();
        double L2 = v20.squaredLength();
        double L3 = v30.squaredLength();
        Vector_3 d =  L1*cross(v20,v30) + L2*cross(v30, v10) + L3*cross(v10, v20);
        return d;
    }

    double volume()
    {
        Vector_3 v01 = data[1] - data[0];
        Vector_3 v02 = data[2] - data[0];
        Vector_3 v03 = data[3] - data[0];
        return  dot(cross(v01, v02),v03)/6.0;
    }

    double faceArea(int i)
    {
         Vector_3 v01 = data[face[i][1]] - data[face[i][0]];
         Vector_3 v02 = data[face[i][2]] - data[face[i][0]];
         return std::sqrt(cross(v01, v02).squaredLength())/2.0;
    }

    double surfaceArea()
    {
        double sum = 0.0;
        for(int i = 0; i < 4; i++)
        {
            sum += faceArea(i);
        }
        return sum;
    }

    double sideLength(int i)
    {
        return std::sqrt((data[edge[i][1]] - data[edge[i][0]]).squaredLength() );
    }

    template<class Mat, class Vec>
    void buildLocalLinearSystem(Mat & A, 
            Vec & rhs_x, Vec & rhs_y, Vec & rhs_z)
    {
        //Vector_3 d[4];
        double s[4];
        double sn = 0.0;
        for(int i = 0; i < 4; i++)
        {
         //   d[i] = direction(i);
            s[i] = faceArea(i);
            sn += s[i];
        }

        double l[6];
        double ln = 0.0;
        for(int i = 0; i < 6; i++)
        {
            l[i] = sideLength(i);
            ln += l[i];
        }
        
        //double dd = d[0].squaredLength();
        //double R =std::sqrt(dd)/vol/12.0;
        //double r = 3.0*vol/sum;
        //double mu = R/r/3.0;

        double vol = volume();
        double mu = ln*sn/vol/36.0/std::sqrt(6.0);
        std::vector<Vector_3> rhs(4, Vector_3(0.0, 0.0, 0.0));
        for(int idx =0; idx < 12; idx ++)
        {
            int i = index[idx][0];  int j = index[idx][1];
            int k = index[idx][2]; int m = index[idx][3];
            A[i][j] =  -0.25*mu*(dot((data[k] - data[i]), (data[k] - data[j]))/s[m] +
                                dot(data[m] - data[i], data[m] - data[j])/s[k])/sn;
            A[i][j] -= mu/ln/l[edge_idx[i][j]];
            A[i][i] += std::abs(A[i][j]);
            rhs[i] += cross((data[k] - data[m]) + (data[j] - data[m]), data[i] - data[j])/vol/18.0;


            //A[i][j] = 2.0*dot(d[i], cross(data[i] - data[k], data[i] - data[m]))/dd;
            //A[i][j] = -0.25*mu*(dot(data[i] - data[m], data[j] - data[m])/s[k] +
            //                    dot(data[i] - data[k], data[j] - data[k])/s[m])/sum;
            //A[i][i] += std::abs(A[i][j]);
            
            //double t = (data[i] - data[m]).squaredLength() - (data[i] - data[k]).squaredLength();
            //rhs[i] += t*cross(d[i], data[i] - data[j])/dd; 
            //rhs[i] += cross((data[k] - data[m]) + (data[j] - data[m]), data[i] - data[j])/(9*vol);
            //rhs[i]+=2.0/dd*dot(d[i],cross(data[i]-data[k], data[i]-data[m]))*(data[i]-data[j]);
        }

        for( int i = 0; i < 4; i++)
        {
            rhs[i] = -mu*rhs[i];
            rhs_x[i] = rhs[i][0];
            rhs_y[i] = rhs[i][1];
            rhs_z[i] = rhs[i][2];
        }
    }

    
    Point_3 & vertex(int i)
    {
        assert( i >=0 && i < 4);
        return data[i];
    }

    Point_3 barycenter()
    {
        double x = data[0][0] + data[1][0] + data[2][0] + data[3][0];
        double y = data[0][1] + data[1][1] + data[2][1] + data[3][1];
        double z = data[0][2] + data[1][2] + data[2][2] + data[3][2];
        return Point_3(x/4.0, y/4.0, z/4.0);
    }

    void computeCPTGradient(std::vector<Vector_3> & grads, 
            std::vector<double> & weight,
            bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(4, Vector_3(0.0, 0.0, 0.0));
            weight.resize(4, 0.0);
        }
        
        Point_3 bp = barycenter();  
        for(int i = 0; i < 4; i++)
        {
            weight[i] = 1.0;
            grads[i] = data[i] - bp;
        }

    }


    void computeUTGradient(std::vector<Vector_3> & grads,
            std::vector<double> & weight,
            double alpha = 2,
            bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(4, Vector_3(0.0, 0.0, 0.0));
            weight.resize(4, 0.0);
        }

        double s[4];
        double l[6];
        double sn = 0.0;
        double ln = 0.0;

        for(int i = 0; i < 4; i++)
        {
            s[i] = faceArea(i);
            sn += s[i];
        }

        for(int i = 0; i < 6; i++)
        {
            l[i] = sideLength(i);
            ln += l[i];
        }

        double J = 6*volume();
        double d = 6.0*std::sqrt(6.0)/sn/ln; 
        double q = d*J; 
        double c = 2*alpha*std::pow(q - 2, 2*alpha - 1);

        for(int idx = 0; idx < 12; idx += 3)
        {
            int i = index[idx][0];  int j = index[idx][1];
            int k = index[idx][2]; int m = index[idx][3];
            Vector_3 vji = data[i] - data[j];
            Vector_3 vki = data[i] - data[k];
            Vector_3 vmi = data[i] - data[m];

            grads[i] += d*(cross(vji, vmi) 
                + cross(vmi, vki)
                + cross(vki, vji));

            double wji = 0.25*(dot(data[i] - data[m], data[j] - data[m])/s[k] +
                          dot(data[i] - data[k], data[j] - data[k])/s[m])/sn;
            double wki = 0.25*(dot(data[i] - data[m], data[k] - data[m])/s[j] +
                          dot(data[i] - data[j], data[k] - data[j])/s[m])/sn;
            double wmi = 0.25*(dot(data[i] - data[j], data[m] - data[j])/s[k] +
                          dot(data[i] - data[k], data[m] - data[k])/s[j])/sn;

            grads[i] -= q/sn*(wji*vji + wki*vki + wmi*vmi);
            grads[i] -= q/ln*(vji/l[edge_idx[j][i]] + vki/l[edge_idx[k][i]] + vmi/l[edge_idx[m][i]]);
            grads[i] *= c;

            weight[i] += d*(l[edge_idx[j][i]] + l[edge_idx[k][i]] + l[edge_idx[m][i]]);
            weight[i] += std::abs(q)*(wji + wki + wmi);
            weight[i] += std::abs(q)*(1.0/l[edge_idx[j][i]] + 1.0/l[edge_idx[k][i]] + 1.0/l[edge_idx[m][i]])/ln;
            weight[i] *= std::abs(c);
        }
    }

    void computeRRGradient(std::vector<Vector_3 > & grads, std::vector<double> & weight, bool init_zero = true)
    {
        if(init_zero)
        {
            grads.resize(4, Vector_3(0.0, 0.0, 0.0));
            weight.resize(4, 0.0);
        }

        Vector_3 d[4];
        double s[4];
        double sum= 0.0;
        for(int i = 0; i < 4; i++)
        {
            d[i] = direction(i);
            s[i] = faceArea(i);
            sum += s[i];
        }
        double dd = d[0].squaredLength();
        double vol = volume();

        double R = std::sqrt(dd)/vol/12.0;
        double r = 3.0*vol/sum;
        double Rr = R/r/3.0;

        for(int idx = 0; idx < 12; idx ++)
        {
            int i = index[idx][0];  int j = index[idx][1];
            int k = index[idx][2]; int m = index[idx][3];
            Vector_3 vji = data[i] - data[j];

            double w0 = 2.0*dot(cross(data[i] - data[k], data[i] - data[m]), d[i])/dd;
            double w1 = 0.25*(dot(data[i] - data[m], data[j] - data[m])/s[k] +
                          dot(data[i] - data[k], data[j] - data[k])/s[m])/sum;
            grads[i] += (w0 + w1)*vji;

            weight[i] += w0 + w1;

            double w2 = ((data[i] - data[m]).squaredLength() - (data[i] - data[k]).squaredLength())/dd;
            grads[i] += w2*cross(d[i], vji);

            grads[i] += cross((data[k] - data[m]) + (data[j] - data[m]), vji)/9.0/vol;
        }

        for(int i = 0; i < 4; i++)
        {
            grads[i] = Rr*grads[i];
            weight[i] = Rr*weight[i];
        }
    }

    void computeGradient(int idx, Vector_3 & grad, double & weight)
    {
        Vector_3 d = direction(idx);
        double s[4];
        double sum= 0.0;
        for(int i = 0; i < 4; i++)
        {
            s[i] = faceArea(i);
            sum += s[i];
        }
        double dd = d.squaredLength();
        double vol = volume();

        double R = std::sqrt(dd)/vol/12.0;
        double r = 3.0*vol/sum;
        double Rr =  R/r/3.0;

        for(int w = 3*idx; w < 3*(idx+1); w++ )
        {
            int i = index[w][0];  int j = index[w][1];
            int k = index[w][2]; int m = index[w][3];
            Vector_3 vji = data[i] - data[j];

            double w0 = 2.0*dot(cross(data[i] - data[k], data[i] - data[m]), d)/dd;
            double w1 = 0.25*(dot(data[i] - data[m], data[j] - data[m])/s[k] +
                          dot(data[i] - data[k], data[j] - data[k])/s[m])/sum;
            grad += (w0 + w1)*vji;

            weight += w0 + w1;

            double w2 = ((data[i] - data[m]).squaredLength() - (data[i] - data[k]).squaredLength())/dd;
            grad += w2*cross(d, vji);

            grad += cross((data[k] - data[m]) + (data[j] - data[m]), vji)/9.0/vol;
        }

        grad = Rr*grad;
        weight = Rr*weight;
    }

    void computeAngle(double &max_angle, double &min_angle)
    {

        std::vector<Vector_3 > ns(4);
        for(int i = 0; i < 4; i++)
        {
                ns[i] = cross(data[face[i][1]] - data[face[i][0]],
                        data[face[i][2]] - data[face[i][0]]);
                ns[i] /= std::sqrt(ns[i].squaredLength());
        }

        for(int i = 0; i < 4; i++)
            for(int j = i+1; j < 4; j++)
            {
                double angle = PI - std::acos(ns[i]*ns[j]);
                angle = angle/PI*180;
                if(max_angle < angle)
                {
                    max_angle = angle;
                }
                if(min_angle > angle)
                {
                    min_angle = angle;
                }
            }
    }
private:
    Point_3 data[4];
};

}

}

#endif // end of TETRAHEDRON_H
