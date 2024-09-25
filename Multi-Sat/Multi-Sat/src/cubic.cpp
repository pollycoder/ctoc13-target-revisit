/**
 * @brief     implementation of Cardano's method for solving cubic equations
 * @author    Thomas Atwood (tatwood.net)
 * @date      2011
 * @copyright unlicense / public domain
 ****************************************************************************/
#include <math.h>

//20200925改成了双精度浮点数，精度也提高到1e-12
const double EPSILON = 1.0e-12;

int solve_cubic(double a, double b, double c, double d, double* xout)
{
    static const double cos120 = -0.5;
    static const double sin120 = 0.86602540378443864676372317075294;
    int n = 0;
    if(fabs(d) < EPSILON)
    {
        // first solution is x = 0
        *xout = 0.0;
        ++n;
        ++xout;
        // divide all terms by x, converting to quadratic equation
        d = c;
        c = b;
        b = a;
        a = 0.0;
    }
    if(fabs(a) < EPSILON)
    {
        if(fabs(b) < EPSILON)
        {
            // linear equation
            if(fabs(c) > EPSILON)
            {
                *xout = -d/c;
                n += 1;
            }
        }
        else
        {
            // quadratic equation
            double yy = c*c - 4.0*b*d;
            if(yy >= 0)
            {
                double inv2b = 1.0/(2.0*b); 
                double y = sqrt(yy);
                xout[0] = (-c + y) * inv2b;
                xout[1] = (-c - y) * inv2b;
                n += 2;
            }
        }
    }
    else
    {
        // cubic equation
        double inva = 1.0/a;
        double invaa = inva*inva;
        double bb = b*b;
        double bover3a = b*(1/3.0f)*inva;
        double p = (3.0*a*c - bb)*(1.0/3.0)*invaa;
        double halfq = (2.0*bb*b - 9.0*a*b*c + 27.0*a*a*d)*(0.5/27.0)*invaa*inva;
        double yy = p*p*p/27.0 + halfq*halfq;
        if(yy > EPSILON)
        {
            // sqrt is positive: one real solution
            double y = sqrt(yy);
            double uuu = -halfq + y;
            double vvv = -halfq - y;
            double www = fabs(uuu) > fabs(vvv) ? uuu : vvv;
            double w = (www < 0.0) ? -pow(fabs(www),1.0/3.0) : pow(www,1.0/3.0);
            *xout = w - p/(3.0*w) - bover3a;
            n = 1;
        }
        else if(yy < -EPSILON)
        {
            // sqrt is negative: three real solutions
            double x = -halfq;
            double y = sqrt(-yy);
            double theta;
            double r;
            double ux;
            double uyi;
            // convert to polar form
            if(fabs(x) > EPSILON)
            {
                theta = (x > 0) ? atan(y/x) : (atan(y/x) + 3.14159625f);
                r = sqrt(x*x - yy);
            }
            else
            {
                // vertical line
                theta = 3.1415926535897932384626433832795/2;
                r = y;
            }
            // calc cube root
            theta /= 3.0;
            r = pow(r, 1.0/3.0);
            // convert to complex coordinate
            ux = cos(theta)*r;
            uyi = sin(theta)*r;
            // first solution
            xout[0] = ux+ux - bover3a;
            // second solution, rotate +120 degrees
            xout[1] = 2.0*(ux*cos120 - uyi*sin120) - bover3a;
            // third solution, rotate -120 degrees
            xout[2] = 2.0*(ux*cos120 + uyi*sin120) - bover3a;
            n = 3;
        }
        else
        {
            // sqrt is zero: two real solutions
            double www = -halfq;
            double w = (www < 0) ? -pow(fabs(www),1.0/3.0) : pow(www,1.0/3.0); 
            // first solution           
            xout[0] = w+w - bover3a;
            // second solution, rotate +120 degrees
            xout[1] = 2.0*w*cos120 - bover3a;
            n = 2;
        }
    }
    return n;
}