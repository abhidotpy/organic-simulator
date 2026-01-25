#pragma once

#include <cmath>
#include <array>

extern double PI;

// FIRST WE DEFINE THE NECESSARY DATA STRUCTURES ESSENTIAL FOR MOLECULAR SIMULATION

struct Vector3
{
    double x, y, z;
    Vector3() : x(0), y(0), z(0) {}
    Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

    double mag() const;
    double mag2() const;
    Vector3 normalize();

    double dot(const Vector3& other) const;
    inline Vector3 cross(const Vector3& other) const;
};

struct Matrix3
{
    std::array<double, 9> arr;
    Matrix3() : arr() {}

    double& operator()(int i, int j);
    double operator()(int i, int j) const;
    double& operator()(int i);
    double operator()(int i) const;

    double determinant() const;
    Matrix3 transpose() const;
    Matrix3 inverse() const;
};

struct Quaternion

{
    double w, x, y, z;
    Quaternion() : w(1), x(0), y(0), z(0) {}
    Quaternion(double w, double x, double y, double z) : w(w), x(x), y(y), z(z) {}
    Quaternion(Vector3 v) : w(0), x(v.x), y(v.y), z(v.z) {}

    Quaternion operator~() const;
    double mag() const;
    Quaternion normalize() const;
    Matrix3 as_matrix() const;
};


// DEFINE ADDITION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vector3 operator+(const Vector3& a, const Vector3& b);
Matrix3 operator+(const Matrix3& a, const Matrix3& b);

// DEFINE SUBTRACTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vector3 operator-(const Vector3& a, const Vector3& b);
Matrix3 operator-(const Matrix3& a, const Matrix3& b);

// DEFINE SCALAR MULTIPLICATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vector3 operator*(const double scalar, const Vector3& other);
Vector3 operator*(const Vector3& other, const double scalar);
Matrix3 operator*(const double scalar, const Matrix3& other);
Matrix3 operator*(const Matrix3& other, const double scalar);
Quaternion operator*(const double scalar, const Quaternion& other);
Quaternion operator*(const Quaternion& other, const double scalar);

// DEFINE SCALAR DIVISION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vector3 operator/(const Vector3& other, const double scalar);
Matrix3 operator/(const Matrix3& other, const double scalar);
Quaternion operator/(const Quaternion& other, const double scalar);


// DEFINE INTER OPERATOR MULTIPLICATION RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vector3 operator*(const Vector3& A, const Vector3& B);              // VECTOR COMPONENT MULTIPLICATION
Vector3 operator*(const Vector3& A, const Matrix3& B);              // VECTOR MATRIX MULTIPLICATION
Vector3 operator*(const Matrix3& A, const Vector3& B);              // MATRIX VECTOR MULTIPLICATION
Matrix3 operator*(const Matrix3& A, const Matrix3& B);              // MATRIX MATRIX MULTIPLICATION
Quaternion operator*(const Quaternion& A, const Quaternion& B);     // QUATERNION MULTIPLICATION

// DEFINE UTILITY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double degrees_to_radians(double degrees);
double radians_to_degrees(double radians);

double dot_product ( const Vector3& a, const Vector3& b );
Vector3 cross_product ( const Vector3& a, const Vector3& b );


template <typename Function>
inline double golden_section_search(Function func)
{
    const int max_iter = 10000;
    float eps = 1e-8, phi = 0.5 * ( 3.0 - sqrt(5.0) );
    float a, b, x1, x2, fx1, fx2;
    int I, count;

    a = 0.0; 
    b = 1.0;
    x1 = (1 - phi) * a + phi * b;
    x2 = phi * a + (1 - phi) * b;
    fx1 = func(x1); 
    fx2 = func(x2);
    count = 0;

    
    while ( std::abs(x2 - x1) > eps && count < max_iter )
    {
        count = count + 1;

        if ( fx1 < fx2 )
        {
            b = x2;
            x2 = x1;  
            fx2 = fx1;
            x1 = (1 - phi) * a + phi * b;
            fx1 = func(x1);
        }
        
        else
        {
            a = x1;
            x1 = x2; 
            fx1 = fx2;
            x2 = phi * a + (1 - phi) * b;
            fx2 = func(x2);
        }
    }

    return (x1 + x2) / 2.0;
}

template <typename Function>
inline double brent(Function func)
{
    const int max_iter = 10000;
    float eps = 1e-8, phi = 0.5 * ( 3.0 - sqrt(5.0) );
    float a, b, c, x, w, v, u;
    float fa, fb, fc, fx, fw, fv, fu;
    float deltax = 0.0, atol = 1e-5, rtol = 1e-3;
    float tol1, tol2, xmid, tmp1, tmp2, rat, p, dx_temp;
    float I, iter = 0;

    a = 0.0; b = 1.0;
    c = (1 - phi) * a + phi * b;
    fa = func(a); fb = func(b); fc = func(c);
    v = c; w = v; x = w;
    fv = fc; fw = fv; fx = fw;

    while (iter < max_iter)
    {
        tol1 = rtol * std::abs(x) + atol;
        tol2 = 2.0 * tol1;
        xmid = 0.5 * (a + b);

        if (std::abs(x - xmid) < (tol2 - 0.5 * (b - a))) 
        {
            break; 
        }

        if (std::abs(deltax) <= tol1)
        {
            if (x >= xmid)
            {
                deltax = a - x;
            }
            else
            {
                deltax = b - x;

            }
            rat = phi * deltax;
        }
        else
        {
            tmp1 = (x - w) * (fx - fv);
            tmp2 = (x - v) * (fx - fw);
            p = (x - v) * tmp2 - (x - w) * tmp1;
            tmp2 = 2.0 * (tmp2 - tmp1);

            if (tmp2 > 0.0) {p = -p;}

            tmp2 = std::abs(tmp2);
            dx_temp = deltax;
            deltax = rat;
            
            if ((p > tmp2 * (a - x)) && (p < tmp2 * (b - x)) && (std::abs(p) < std::abs(0.5 * tmp2 * dx_temp)))
            {
                rat = p * 1.0 / tmp2;
                u = x + rat;
                if ((u - a) < tol2 || (b - u) < tol2)
                {
                    if (xmid - x >= 0)
                    {
                        rat = tol1;
                    }
                    else
                    {
                        rat = -tol1;
                    }
                }
            }
            else
            {
                if (x >= xmid)
                {
                    deltax = a - x;
                }
                else
                {
                    deltax = b - x;
                }
                rat = phi * deltax;
            }
        }

        if (std::abs(rat) < tol1)
        {
            if (rat >= 0)
            {
                u = x + tol1;
            }
            else
            {
                u = x - tol1;
            }
        }
        else
        {
            u = x + rat;
        }

        fu = func(u);

        if (fu > fx)
        {
            if (u < x)
            {
                a = u;
            }
            else
            {
                b = u;
            }

            if ((fu <= fw) || (w == x))
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }

            else if ((fu <= fv) || (v == x) || (v == w))
            {
                v = u;
                fv = fu;
            }
        }
        else
        {
            if (u >= x)
            {
                a = x;
            }
            else
            {
                b = x;
            }

            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        iter = iter + 1;
    }

    return x;
}