#include <cmath>
#include <array>
#include "xmath.h"

double PI = 4.0 * atan(1.0);

double Vector3 :: mag() const
{
    return sqrt(x * x + y * y + z * z);
}

double Vector3 :: mag2() const
{
    return x * x + y * y + z * z;
}

Vector3 Vector3 :: normalize()
{
    double mag = this->mag();
    return Vector3(x / mag, y / mag, z / mag);
}

double Vector3 :: dot(const Vector3& other) const
{
    return x * other.x + y * other.y + z * other.z;
}

Vector3 Vector3 :: cross(const Vector3& other) const
{
    return Vector3(y * other.z - z * other.y,
                    z * other.x - x * other.z,
                    x * other.y - y * other.x);
}

double& Matrix3 :: operator()(int i, int j)
{
    return arr[3 * i + j];
}

double Matrix3 :: operator()(int i, int j) const
{
    return arr[3 * i + j];
}

double& Matrix3 :: operator()(int i) 
{
    return arr[i]; 
}

double Matrix3 :: operator()(int i) const 
{ 
    return arr[i]; 
}

double Matrix3 :: determinant() const
{
    return (*this)(0, 0) * (*this)(1, 1) * (*this)(2, 2) +
            (*this)(0, 1) * (*this)(1, 2) * (*this)(2, 0) +
            (*this)(0, 2) * (*this)(1, 0) * (*this)(2, 1) -
            (*this)(0, 2) * (*this)(1, 1) * (*this)(2, 0) -
            (*this)(0, 0) * (*this)(1, 2) * (*this)(2, 1) -
            (*this)(0, 1) * (*this)(1, 0) * (*this)(2, 2);
}

Matrix3 Matrix3 :: transpose() const
{
    Matrix3 result;
    result(0, 1) = (*this)(1, 0);
    result(0, 2) = (*this)(2, 0);
    result(1, 0) = (*this)(0, 1);
    result(1, 2) = (*this)(2, 1);
    result(2, 0) = (*this)(0, 2);
    result(2, 1) = (*this)(1, 2);
    return result;
}

Matrix3 Matrix3 :: inverse() const
{
    Matrix3 result;
    float det = determinant();
    result(0, 0) = ((*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1)) / det;
    result(0, 1) = ((*this)(0, 2) * (*this)(2, 1) - (*this)(0, 1) * (*this)(2, 2)) / det;
    result(0, 2) = ((*this)(0, 1) * (*this)(1, 2) - (*this)(0, 2) * (*this)(1, 1)) / det;
    result(1, 0) = ((*this)(1, 2) * (*this)(2, 0) - (*this)(1, 0) * (*this)(2, 2)) / det;
    result(1, 1) = ((*this)(0, 0) * (*this)(2, 2) - (*this)(0, 2) * (*this)(2, 0)) / det;
    result(1, 2) = ((*this)(0, 2) * (*this)(1, 0) - (*this)(0, 0) * (*this)(1, 2)) / det;
    result(2, 0) = ((*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0)) / det;
    result(2, 1) = ((*this)(0, 1) * (*this)(2, 0) - (*this)(0, 0) * (*this)(2, 1)) / det;
    result(2, 2) = ((*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0)) / det;

    return result;
}

Quaternion Quaternion :: operator~() const
{
    return Quaternion(w, -x, -y, -z);
}

inline double Quaternion :: mag() const
{
    return sqrt(w * w + x * x + y * y + z * z);
}

inline Quaternion Quaternion :: normalize() const
{
    double mag = this->mag();
    return Quaternion(w / mag, x / mag, y / mag, z / mag);
}

inline Matrix3 Quaternion :: as_matrix() const
{
    Matrix3 result;
    result(0, 0) = 1 - 2 * y * y - 2 * z * z;
    result(0, 1) = 2 * x * y + 2 * w * z;
    result(0, 2) = 2 * x * z - 2 * w * y;
    result(1, 0) = 2 * x * y - 2 * w * z;
    result(1, 1) = 1 - 2 * x * x - 2 * z * z;
    result(1, 2) = 2 * y * z + 2 * w * x;
    result(2, 0) = 2 * x * z + 2 * w * y;
    result(2, 1) = 2 * y * z - 2 * w * x;
    result(2, 2) = 1 - 2 * x * x - 2 * y * y;
    return result;
}

Vector3 operator+(const Vector3& a, const Vector3& b) 
{
    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Matrix3 operator+(const Matrix3& a, const Matrix3& b)
{
    Matrix3 result;
    result(0, 0) = a(0, 0) + b(0, 0);
    result(0, 1) = a(0, 1) + b(0, 1);
    result(0, 2) = a(0, 2) + b(0, 2);
    result(1, 0) = a(1, 0) + b(1, 0);
    result(1, 1) = a(1, 1) + b(1, 1);
    result(1, 2) = a(1, 2) + b(1, 2);
    result(2, 0) = a(2, 0) + b(2, 0);
    result(2, 1) = a(2, 1) + b(2, 1);
    result(2, 2) = a(2, 2) + b(2, 2);
    return result;
}

Vector3 operator-(const Vector3& a, const Vector3& b) 
{
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Matrix3 operator-(const Matrix3& a, const Matrix3& b)
{
    Matrix3 result;
    result(0, 0) = a(0, 0) - b(0, 0);
    result(0, 1) = a(0, 1) - b(0, 1);
    result(0, 2) = a(0, 2) - b(0, 2);
    result(1, 0) = a(1, 0) - b(1, 0);
    result(1, 1) = a(1, 1) - b(1, 1);
    result(1, 2) = a(1, 2) - b(1, 2);
    result(2, 0) = a(2, 0) - b(2, 0);
    result(2, 1) = a(2, 1) - b(2, 1);
    result(2, 2) = a(2, 2) - b(2, 2);
    return result;
}

Vector3 operator*(const double scalar, const Vector3& other)
{
    return Vector3(scalar * other.x, scalar * other.y, scalar * other.z);
}

Vector3 operator*(const Vector3& other, const double scalar)
{
    return Vector3(scalar * other.x, scalar * other.y, scalar * other.z);
}

Matrix3 operator*(const double scalar, const Matrix3& other)
{
    Matrix3 result;
    for (int i = 0; i < 9; i++)
    {
        result(i) = scalar * other(i);
    }
    return result;
}

Matrix3 operator*(const Matrix3& other, const double scalar)
{
    Matrix3 result;
    for (int i = 0; i < 9; i++)
    {
        result(i) = scalar * other(i);
    }
    return result;
}

Quaternion operator*(const double scalar, const Quaternion& other)
{
    return Quaternion(scalar * other.w, scalar * other.x, scalar * other.y, scalar * other.z);
}

Quaternion operator*(const Quaternion& other, const double scalar)
{
    return Quaternion(scalar * other.w, scalar * other.x, scalar * other.y, scalar * other.z);
}

Vector3 operator/(const Vector3& other, const double scalar)
{
    return Vector3(other.x / scalar, other.y / scalar, other.z / scalar);
}

Matrix3 operator/(const Matrix3& other, const double scalar)
{
    Matrix3 result;
    for (int i = 0; i < 9; i++)
    {
        result(i) = other(i) / scalar;
    }
    return result;
}

Quaternion operator/(const Quaternion& other, const double scalar)
{
    return Quaternion(other.w / scalar, other.x / scalar, other.y / scalar, other.z / scalar);
}

Vector3 operator*(const Vector3& A, const Vector3& B)    // VECTOR COMPONENT MULTIPLICATION
{
    return Vector3(A.x * B.x, A.y * B.y, A.z * B.z);
}

Vector3 operator*(const Vector3& A, const Matrix3& B)    // VECTOR MATRIX MULTIPLICATION
{
    Vector3 result;
    
    result.x = A.x * B(0, 0) + A.y * B(1, 0) + A.z * B(2, 0);
    result.y = A.x * B(0, 1) + A.y * B(1, 1) + A.z * B(2, 1);
    result.z = A.x * B(0, 2) + A.y * B(1, 2) + A.z * B(2, 2);

    return result;
}

Vector3 operator*(const Matrix3& A, const Vector3& B)    // MATRIX VECTOR MULTIPLICATION
{
    Vector3 result;

    result.x = A(0, 0) * B.x + A(0, 1) * B.y + A(0, 2) * B.z;
    result.y = A(1, 0) * B.x + A(1, 1) * B.y + A(1, 2) * B.z;
    result.z = A(2, 0) * B.x + A(2, 1) * B.y + A(2, 2) * B.z;

    return result;
}

Matrix3 operator*(const Matrix3& A, const Matrix3& B)    // MATRIX MATRIX MULTIPLICATION
{
    Matrix3 result;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result(i, j) = A(i, 0) * B(0, j) + A(i, 1) * B(1, j) + A(i, 2) * B(2, j);
        }
    }
    return result;
}

Quaternion operator*(const Quaternion& A, const Quaternion& B)    // QUATERNION MULTIPLICATION
{
    return Quaternion(A.w * B.w - A.x * B.x - A.y * B.y - A.z * B.z,
                      A.w * B.x + A.x * B.w + A.y * B.z - A.z * B.y,
                      A.w * B.y - A.x * B.z + A.y * B.w + A.z * B.x,
                      A.w * B.z + A.x * B.y - A.y * B.x + A.z * B.w);
}

double degrees_to_radians(double degrees) {
    return degrees * PI / 180.0;
}

double radians_to_degrees(double radians) {
    return radians * 180.0 / PI;
}

double dot_product ( const Vector3& a, const Vector3& b )
{
    return a.dot(b);
}

Vector3 cross_product ( const Vector3& a, const Vector3& b )
{
    Vector3 c = a.cross(b);
    return c;
}