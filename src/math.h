#ifndef SPHSOL_MATH_H_
#define SPHSOL_MATH_H_

#include"global.h"


//--------------------------------------------------------------------------------------------------
// Function power
//--------------------------------------------------------------------------------------------------
inline double pow2(double val) { return val * val; }
inline double pow3(double val) { return val * val * val; }
inline double pow4(double val) { return val * val * val * val; }


//--------------------------------------------------------------------------------------------------
// Funtion min / max
//--------------------------------------------------------------------------------------------------
template <typename T>
inline T min(const T& x, const T& y) { return (x < y) ? x : y; }
template <typename T>
inline T max(const T& x, const T& y) { return (x > y) ? x : y; }
template <typename T, typename... Ts>
inline T min(const T& x, const T& y, const Ts&... xs) { return min(min(x, y), xs...); }
template <typename T, typename... Ts>
inline T max(const T& x, const T& y, const Ts&... xs) { return max(max(x, y), xs...); }


//--------------------------------------------------------------------------------------------------
// Search a value's position in a given list using bisection method
//--------------------------------------------------------------------------------------------------
int BinarySearch(const std::vector<double>& list, double val);


//--------------------------------------------------------------------------------------------------
// Judge if a point is ouside the volume limit
//--------------------------------------------------------------------------------------------------
inline bool isOutsideLimit(const Vec3d& pos, const Vec3d& vol_min, const Vec3d& vol_max)
{
  if (pos(0) < vol_min(0) || pos(0) > vol_max(0) || pos(1) < vol_min(1) ||
      pos(1) > vol_max(1) || pos(2) < vol_min(2) || pos(2) > vol_max(2))
  {
    return true;
  }
  else
  {
    return false;
  }
};


//--------------------------------------------------------------------------------------------------
// Poisson distribution: Knuth's algorithm is recomended for lambda < 30
//--------------------------------------------------------------------------------------------------
int PoissonDistribution(double lambda);


//--------------------------------------------------------------------------------------------------
// Function to return median of three numbers
//--------------------------------------------------------------------------------------------------
template <typename T>
T middleOfThree(T a, T b, T c)
{
  // Checking for b
  if ((a < b && b < c) || (c < b && b < a))
    return b;
  // Checking for a
  else if ((b < a && a < c) || (c < a && a < b))
    return a;
  else
    return c;
}


//--------------------------------------------------------------------------------------------------
// method for calculating the pseudo-Inverse as recommended by Eigen developers
//--------------------------------------------------------------------------------------------------
template<typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
{
	Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
  // For a non-square matrix
  // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}


#endif
