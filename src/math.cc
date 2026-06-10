#include"math.h"


//--------------------------------------------------------------------------------------------------
// Search a value's position in a given list using bisection method
//--------------------------------------------------------------------------------------------------
int BinarySearch(const std::vector<double>& list, double val)
{
  int left, right, mid;

  left = 0;
  right = list.size() - 1;

  if (list[left] > val)
  {
    return left;
  }
  if (list[right] < val)
  {
    return right;
  }

  while (right - left > 1)
  {
    mid = (left + right) / 2;
    if (list[mid] >= val)
    {
      right = mid;
    }
    else
    {
      left = mid;
    }
  }

  return left;
};


//--------------------------------------------------------------------------------------------------
// Poisson distribution: Knuth's algorithm is recomended for lambda < 30
//--------------------------------------------------------------------------------------------------
int PoissonDistribution(double lambda)
{
  double L = exp(-lambda);
  int k = 0;
  double p = 1.0;
  do {
    k = k + 1;
    double u = static_cast<double>(rand()) / RAND_MAX;
    p = p * u;
  } while (p > L);
  
  return k - 1;
 }
 
