#include <iostream>
#include "legendre.h"

using namespace std ;
using namespace Legendre ;

int main()
{
  double pn ;

  cout.precision(5) ;
  for (unsigned int n = 0 ; n <= 5 ; n++)
  {
    for (double x = -1.0 ; x <= 1.0 ; x = x + 0.1)
    { 
      pn = Pn<double>(n, x) ;
      cout << "P" << n << "(" << x << ") = " << pn << endl ;
    }
    cout << endl ;
  }

  return 0 ;
}
