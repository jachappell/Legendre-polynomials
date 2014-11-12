//==================================================================
/**
 *  legendre.h -- C++ functions to evaluate Legendre polynomials
 *
 *  Copyright (C) 2005 by James A. Chappell
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or
 *  sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  condition:
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */
//=================================================================
/*
 * legendre.h:  Version 0.01
 * Created by James A. Chappell <rlrrlrll@gmail.com>
 * http://www.storage-b.com/math-numerical-analysis/18
 * Created 29 September 2005
 *
 * History:
 * 29-sep-2005  created
 */
//==============

#ifndef __LEGENDRE_H__
#define __LEGENDRE_H__
/*
 *	Function calculates Legendre Polynomials Pn(x)
 */

namespace Legendre
{
  // n = 0
  inline double P0(double x)
  {
    return 1.0 ;
  }

  // n = 1
  inline double P1(double x)
  {
    return x ;
  }

  // n = 2
  inline double P2(double x)
  {
    return ((3.0 * x*x) - 1.0) * 0.5 ;
  }

/*
 *	Pn(x)
 */
  inline double Pn(unsigned int n, double x)
  {
    if (n == 0)
    {
      return P0(x) ;
    }
    else if (n == 1)
    {
      return P1(x) ;
    }
    else if (n == 2)
    {
      return P2(x) ;
    }
    
    if (x == 1.0)
    {
      return 1.0 ;
    }

    if (x == -1.0)
    {
      return ((n % 2 == 0) ? 1.0 : -1.0) ;
    }

    if ((x == 0.0) && (n % 2))
    {
      return 0.0 ;
    }

/* We could simply do this:
    return (double(((2 * n) - 1)) * x * Pn(n - 1, x) -
          (double(n - 1)) * Pn(n - 2, x)) / (double)n ;
   but it could be slow for large n */
  
    double pnm1(P2(x)) ;
    double pnm2(P1(x)) ;
    double pn(pnm1) ;

    for (unsigned int l = 3 ; l <= n ; l++)
    { 
      pn = (((2.0 * (double)l) - 1.0) * x * pnm1 - 
            (((double)l - 1.0) * pnm2)) / (double)l ;
      pnm2 = pnm1;
      pnm1 = pn;
    }

    return pn ;
  }
}
#endif
