//==================================================================
/**
 *  legendre.h -- C++ functions to evaluate Legendre polynomials
 *
 *  Copyright (C) 2019 by James A. Chappell
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
 * legendre.h:  Version 0.03
 * Created by James A. Chappell <rlrrlrll@gmail.com>
 * http://www.storage-b.com/math-numerical-analysis/18
 * Created 29 September 2005
 *
 * History:
 * 29-sep-2005  created
 * 14-nov-2014  templates
 * 01-nov-2019  deduced types
 */
//==============

#ifndef __LEGENDRE_H__
#define __LEGENDRE_H__
/*
 *	Function calculates Legendre Polynomials Pn(x)
 */
namespace Storage_B
{
  namespace Legendre
  {
    // n = 0
    template <class T> inline auto P0(const T& x)
    {
      return static_cast<T>(1);
    }

    // n = 1
    template <class T> inline auto P1(const T& x)
    {
      return x;
    }

    // n = 2
    template <class T> inline auto P2(const T& x)
    {
      return ((static_cast<T>(3) * x*x) - static_cast<T>(1)) /
        static_cast<T>(2);
    }

/*
 *	  Pn(x)
 */
    template <class T> inline auto Pn(unsigned int n, const T& x)
    {
      switch(n)
      {
        case 0:
          return P0<T>(x);

        case 1:
          return P1<T>(x);

        case 2:
          return P2<T>(x);

        default:
          break;
      }
      
/*  We could simply do this:
      return (static_cast<T>(((2 * n) - 1)) * x * Pn(n - 1, x) -
          (static_cast<T>(n - 1)) * Pn(n - 2, x)) / static_cast<T>(n);
    but it could be slow for large n */

      auto pnm1(P2<T>(x));
      auto pnm2(P1<T>(x));
      T pn;

      for (auto m = 3u ; m <= n ; ++m)
      { 
        pn = (((static_cast<T>(2) * static_cast<T>(m)) - static_cast<T>(1))
            * x * pnm1 - (static_cast<T>(m - 1) * pnm2))
              / static_cast<T>(m);
        pnm2 = pnm1;
        pnm1 = pn;
      }

      return pn;
    }
  }
}
#endif
