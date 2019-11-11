#include "legendre.h"

#include "Random.h"

#include <memory>

#define BOOST_TEST_MODULE Legendre
#include <boost/test/included/unit_test.hpp>

using namespace Storage_B;

namespace
{
  Random<double> dRan(-1.0, 1.0);
  Random<float> fRan(-1.0f, 1.0f);
  Random<unsigned int> iRan(2, 10);

  auto tol = 0.00049;
}

// P0
BOOST_AUTO_TEST_CASE(P0)
{
  auto val1 = Legendre::P0<double>(dRan());
  auto val2 = Legendre::Pn<double>(0, dRan());

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);
}

// P1
BOOST_AUTO_TEST_CASE(P1)
{
  auto x = fRan();

  auto val1 = Legendre::P1<float>(x);
  auto val2 = Legendre::Pn<float>(1, x);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == x);
}

// Pn(x) when x = 1
BOOST_AUTO_TEST_CASE(PnXequalsOne)
{
  auto n = iRan();

  auto val = Legendre::Pn<double>(n, 1.0);

  BOOST_TEST(val == 1.0);
}

// Pn(x) when x = -1
BOOST_AUTO_TEST_CASE(PnXequalsMinusOne)
{
  auto n1 = iRan();
  auto n2 = n1 + 1;

  auto val1 = Legendre::Pn<long double>(n1, -1.0);
  auto val2 = Legendre::Pn<long double>(n2, -1.0);

  BOOST_TEST(val1 == (n1 % 2 == 0 ? 1.0 : -1.0));
  BOOST_TEST(val2 == (n2 % 2 == 0 ? 1.0 : -1.0));
}

// Pn(x) when x = 0
BOOST_AUTO_TEST_CASE(PnXequalsZero)
{
  auto n = iRan();
  if ((n % 2) == 0) n++;

  auto val = Legendre::Pn<float>(n, 0.0);

  BOOST_TEST(val == 0.0);
}


BOOST_AUTO_TEST_CASE(Pn0)
{
  BOOST_TEST(Legendre::Pn<float>(2, 0.0) == -0.5);
  BOOST_TEST(Legendre::Pn<float>(4, 0.0) == 0.3750);
}

BOOST_AUTO_TEST_CASE(Pn05)
{
  BOOST_TEST(Legendre::Pn<float>(2, 0.05) == -0.4963,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(3, 0.05) == -0.0747,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(4, 0.05) == 0.3657,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(5, 0.05) == 0.0927,
                                  boost::test_tools::tolerance(tol));
}

BOOST_AUTO_TEST_CASE(Pn50)
{
  BOOST_TEST(Legendre::Pn<float>(2, 0.5) == -0.1250);
  BOOST_TEST(Legendre::Pn<float>(3, 0.5) == -0.4375);
  BOOST_TEST(Legendre::Pn<float>(4, 0.5) == -0.2891,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(5, 0.5) == 0.0898,
                                  boost::test_tools::tolerance(tol));
}

BOOST_AUTO_TEST_CASE(Pn95)
{
  BOOST_TEST(Legendre::Pn<float>(2, 0.95) == 0.8538,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(3, 0.95) == 0.7184,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(4, 0.95) == 0.5541,
                                  boost::test_tools::tolerance(tol));
  BOOST_TEST(Legendre::Pn<float>(5, 0.95) == 0.3727,
                                  boost::test_tools::tolerance(tol));
}
