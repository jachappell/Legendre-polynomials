#include "legendre.h"

#include <memory>
#include <random>

#define BOOST_TEST_MODULE Legendre
#include <boost/test/included/unit_test.hpp>

using namespace Storage_B;

namespace
{
  class Random 
  {
  public:
    Random()
    {
      std::random_device rd;
      _gen = std::make_shared<std::mt19937>(rd());
    }
  
    // no copy
    Random(const Random&) = delete;
    Random& operator=(const Random&) = delete;

    ~Random() = default;

    template <class T> auto one_to_one()
    {
      return random(static_cast<T>(-1), static_cast<T>(1));
    }

    template <class T> auto zero_to_one()
    {
      return random(static_cast<T>(0), static_cast<T>(1));
    }

    template <class T> auto random(T low, T high)
    {
      std::uniform_real_distribution<T> dis(low, high);

      return dis(*_gen);
    }

    template <class T> auto random_int(T low, T high)
    {
      std::uniform_int_distribution<T> dis(low, high);

      return dis(*_gen);
    }

  private:
    std::shared_ptr<std::mt19937> _gen;
  };

  Random ran;

  auto tol = 0.00049;
}

// P0
BOOST_AUTO_TEST_CASE(P0)
{
  auto val1 = Legendre::P0<double>(ran.one_to_one<double>());
  auto val2 = Legendre::Pn<double>(0, ran.one_to_one<double>());

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);
}

// P1
BOOST_AUTO_TEST_CASE(P1)
{
  auto x = ran.one_to_one<float>();

  auto val1 = Legendre::P1<float>(x);
  auto val2 = Legendre::Pn<float>(1, x);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == x);
}

// Pn(x) when x = 1
BOOST_AUTO_TEST_CASE(PnXequalsOne)
{
  auto n = ran.random_int<unsigned int>(2, 10);

  auto val = Legendre::Pn<double>(n, 1.0);

  BOOST_TEST(val == 1.0);
}

// Pn(x) when x = -1
BOOST_AUTO_TEST_CASE(PnXequalsMinusOne)
{
  auto n1 = ran.random_int<unsigned int>(2, 10);
  auto n2 = n1 + 1;

  auto val1 = Legendre::Pn<long double>(n1, -1.0);
  auto val2 = Legendre::Pn<long double>(n2, -1.0);

  BOOST_TEST(val1 == (n1 % 2 == 0 ? 1.0 : -1.0));
  BOOST_TEST(val2 == (n2 % 2 == 0 ? 1.0 : -1.0));
}

// Pn(x) when x = 0
BOOST_AUTO_TEST_CASE(PnXequalsZero)
{
  auto n = ran.random_int<unsigned int>(2, 10);
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
