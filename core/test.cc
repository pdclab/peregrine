#include <UnitTest++/UnitTest++.h>

#include "unittests/Graph_test.hh"
#include "unittests/PatternGenerator_test.hh"
#include "unittests/PatternMatching_test.hh"

#include "integrationtests/Counting_test.hh"
#include "integrationtests/Matching_test.hh"
#include "integrationtests/DataConverter_test.hh"
#include "integrationtests/EarlyTermination_test.hh"
#include "integrationtests/OnTheFly_test.hh"

int main()
{
  std::cout << "#THREADS: " << std::thread::hardware_concurrency() << std::endl;
  return UnitTest::RunAllTests();
}
