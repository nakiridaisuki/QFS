#pragma once

#include <string>
class Test {
  public:
    Test(std::string name) : name(name) {}

    void SayHey();

  private:
    std::string name;
};
