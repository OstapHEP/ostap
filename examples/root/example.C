//usr/bin/env root.exe -l ${0}\(\""${0}"\",\""${*}"\"\); exit $?

// #include "Ostap/ValueWithError.h"
{
  Ostap::Math::ValueWithError b(2,1);
  std::cout << b << std::endl;
  return 0;
}
