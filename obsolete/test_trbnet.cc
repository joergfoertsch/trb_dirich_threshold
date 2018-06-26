#include "trbnet.h"
#include "trberror.h"
#include <iostream>

int main (int arc, char** argv){
  std::cout << init_ports() << trb_strerror()  << std::endl;
}
