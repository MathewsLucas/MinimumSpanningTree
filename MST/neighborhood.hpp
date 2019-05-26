#pragma once
//#include <list> 
//#include <iterator>
#include "edge.hpp"

using namespace std; 
/*  DESCRIÇÃO  */
class neighborhood {
  public: 
  list <edge> neighbors;
  
  int size;
  neighborhood();
  ~neighborhood();
  void addEdge(edge * _e);
  int get_size();
  void print(); //for testing
};