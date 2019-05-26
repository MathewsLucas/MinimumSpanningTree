#pragma once
//#include <list> 
//#include <iterator>
/* edge -- pairs of cloud point indices with methods to
 *  - compute the distance from source to target
 *  - compare the lengths of two edges -- needed for sorting
 */

class edge {
  int  s, t;
  
public:
  double weight; 
  edge(int _s, int _t, int _weight) {
    s = _s;
    t = _t;
    weight = _weight;
  }

  int get_source() { return s; }
  int get_target() { return t; }
  int get_source() const { return s; }
  int get_target() const { return t; }
  double get_weight(){return weight;}
  static bool compare (edge e1, edge e2);
  static void print(std::list<edge> l);
};
struct Comp{
    bool operator()(const edge& a, const edge& b){
        return a.weight>b.weight;
    }
};