#include "edge.hpp"

/* edge -- method implementations */

bool edge::compare(edge e1, edge e2)
{
  if (e1.weight < e2.weight)
    return true;
  else if (e1.weight == e2.weight)
  {
    if (e1.s < e2.s)
      return true;
    else if (e1.s == e2.s)
    {
      if (e1.t < e2.t)
        return true;
    }
  }
  return false;
}
