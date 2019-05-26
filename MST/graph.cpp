#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <queue>
//#include <mpi.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "graph.hpp"
#include "edge.hpp"
#include "neighborhood.hpp"
#include "unionFind.hpp"
#include "polyfit.h"

using namespace std;

/* graph -- method implementations */

graph::graph(int _size, ifstream *inFile, ifstream *weightFile)
{
  //Graph's constructor for exercise 1 and 2
  //it's input is the matrix of adjacencies and a file with the weights
  size = _size;
  string line;
  isEdge = new bool *[size];
  numEdges = 0;
  for (int i = 0; i < size; i++)
  {
    isEdge[i] = new bool[size];
    for (int j = 0; j < size; j++)
    {
      isEdge[i][j] = false;
    }
  }
  Adjlist = new neighborhood[_size];
  string wl;
  getline(*weightFile, wl); //weighline
  double _weight;
  stringstream sw(wl); //stringstream of the weight file
  if ((*inFile).is_open())
  {
    int countEdge = 0;
    getline(*inFile, line);
    while (getline(*inFile, line))
    {
      stringstream s(line);
      string word;
      int count = 0;
      int _s, _t;
      while (s >> word)
      {
        if (count == 0)
        {
          stringstream(word) >> _s;
          _s--;
          //std::cout << _s << "\n";
        }
        if (count >= 2)
        { //We ignore the 2 first elements of the string ( "i (level)")
          //stringstream(wl) >> _weight;
          if (count == 2)
            word = word.substr(1);
          stringstream(word) >> _t;
          _t--;
          if (!isEdge[_s][_t])
          {
            sw >> _weight;
            edge *e = new edge(_s, _t, _weight);
            edge *e1 = new edge(_t, _s, _weight);
            Adjlist[_s].addEdge(e);
            Adjlist[_t].addEdge(e1);
            isEdge[_s][_t] = true;
            isEdge[_t][_s] = true;
            numEdges++;
          }
        }
        count++;
      }
    }
    (*weightFile).close();
    (*inFile).close();
  }
}

graph::graph(cloud &_c)
{
  //Graph's constructor from a cloud
  //It constructs a complete graph with the points of the cloud _C

  c = &_c;
  //edge::set_cloud(c);
  size = c->get_n();
  Adjlist = new neighborhood[size];

  isEdge = new bool *[size];
  for (int i = 0; i < size; i++)
    isEdge[i] = NULL;
  for (int i = 0; i < size; i++)
  {
    isEdge[i] = new bool[size];
    for (int j = i + 1; j < size; j++)
    {
      if (isEdge[j] == NULL)
        isEdge[j] = new bool[size];
      edge *e = new edge(i, j);
      edge *e1 = new edge(j, i);
      Adjlist[i].addEdge(e);
      Adjlist[j].addEdge(e1);
      isEdge[i][j] = true;
      isEdge[j][i] = true;
    }
    Adjlist[i].neighbors.sort(edge::compare);
  }
}

double graph::getTotalWeight(std::list<edge> l)
{
  std::list<edge>::iterator it;
  double w = 0;
  for (it = l.begin(); it != l.end(); it++)
    w += it->get_weight();
  return w;
}

//Task 1

std::list<edge> graph::Boruvska()
{
  //Implementation in C++ of Boruvska
  //The return is a list of the edges that compounds the MST

  std::list<edge> Edges;
  std::list<int> leaders;
  for (int i = 0; i < size; i++)
  {
    leaders.push_back(i);
  }
  unionFind *parent = new unionFind(get_size());
  edge **minEdges;
  int aux;
  bool flag = true;                        //if we  7.7887n't have edges that can expand our forest
  while (parent->nbOfGroups() > 1 && flag) //if its 7.7887quals to 1 we already have a ST
  {
    //initialise minEdges
    minEdges = new edge *[get_size()];
    for (int i = 0; i < get_size(); i++)
      minEdges[i] = NULL;
    flag = false;

    for (int i = 0; i < get_size(); i++)
    {
      std::list<edge>::iterator it = Adjlist[i].neighbors.begin();
      for (; it != Adjlist[i].neighbors.end(); it++)
      {
        if (parent->find(it->get_target()) != parent->find(it->get_source()))
        {
          flag = true;
          if (minEdges[parent->find(it->get_source())] == NULL)
          {
            minEdges[parent->find(it->get_source())] = &(*it);
          }
          else if (minEdges[parent->find(it->get_source())]->get_weight() > it->get_weight())
          {
            //add direct
            minEdges[parent->find(it->get_source())] = &(*it);
          }
        }
      }
    }

    //add all minimums to the answer
    //and do union
    aux = parent->nbOfGroups();
    std::list<int>::iterator leaders_it = leaders.begin();
    for (; leaders_it != leaders.end(); leaders_it++)
    {
      if ((minEdges[*leaders_it] != NULL) && parent->find(minEdges[*leaders_it]->get_target()) != parent->find(minEdges[*leaders_it]->get_source()))
      {
        leaders.remove(parent->find(minEdges[*leaders_it]->get_target()));
        parent->merge(minEdges[*leaders_it]->get_source(), minEdges[*leaders_it]->get_target());
        Edges.push_back(*minEdges[*leaders_it]);
      }
    }

    //manual garbage collector
    delete minEdges;
  }
  delete parent;
  Edges.sort(edge::compare);
  return Edges;
}

int graph::existAFalse(bool *visited)
{
  for (int i = 0; i < get_size(); i++)
    if (!visited[i])
      return i;
  return -1;
}

std::list<edge> graph::Prim()
{
  /*
  * escolhemos um vertice, pegamos todas as arestas q saem desse vertice em ordem crescente, marcar o vertice como visitado
  * tiramos a menos aresta, se a target dela n estiver na union do vertice essa aresta pertence a ST,
  * mergemos o novo vertice a union, add as arestas q saem do novo vertice a pq e marcamos ele como visitado
  * condicao de parada, todos os vertices terem sido visitados
  */

  //Variables
  std::list<edge> Edges;
  unionFind *parent = new unionFind(get_size());
  priority_queue<edge, vector<edge>, Comp> pq;
  std::list<int> l;
  int aux;
  edge auxEdge;

  bool *visited = new bool[get_size()];
  for (int i = 0; i < get_size(); i++)
    visited[i] = false;

  aux = existAFalse(visited);

  while (aux != -1)
  {
    //add all sorting edges from existAFalse to pq
    if (!l.empty())
    {
      aux = l.front();
      l.pop_front();
    }
    std::list<edge>::iterator it = Adjlist[aux].neighbors.begin();

    for (; it != Adjlist[aux].neighbors.end(); it++)
      pq.push(*it);
    visited[aux] = true;

    //puts the smaller in the collection AND refresh visited
    while (!pq.empty())
    {
      auxEdge = pq.top();
      pq.pop();
      if (parent->find(auxEdge.get_source()) != parent->find(auxEdge.get_target()))
      {
        Edges.push_back(auxEdge);
        parent->merge(auxEdge.get_source(), auxEdge.get_target());
        l.push_back(auxEdge.get_target());
        break;
      }
    }
    aux = existAFalse(visited);
  }
  Edges.sort(edge::compare);
  return Edges;
}

//end Task 1

//Task 2 : Kruskal implementation

std::list<edge> *graph::getAllEdges()
{

  //Return a list with all the edges of the graph

  list<edge> *res = new list<edge>();
  for (int i = 0; i < size; i++)
  {
    std::list<edge>::iterator it;
    for (it = Adjlist[i].neighbors.begin(); it != Adjlist[i].neighbors.end(); it++)
      res->push_back(*it);
  }
  return res;
}
std::list<edge> *graph::Kruskal()
{
  unionFind *parent = new unionFind(get_size());
  //Implementation in C++ of Kruskal
  //The return is a list of the edges of the MST

  list<edge> *le = getAllEdges();

  list<edge> *ce = new list<edge>(); //Chosen edges
  le->sort(edge::compare);

  std::list<edge>::iterator it;
  for (it = le->begin(); it != le->end(); it++)
    if (parent->find(it->get_source()) != parent->find(it->get_target()))
    {
      parent->merge(it->get_source(), it->get_target());
      if (it->get_source() < it->get_target())
        ce->push_back(*it);
      else
        ce->push_back(edge(it->get_target(), it->get_source(), it->get_weight()));
    }
  ce->sort(edge::compare);
  return ce;
}

//End task 2.

//Task 5

std::list<std::list<edge>>* graph::kClustering(int k, cloud *c)
{
  std::list<edge> *poped_edges = new list<edge>();
  std::list<edge> *ce = Kruskal();
  popEdgesGetClusters(ce, poped_edges, k, c);
  return getClusters(ce);
}

std::list<std::list<edge>> *graph::getClusters(list<edge> *ce)
{
  //Given a Forest (aka ce) getClusters returns the list of the connected elements (clusters)
  //THe return tipe is a list of (list of edges), each elements of the list is a connected element

  std::list<std::list<edge>> *ret = new std::list<std::list<edge>>();
  unionFind *parent = new unionFind(get_size());
  std::list<edge>::iterator it;
  for (it = ce->begin(); it != ce->end(); it++)
    if (parent->find(it->get_source()) != parent->find(it->get_target()))
      parent->merge(it->get_source(), it->get_target());
  std::map<int, std::list<edge>> edgelist;
  std::list<edge>::iterator it2;
  for (it2 = ce->begin(); it2 != ce->end(); it2++)
  {
    int root = parent->find(it2->get_source());
    //check if key is not present
    if (edgelist.find(root) == edgelist.end())
    {
      std::list<edge> aux;
      aux.push_back(*it2);
      edgelist.insert(std::pair<int, std::list<edge>>(root, aux));
    }
    else
      edgelist.find(root)->second.push_back(*it2);
  }

  std::map<int, std::list<edge>>::iterator itr;
  for (itr = edgelist.begin(); itr != edgelist.end(); itr++)
    if (itr->second.size() > 5)
      ret->push_front(itr->second);
  return ret;
}

bool graph::popEdgesGetClusters(list<edge> *ce, list<edge> *poped_edges, int k, cloud *c)
{
  //Receives a list of edges (previously ordered) ce and the number of cluster desired.
  //Pop the edges of ce (that must be ordered) to get k clusters (eliminating the not relevant points)
  //and keeps the poped elements of ce in poped_edges

  std::list<std::list<edge>> *ret = getClusters(ce);

  while (ret->size() < k)
  {
    if (ce->size() > 0)
    {
      poped_edges->push_back(ce->back());
      ce->pop_back();
      ret->clear(); //We clear the space allocated that is no longer used
      ret = getClusters(ce);
    }
    else
      return false;
  }
  return true;
}
void graph::pushBackAllEdges(list<edge> *ce, list<edge> *poped_edges)
{
  //push back all the edges of ce that were poped(kept in poped_edges)
  std::list<std::list<edge>> *ret = getClusters(ce);
  while (poped_edges->size() > 0)
  {
    ce->push_back(poped_edges->back());
    poped_edges->pop_back();
  }
}

std::list<std::list<edge>> *graph::getBestk(cloud *c, int max_k)
{
  std::list<edge> *ce = Kruskal();
  //Receives a cloud of points (c) and returns a list of connected elements that optimizes the
  //utilityFunction.
  std::list<std::list<edge>> *ret = new std::list<std::list<edge>>();
  std::list<edge> *poped_edges = new list<edge>();
  double U = __DBL_MAX__;
  double u;
  int best_k = 1;
  c->set_label();
  for (int i = 1; i <= max_k; i++)
  {
    if (!popEdgesGetClusters(ce, poped_edges, i, c))
    {
      break; //we couldn't get any more clusters
    }
    ret->clear(); //We clear the space allocated that is no longer used
    ret = getClusters(ce);
    c->set_k(i);
    setLabel(ret, c);

    u = c->utilityFunction();
    if (u < U)
    {
      U = u;
      best_k = i;
    }
  }
  pushBackAllEdges(ce, poped_edges);
  popEdgesGetClusters(ce, poped_edges, best_k, c);
  ret->clear(); //We clear the space allocated that is no longer used
  ret = getClusters(ce);

  c->set_k(best_k);

  c->set_label(); //Clear the labels just in case
  setLabel(ret, c);
  c->set_centroid_centers();
  return ret;
}

void graph::setLabel(std::list<std::list<edge>> *ret, cloud *c)
{
  //Receives the list of connected elements (ret) of a cloud (c) and set the labels of each connected
  //elements of c

  int l = 1; // l=0 will caracterize isolate points
  int k = ret->size();
  std::list<std::list<edge>>::iterator it;
  std::list<edge>::iterator itr;
  for (it = ret->begin(); it != ret->end(); it++)
  {
    for (itr = it->begin(); itr != it->end(); itr++)
    {
      c->get_point(itr->get_source()).label = l;
      c->get_point(itr->get_target()).label = l;
    }
    l = (l + 1);
  }
}

void graph::refineClustering(cloud *c, list<edge> *ce)
{
  //Can be used to refine the clustering
  //this function will try to find pairs of clusters really close to each others (that should be only
  //one cluster) and will merge them together.

  double dist_moyenne = c->intercluster_distance_m();
  double intercluster_std = c->intercluster_std();
  std::cout << "Dist " << dist_moyenne << " STD " << intercluster_std << "\n";
  c->set_representants();
  std::cout << dist_moyenne - (0.8) * intercluster_std << "\n";
  for (int i = 1; i < c->get_k(); i++)
    for (int j = i + 1; j <= c->get_k(); j++)
      if ((c->get_center(i)).dist(c->get_center(j)) < dist_moyenne - (0.8) * intercluster_std)
        //They actually got to be in the same cluster
        ce->push_back(edge(c->get_representant(i), c->get_representant(j)));
}

std::list<std::list<edge>> *graph::MSDR(cloud *c, int max_k, std::list<edge> *ce)
{
  //This method is a implementation of the Maximum Standard Deviation Reduction Clustering Algorithm (MSDR)
  //that can be found in the article: https://ieeexplore.ieee.org/document/4031882

  std::list<std::list<edge>> *S = new std::list<std::list<edge>>();
  std::list<edge> *removed_edges = new list<edge>();
  std::list<edge>::iterator it, it_aux;

  edge e;
  int k = 0, flag = 0;
  double sigmak = 0.0;
  double eps = 0.0001;
  double temp = 0.0;
  double StdDevRed[100];

  for (int i = 0; i < 100; i++)
    StdDevRed[i] = 0;
  /*double graph::getTotalWeight(std::list<edge> l)
{
  std::list<edge>::iterator it;
  double w = 0;
  for (it = l.begin(); it != l.end(); it++)
    w += it->get_weight();
  return w;
}*/
  int Nclusters;
  int i = 0;
  S = getClusters(ce);
  sigmak = edge::std_deviation(*S);
  Nclusters = S->size();
  do
  {
    i++;
    temp = sigmak;

    //Choose an edge that leads to max StdDev reduction and remove it
    e = edge(0, 0);
    int count = 0;
    int N = ce->size();
    flag = 0;
    for (it = ce->begin(); count < N; count++)
    {
      edge beRemoved = edge(it->get_source(), it->get_target());
      it = ce->erase(it);
      S->clear();
      S = getClusters(ce);

      sigmak = edge::std_deviation(*S);
      if (StdDevRed[i] < temp - sigmak)
      {
        StdDevRed[i] = temp - sigmak; //temp > sigmak
        e = beRemoved;
        flag = 1;
        ce->push_back(beRemoved);
        it_aux = ce->end();
        it_aux--;
      }
      else
      {
        ce->push_back(beRemoved);
      }
    }
    ce->erase(it_aux);
    removed_edges->push_back(e);
    sigmak = temp - StdDevRed[i];
    Nclusters++;
    //std::cout << "StdDevRed[" << i << "]: " << StdDevRed[i] << "\n";

  } while (abs(StdDevRed[i] - StdDevRed[i - 1]) > abs(eps * (StdDevRed[i] + 1)));
  std::cout << removed_edges->size() << " " << i << "\n";
  k = 0;
  if (i == 1)
    k = 2;
  else if (i == 2)
  {
    if (StdDevRed[1] < StdDevRed[2])
      k = 2;
    else
      k = 3;
  }
  else
    k = computeMinima(StdDevRed, i);
  if (k == 0)
    k = i;

  int count = 1;
  for (it = removed_edges->begin(); it != removed_edges->end(); it++, count++)
    if (count >= k)
      ce->push_back(*it);

  S->clear();
  S = getClusters(ce);
  k = S->size();

  c->set_k(k);

  c->set_label(); //Clear the labels just in case
  setLabel(S, c);

  c->set_centroid_centers();
  return S;
}

std::list<std::list<edge>> *graph::MSDR(cloud *c, int max_k)
{
  std::list<edge> *ce = Kruskal();
  return MSDR(c, max_k, ce);
}

int graph::computeMinima(double *y, int N)
{
  //Implementation of polyfit that returns the indice correspondind to the first minimum of
  //the polynomio function

  int order = 3;
  double eps = 0.001;
  double xData[N];
  double yData[N];
  for (int i = 0; i < N; i++)
  {
    yData[i] = y[i + 1];
    xData[i] = i;
  }
  double coefficients[order + 1];
  polyfit(xData, yData, N, order, coefficients);

  /*for (int i = 0; i < order + 1; i++)
  {
    std::cout << coefficients[i] << " ";
  }
  std::cout << "\n";*/

  double dx = 0.01;
  int Npts = (N - 1) / dx + 1;
  double xFit[Npts];
  double yFit[Npts];

  xFit[0] = 0;
  for (int i = 0; i < Npts; i++)
  {
    if (i > 0)
      xFit[i] = xFit[i - 1] + dx;
    yFit[i] = 0;
    for (int r = 0; r < order + 1; r++)
      yFit[i] += coefficients[r] * pow(xFit[i], r);
  }

  for (int i = 0; i < Npts - 2; i++)
    if ((yFit[i] > yFit[i + 1]) && (yFit[i + 1] < yFit[i + 2]))
    {
      i++;
      std::cout << xFit[i] << "\n";
      if (abs(ceil(xFit[i]) - xFit[i]) < abs(floor(xFit[i]) - xFit[i]))
      {
        return ceil(xFit[i]) + 1;
      }
      else
      {
        return floor(xFit[i]) + 1;
      }
    }
  return N;
}

void graph::graphToFile(std::ofstream *outfile)
{
  (*outfile) << numEdges << '\n';
  for (int i = 0; i < size; i++)
  {
    std::list<edge>::iterator it;
    it = Adjlist[i].neighbors.begin();
    for (; it != Adjlist[i].neighbors.end(); it++)
    {
      (*outfile) << it->get_source() << " " << it->get_target() << " " << it->get_weight() << "\n";
    }
  }
}
/*
void graph::graphToFile(std::ofstream *outfile){
  for (int i = 0; i < size; i++)
  {
    std::list<edge>::iterator it;
    it = Adjlist[i].neighbors.begin();
    std::cout << it->get_source() << " [" << it->get_target();
    it++;
    for (; it != Adjlist[i].neighbors.end(); it++)
    {
      std::cout << ", " <<it->get_target(); 
    }
    std::cout << "]\n";

  }
}*/
//End task 5.
