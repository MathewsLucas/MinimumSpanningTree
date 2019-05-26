#include "unionFind.hpp"

unionFind::unionFind(int n)
{
    size = n;
    parents = new int[n];
    isLeader = new bool[n];
    for (int i = 0; i < n; i++)
    {
        parents[i] = i;
        isLeader[i] = true;
    }
    groups = n;
}

int unionFind::getSize(){
    return size;
}

int unionFind::find(int src)
{
    //Find must not return null
    //if(src == null) throw error
    int p = parents[src];

    while (p != parents[p])
        p = parents[p];

    parents[src] = p;
    return p;
}

void unionFind::merge(int v0, int v1)
{
    //TODO
    int Pv0 = find(v0);
    int Pv1 = find(v1);

    if (Pv0 != Pv1)
    {
        isLeader[Pv1] = false;
        parents[Pv1] = Pv0;
        groups--;
    }
}

cell** unionFind::leaders(){
    cell** leaders = new cell*[nbOfGroups()];
    int cont = 0;
    for(int i = 0; i < size; i++){
        if(isLeader[i]){
            leaders[cont] = new cell(i);
            cont++;
        }
    }

    return leaders;
}

int unionFind::nbOfGroups()
{
    return groups;
}

unionFind::~unionFind()
{
    delete parents;
}
