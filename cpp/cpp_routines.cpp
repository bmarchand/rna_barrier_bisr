#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <map> 
#include <list> 
#include <utility> 
#include <iostream> 

using namespace std;

int LARGE_VALUE = 1 << 20;

bool BFS_hopcroft_karp(map<int, list<int>> &ngbh,
                       list<int> &left,
                       list<int> &right,
                       map<int, int> &match_of_left,
                       map<int, int> &match_of_right,
                       map<int, int> &distances)
{
    

    list<int> queue;

    for (auto const & u : left)
    {
        if (match_of_left[u]==-1)
        {
            distances[u] = 0;
            queue.push_back(u);
        }
        else
        {
            distances[u] = LARGE_VALUE;
        }
    }
    distances[-1] = LARGE_VALUE;

    while (queue.size() > 0)
    {
        int u = queue.front();
        queue.pop_front();
        if (distances[u] < distances[-1])
        {
            for (auto const & v : ngbh[u])
            {
                if (distances[match_of_right[v]] == LARGE_VALUE)
                {
                    distances[match_of_right[v]] = distances[u] + 1;
                    queue.push_back(match_of_right[v]);
                }
            }
        }
    }

    if (distances[-1] < LARGE_VALUE)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool DFS_hopcroft_karp(int u,
                       map<int, list<int>> ngbh,
                       map<int, int> &distances,
                       map<int, int> &match_of_right,
                       map<int, int> &match_of_left)

{
    if (u >= 0)
    {
        for (auto const & v : ngbh[u])
        {
            if (distances[match_of_right[v]] == distances[u] + 1)
            {
                if (DFS_hopcroft_karp(match_of_right[v],
                                      ngbh,
                                      distances,
                                      match_of_right,
                                      match_of_left))
                {
                    match_of_right[v] = u;
                    match_of_left[u] = v;
                    return true;
                }
            }
        }
        distances[u] = LARGE_VALUE;
        return false;
    }
    return true;
}

list<pair<int, int>> hopcroft_karp(map<int, list<int>> ngbh,   // adjacency
                                   list<int> left,             // one side
                                   list<int> right)            // other side
{

    map<int, int> match_of_left;  // keys: left, values: the match in right
    map<int, int> match_of_right; // keys: right, values: the match in left

    for (auto const& u : left) match_of_left[u] = -1; // init as unmatched
    for (auto const& u : right) match_of_right[u] = -1; 

    map<int, int> distances;

    bool keep_going = BFS_hopcroft_karp(ngbh, 
                                        left,
                                        right,
                                        match_of_left,
                                        match_of_right,
                                        distances);

    // main loop
    while (keep_going)
    {
        for (auto const& u : left)
        {
            if (match_of_left[u] == -1) 
            {
                DFS_hopcroft_karp(u,
                                  ngbh,
                                  distances,
                                  match_of_right,
                                  match_of_left);
            }            
        }
        keep_going = BFS_hopcroft_karp(ngbh, 
                                       left,
                                       right,
                                       match_of_left,
                                       match_of_right,
                                       distances);
    }
    
    list<pair<int,int>> matching; // the object that will be returned

    for (auto const & u : left)
    {
        if (match_of_left[u] != -1)
        {
            pair<int, int> match_edge;
            match_edge.first = u;
            match_edge.second = match_of_left[u];
            matching.push_back(match_edge);
        } 
    }


    return matching;
}

int add(int i, int j) {
    return i + j;
}

void DFS_visit(map<int, list<int>> out_ngbh,
               int u,
               list<int> &container,
               map<int, int> &color)
{
    color[u] = 1; 
    for (auto const& v : out_ngbh[u])
    {
        if (color[v] == 0) 
        {
            DFS_visit(out_ngbh, 
                      v,
                      container,
                      color);
        }
    }
    container.push_back(u);
}

void DFS(map<int, list<int>> out_ngbh,
         list<int> &finishing_order,
         list<list<int>> &SCC,
         bool transpose=false)
{
    // In both cases, start by coloring all vertives white
    map<int, int> color;

    for (auto const& pair : out_ngbh)
    {
        color[pair.first] = 0; // white
    }

    // If we are in the first, forward "finishing_times computing" case:
    if (!transpose)
    {
        for (auto const& pair : out_ngbh)
        {
            if (color[pair.first] == 0) 
            {
                DFS_visit(out_ngbh, 
                          pair.first, 
                          finishing_order,
                          color);
            }
        }
    }
    // Else we are in the "reverse, finishing_order-based SCC computing case"
    else
    {
        list<int>::reverse_iterator rev_iter;
        for (rev_iter = finishing_order.rbegin(); 
             rev_iter != finishing_order.rend();
             rev_iter++)
        {

            list<int> scc;

            if (color[*rev_iter] == 0) 
            {
                DFS_visit(out_ngbh, 
                          *rev_iter, 
                          scc,
                          color);
            }
   
            if (scc.size() > 0) SCC.push_back(scc);
        }

    }
}



list<list<int>> strongly_connected_components(map<int, list<int>> out_ngbh,
                                   map<int, list<int>> in_ngbh)
{
    
    list<int> finishing_order;
    list<list<int>> SCC; 

    // first call: computes finishing_order, SCC untouched.
    DFS(out_ngbh, finishing_order, SCC); 

    // second call: uses finishing_order to iterate, populate SCC
    DFS(in_ngbh, finishing_order, SCC, true);

    // The SCCs are returned in a topological order.
    return SCC;
}

PYBIND11_MODULE(bisr_dpw_cpp_routines, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function which adds two numbers");

    m.def("strongly_connected_components", 
          &strongly_connected_components, 
          "SCC");

    m.def("hopcroft_karp",
          &hopcroft_karp,
          "maximum matching algorithm");
}
