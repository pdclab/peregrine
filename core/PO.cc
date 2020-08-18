#include <iostream>
#include <set>
#include <numeric>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

#include <stdint.h>

#include "./bliss-0.73/graph.hh"
#include "utils.hh"

namespace PO
{
  std::string to_string(std::vector<uint32_t> p)
  {
    std::string res;
    for (auto v : p)
    {
      res += std::to_string(v);
    }

    return res;
  }

  void addProducts(std::vector<std::vector<uint32_t>> *permgroup)
  {
    // multiplying std::vector<uint32_t>s is just function composition: (p1*p2)[i] = p1[p2[i]]
    std::vector<std::vector<uint32_t>> products;
    // for filtering duplicates
    std::unordered_set<std::string> dups;
    for (auto f : *permgroup)
    {
      dups.insert(to_string(f));
    }

    for (auto k = permgroup->begin(); k != permgroup->end(); k++)
    {
      for (auto l = permgroup->begin(); l != permgroup->end(); l++)
      {
        //if (k >= l) continue;
        std::vector<uint32_t> p1 = *k;
        std::vector<uint32_t> p2 = *l;

        std::vector<uint32_t> product;
        product.resize(p1.size());
        for (unsigned i = 0; i < product.size(); i++)
        {
          product[i] = p1[p2[i]];
        }

        // don't count duplicates
        if (dups.count(to_string(product)) == 0)
        {
          dups.insert(to_string(product));
          products.push_back(product);
        }
      }
    }

    for (auto p : products)
    {
      permgroup->push_back(p);
    }
  }

  std::vector<std::vector<uint32_t>> getAutomorphisms(bliss::Graph H)
  {
    std::vector<std::vector<uint32_t>> result;
    bliss::Stats s;
    H.find_automorphisms(s, [](void *resultp, uint32_t aut_sz, const uint32_t *aut) {
      std::vector<uint32_t> _aut;
      for (uint32_t i = 0; i < aut_sz; i++)
      {
        _aut.push_back(aut[i]);
      }
      ((std::vector<std::vector<uint32_t>> *)resultp)->push_back(_aut);
    }, &result);

    // add products until full group is generated
    uint32_t counter = 0;
    uint32_t lastSize = 0;
    while (result.size() != lastSize)
    {
      lastSize = result.size();

      addProducts(&result);

      counter++;
      if (counter > 100)
      {
        break;
      }
    }

    return result;
  }


  // std::vector<uint32_t>s are positional: node n0 is mapped to std::vector<uint32_t>[0]
  // A is not necessarily Aut(H).
  // Consider the std::vector of std::vector<uint32_t>s as a table, then columns are equivalence classes with respect to isomorphism.
  // An A-equivalence class consists of vertices that occur in each other's columns
  auto getAEquivalenceClasses(const std::vector<std::vector<uint32_t>> &A, const bliss::Graph &H)
  {
    std::map<uint32_t, std::set<uint32_t>> eclasses;

    for (unsigned int id = 0; id < H.get_nof_vertices(); ++id)
    {
      std::set<uint32_t> eclass;
      for (auto &&perm : A)
      {
        eclass.insert(perm[id]);
      }

      uint32_t rep = *std::min_element(eclass.cbegin(), eclass.cend());
      eclasses[rep].insert(eclass.cbegin(), eclass.cend());
    }

    return eclasses;
  }

  std::vector<std::pair<uint32_t, uint32_t>> findPOs(const bliss::Graph &H)
  {
    // compute Aut_H and H_E
    auto Aut_H = getAutomorphisms(H);

    std::vector<std::pair<uint32_t, uint32_t>> result;
    std::vector<std::vector<uint32_t>> C(H.get_nof_vertices());

    auto A = Aut_H;
    auto eclasses = getAEquivalenceClasses(A, H);
    // Don't pick largest: want conditions on what will likely be the core
    // vertices first, so that there are no conditions across sibling groups.
    auto eclass_it = std::find_if(eclasses.cbegin(), eclasses.cend(), [](auto &&a1) { return a1.second.size() > 1; });
    while (eclass_it != eclasses.cend() && eclass_it->second.size() > 1)
    {
      const auto &eclass = eclass_it->second;

      //std::cout << "A:" << std::endl;
      //for (auto &&f : A) utils::print_vector(f);

      //std::cout << "eclasses:" << std::endl;
      //for (auto &&s : eclasses)
      //{
      //  std::cout << s.first << ": ";
      //  for (uint32_t v : s.second) std::cout << v << " ";
      //  std::cout << std::endl;
      //}

      // arbitrary element of the equivalence class
      uint32_t n0 = *eclass.cbegin();

      for (auto &&f : A)
      {
        // there must exist a min != n0 since eclass_it->size() > 1
        uint32_t min = *std::min_element(std::next(eclass.cbegin()), eclass.cend(), [f](uint32_t n, uint32_t m) { return f[n] < f[m]; });

        // AnalyzedPattern expects 1-based vertices
        result.emplace_back(n0+1, min+1);
      }

      A.erase(std::remove_if(A.begin(), A.end(), [n0](auto &&L)
            {
              return L[n0] != n0;
            }), A.end());


      eclasses = getAEquivalenceClasses(A, H);
      eclass_it = std::find_if(eclasses.cbegin(), eclasses.cend(), [](auto &&a1) { return a1.second.size() > 1; });
    }

    // likely many duplicate conditions: get rid of them
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());

    return result;
  }

  std::vector<std::vector<uint32_t>> automorphicSets(const bliss::Graph &H)
  {
    std::vector<std::vector<uint32_t>> result;
    uint32_t n = H.get_nof_vertices();
    std::vector<bool> visited(n, false);
    auto automorphisms = getAutomorphisms(H);

    if (automorphisms.empty())
    {
      result.resize(n);
      for (uint32_t i = 0; i < n; ++i)
      {
        result[i] = {i};
      }
    }
    else
    {
      for (uint32_t i = 0; i < n; ++i)
      {
        if (visited[i]) continue;
        std::set<uint32_t> set;
        for (const auto &p : automorphisms)
        {
          set.insert(p[i]);
          visited[p[i]] = true;
        }
        result.emplace_back(set.cbegin(), set.cend());
      }
    }

    return result;
  }
}
