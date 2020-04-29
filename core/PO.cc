#include <iostream>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

#include <stdint.h>

#include "./bliss-0.73/graph.hh"

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

    void printPerm(std::vector<uint32_t> p)
    {
        //std::string index = "";
        std::string perm = "";
        for (uint32_t i = 0; i < p.size(); i++)
        {
            //index = index + std::to_string(i) + " ";
            perm = perm + std::to_string(p[i] + 1) + " ";
        }

        std::cout << perm << std::endl;
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
        //std::unordered_set<std::vector<uint32_t>, std::function<size_t(const std::vector<uint32_t> &)>> result(10, phash);
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
    auto getAEquivalenceClasses(std::vector<std::vector<uint32_t>> A, bliss::Graph H)
    {
        (void)H;

        // use map instead of unordered_map so the vertices are ordered.
        std::map<uint32_t, std::set<uint32_t>> eclasses;

        // a perm is a row
        for (auto perm : A)
        {
            // id is the column cursor
            uint32_t id = 0;
            // for each item of the current row
            for (auto v : perm)
            {
                if (id > v)
                {
                    break;
                }
                else
                {
                    eclasses[id].insert(v);
                    // steal v's equivalence class as well
                    if (id < v && eclasses.count(v) > 0)
                    {
                        for (auto u : eclasses.at(v))
                        {
                            eclasses[id].insert(u);
                        }
                    }
                }
                id++;
            }
        }

        // remove overlapping classes
        for (auto p : eclasses)
        {
            for (auto v : p.second)
            {
                if (v != p.first)
                {
                    eclasses.erase(v);
                }
            }
        }

        return eclasses;
    }

    // key of each equivalence class in the ordered map
    auto getRepresentatives(std::map<uint32_t, std::set<uint32_t>> eclasses)
    {
        std::vector<uint32_t> result;

        for (auto p : eclasses)
        {
            uint32_t key = p.first;
            result.push_back(key);
        }

        return result;
    }

    std::vector<std::pair<uint32_t, uint32_t>> findPOs(bliss::Graph H)
    {
        // compute Aut_H and H_E
        auto Aut_H = getAutomorphisms(H);

        auto H_E_representatives = getRepresentatives(getAEquivalenceClasses(Aut_H, H));

        std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> M;

        for (uint32_t n : H_E_representatives)
        {
            std::vector<std::pair<uint32_t, uint32_t>> C;
            uint32_t np = n;
            std::vector<std::vector<uint32_t>> A = Aut_H;
            //std::cout << "processing " << labels[n] << std::endl;
            while (A.size() != 1)
            {
                auto A_E_classes = getAEquivalenceClasses(A, H);
                //std::cout << "initial eclasses:" << std::endl;
                //for (auto p : A_E_classes) {
                //  std::cout << labels[p.first] << ": ";
                //  for (auto v : p.second) {
                //    std::cout << labels[v] << " ";
                //  }
                //  std::cout << std::endl;
                //}

                // if only one element in the A-equivalence class, there is no symmetry
                if (A_E_classes[np].size() > 1)
                {
                    // find min{label(m) | m ~A np and m != np}
                    // sorted, so 1st element is np, second is min
                    uint32_t min = *(++A_E_classes[np].begin());

                    // add 'label(np) < min{label(m) | m ~A np and m != np}' to C
                    //std::cout << "adding " << labels[np] << " < " << labels[min] << std::endl;
                    C.push_back({np, min});

                    // A <- {f in A | f^-1(np) < f^-1(min)}
                    std::vector<std::vector<uint32_t>> tmpA;
                    for (auto f : A)
                    {
                        ssize_t npinv = -1, mininv = -1;
                        for (uint32_t v = 0; v < f.size(); v++)
                        {
                            if (f[v] == np)
                            {
                                npinv = v;
                            }
                            else if (f[v] == min)
                            {
                                mininv = v;
                            }
                        }

                        assert(npinv != -1 && mininv != -1);
                        if (npinv < mininv)
                        {
                            tmpA.push_back(f);
                            //printPerm(f);
                        }
                    }
                    A = tmpA;
                }

                //std::cout << "picking largest" << std::endl;
                auto new_A_E_classes = getAEquivalenceClasses(A, H);

                //std::cout << "new eclasses:" << std::endl;
                //for (auto p : new_A_E_classes) {
                //  std::cout << labels[p.first] << ": ";
                //  for (auto v : p.second) {
                //    std::cout << labels[v] << " ";
                //  }
                //  std::cout << std::endl;
                //}

                // E <- largest A-equivalence class
                auto largest = new_A_E_classes[0];
                for (auto p : new_A_E_classes)
                {
                    auto eclass = p.second;
                    if (eclass.size() > largest.size())
                    {
                        largest = eclass;
                    }
                }

                // pick arbitrary np from E
                np = *largest.begin();
            }

            M[n] = C;
        }

        std::unordered_map<uint32_t, std::unordered_set<uint32_t>> orders;
        for (auto p : M)
        {
            uint32_t key;
            std::vector<std::pair<uint32_t, uint32_t>> conditions;
            std::tie(key, conditions) = p;
            //std::cout << "key " << key << std::endl;
            for (auto cond : conditions)
            {
                orders[cond.first].insert(cond.second);
            }
        }

        std::vector<std::pair<uint32_t, uint32_t>> partial_orders;
        for (auto p : orders)
        {
            for (auto v : p.second)
            {
                partial_orders.push_back({p.first + 1, v + 1});
            }
        }

        return partial_orders;
    }

    std::vector<std::vector<uint32_t>> automorphicSets(bliss::Graph H)
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


