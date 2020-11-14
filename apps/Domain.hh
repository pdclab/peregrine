#ifndef DOMAIN_HH
#define DOMAIN_HH

#include "roaring/roaring.hh"

struct Domain
{
  /**
   * Default constructor is needed by Peregrine.
   */
  Domain() : sets(Peregrine::Context::current_pattern->num_aut_sets())
  {}

  /**
   * 'Copy' constructor is needed by Peregrine.
   *
   * Note that this is not a true copy constructor since it copies not from
   * Domain but from what actually gets mapped in the callback in fsm.cc.
   * You could also use Domain instead of std::vector as the type for m, but
   * it might be less efficient.
   */
  Domain(const std::vector<uint32_t> &m) : sets(Peregrine::Context::current_pattern->num_aut_sets())
  {
    for (uint32_t i = 0; i < m.size(); ++i)
    {
      uint32_t u = m[i];
      uint32_t set = Peregrine::Context::current_pattern->aut_map[i];
      sets[set].add(u);
    }
  }

  /**
   * A sum operator is needed by Peregrine.
   */
  Domain &operator+=(const Domain &d)
  {
    uint32_t i = 0;
    for (const auto &set : d.get_sets())
    {
      sets[i++] |= set;
    }

    return *this;
  }

  /**
   * Specialization of the sum operator for efficiency.
   *
   * Avoids building Roaring bitsets with a single element just so they can
   * be added to another Domain. Add new elements directly.
   */
  Domain &operator+=(const std::vector<uint32_t> &m)
  {
    for (uint32_t i = 0; i < m.size(); ++i)
    {
      uint32_t u = m[i];
      uint32_t set = Peregrine::Context::current_pattern->aut_map[i];
      sets[set].add(u);
    }

    return *this;
  }

  /**
   * reset() method is required by Peregrine.
   *
   * Resets the data structure for the next pattern to be matched.
   */
  void reset()
  {
    sets.clear();
    sets.resize(Peregrine::Context::current_pattern->num_aut_sets());
  }

  /**
   * Returns trivially-copyable type for the view function.
   */
  uint64_t get_support()
  {
    uint64_t s = sets.front().cardinality();
    for (uint32_t i = 1; i < sets.size(); ++i)
    {
      s = std::min(s, sets[i].cardinality());
    }

    return s;
  }

  const std::vector<Roaring> &get_sets() const
  {
    return sets;
  }

  std::vector<Roaring> sets;
};


// label discovery with one or two edges: num_edges should be either 1 or 2
template <unsigned int num_edges>
struct DiscoveryDomain
{
  DiscoveryDomain() : sets(num_edges)
  {}

  DiscoveryDomain(const std::pair<std::vector<uint32_t>, uint32_t> &p) : sets(num_edges + p.second)
  {
    auto &m = p.first;

    sets[0].add(m[0]);
    if constexpr (num_edges == 2) sets[1].add(m[1]);
    sets[p.second].add(m[num_edges]); // if merge
  }

  DiscoveryDomain &operator+=(const DiscoveryDomain &d)
  {
    sets.resize(d.sets.size());

    uint32_t i = 0;
    for (const auto &set : d.get_sets())
    {
      sets[i++] |= set;
    }

    return *this;
  }

  DiscoveryDomain &operator+=(const std::pair<std::vector<uint32_t>, uint32_t> &p)
  {
    sets.resize(1+p.second);
    auto &m = p.first;
    sets[0].add(m[0]);
    if constexpr (num_edges == 2) sets[1].add(m[1]);
    sets[p.second].add(m[num_edges]); // if merge

    return *this;
  }

  uint64_t get_support()
  {
    uint64_t s = sets.front().cardinality();
    for (uint32_t i = 1; i < sets.size(); ++i)
    {
      s = std::min(s, sets[i].cardinality());
    }

    return s;
  }

  void reset()
  {
    sets.clear();
    sets.resize(Peregrine::Context::current_pattern->num_aut_sets());
  }

  const std::vector<Roaring> &get_sets() const
  {
    return sets;
  }

  std::vector<Roaring> sets;
};


#endif
