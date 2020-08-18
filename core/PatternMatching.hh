#ifndef PATTERN_MATCHING_HH
#define PATTERN_MATCHING_HH

#include <vector>
#include <execution>

#include "DataGraph.hh"

#include "utils.hh"


#define MAP_UNCHECKED(match, qv, dv) {\
  match.map(qv, dv);\
  if constexpr (L == Graph::PARTIALLY_LABELLED || L == Graph::DISCOVER_LABELS) match.map_label(qv, gpb->label(dv));\
}

#define MAPR(match, qv, dv) {\
  MAP__(match, qv, dv, vgs)\
}

#define MAP(match, qv, dv) {\
  MAP__(match, qv, dv, p)\
}

#define CHECK_LABEL(qv, dv) {\
  if constexpr (L == Graph::LABELLED)\
  {\
    if (gpb->label(dv) != p.label(qv)) continue;\
  }\
  else if constexpr (L == Graph::PARTIALLY_LABELLED)\
  {\
    if ((p.label(qv) == static_cast<uint32_t>(-1) && gpb->known_label(gpb->label(dv)))\
        || (p.label(qv) != static_cast<uint32_t>(-1) && gpb->label(dv) != p.label(qv)))\
    {\
      continue;\
    }\
  }\
}

#define MAP2(match, qv, dv) {\
  if constexpr (L == Graph::PARTIALLY_LABELLED)\
  {\
    if ((p.label(qv) == static_cast<uint32_t>(-1) && gpb->known_label(gpb->label(dv)))\
        || (p.label(qv) != static_cast<uint32_t>(-1) && gpb->label(dv) != p.label(qv)))\
    {\
      continue;\
    }\
  }\
\
  MAP_UNCHECKED(match, qv, dv)\
}

#define MAP__(match, qv, dv, p) {\
  if constexpr (L == Graph::LABELLED)\
  {\
    if (gpb->label(dv) != p.label(qv)) continue;\
  }\
  else if constexpr (L == Graph::PARTIALLY_LABELLED)\
  {\
    if ((p.label(qv) == static_cast<uint32_t>(-1) && gpb->known_label(gpb->label(dv)))\
        || (p.label(qv) != static_cast<uint32_t>(-1) && gpb->label(dv) != p.label(qv)))\
    {\
      continue;\
    }\
  }\
\
  MAP_UNCHECKED(match, qv, dv)\
}

#define PROCESS(m)\
{\
  \
  if constexpr (HAE)\
  {\
    if (!check_conditions(m.mapping)) continue;\
  }\
  \
  if constexpr (HAV)\
  {\
    if (!check_anti_vertices(m)) continue;\
  }\
  \
  if constexpr (L == Graph::DISCOVER_LABELS)\
  {\
    /* This labelling is equivalent to others, so need to canonicalize the labels.\
     * Note that the discovered label in partially labelled patterns will\
     * never create symmetry, so this is only necessary for full label discovery.\
     * XXX: assumes the pattern is a 3-star\
     */\
    uint32_t centre = rbi.vmap[vgsi][qo[0]][0] - 1;\
    uint32_t left = (centre + 1) % 3;\
    uint32_t right = (centre + 2) % 3;\
    if (gpb->label(m.mapping[left]) < gpb->label(m.mapping[right]))\
    {\
      std::swap(m.labels[left], m.labels[right]);\
      std::swap(m.mapping[left], m.mapping[right]);\
      user_process({m.mapping, m.labels});\
      /* undo swap: matcher is still using m!*/\
      std::swap(m.labels[left], m.labels[right]);\
      std::swap(m.mapping[left], m.mapping[right]);\
    }\
    else\
    {\
      user_process({m.mapping, m.labels});\
    }\
  }\
  else if constexpr (L == Graph::UNLABELLED)\
  {\
    user_process(m.mapping);\
  }\
  else if constexpr (L == Graph::LABELLED)\
  {\
    user_process({m.mapping, p.get_labels()});\
  }\
  else if constexpr (L == Graph::PARTIALLY_LABELLED)\
  {\
    user_process({m.mapping, m.labels});\
  }\
}

namespace Peregrine
{
  uint64_t binom(uint32_t n, uint32_t k)
  {
    uint64_t res = 1;
    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
    {
        k = n - k;
    }

    // Calculate value of
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (uint32_t i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
  }

  // XXX: this is a weak alias. is there any way to make this a strong typedef
  // without a performance degradation?
  using Pattern = std::vector<uint32_t>;
  Pattern unused_labels;

  struct CompleteMatch
  {
    CompleteMatch(const std::vector<uint32_t> &m)
      : mapping(m), pattern(unused_labels) {}

    CompleteMatch(const std::vector<uint32_t> &m, const Pattern &labels)
      : mapping(m), pattern(labels) {}

    const std::vector<uint32_t> &mapping;
    const Pattern &pattern;
  };

  template <Graph::Labelling L>
  struct partial_match
  {
    partial_match(size_t n) : mapping(n, 0), labels(n, 0) {}

    uint32_t at(uint32_t u) const {
      return mapping[u - 1];
    }

    void map(uint32_t u, uint32_t v, uint32_t l) {
      mapping[u-1] = v;
      labels[u-1] = l;
    }

    void map_label(uint32_t u, uint32_t l)
    {
      labels[u-1] = l;
    }

    void map(uint32_t u, uint32_t v) {
      mapping[u-1] = v;
    }

    uint32_t label(uint32_t u) const
    {
      return labels[u-1];
    }

    void unmap(uint32_t u) {
      mapping[u-1] = 0;
    }

    size_t max_size() const {
      return mapping.size();
    }

    size_t size() const {
      size_t c = 0;
      for (const auto &v : mapping) {
        if (v > 0) ++c;
      }
      return c;
    }

    bool mapped(uint32_t u) const {
      return mapping[u-1] != 0;
    }

    bool data_mapped(uint32_t v) const {
      for (const auto &m : mapping) {
        if (m == v) return true;
      }
      return false;
    }

    std::vector<uint32_t> mapping;
    std::vector<uint32_t> labels;
  };

  template <>
  struct partial_match<Graph::LABELLED>
  {
    partial_match(size_t n) : mapping(n, 0) {}
    partial_match(const partial_match &other) : mapping(other.mapping) {}
    partial_match(const partial_match &other, uint32_t qv, uint32_t dv) : mapping(other.mapping) {
      mapping[qv - 1] = dv;
    }

    uint32_t at(uint32_t u) const {
      return mapping[u - 1];
    }

    void map(uint32_t u, uint32_t v) {
      mapping[u-1] = v;
    }

    void map_label(uint32_t , uint32_t)
    {
      assert(false);
    }

    uint32_t label(uint32_t) const
    {
      assert(false);
      return 0;
    }

    void unmap(uint32_t u) {
      mapping[u-1] = 0;
    }

    size_t max_size() const {
      return mapping.size();
    }

    size_t size() const {
      size_t c = 0;
      for (const auto &v : mapping) {
        if (v > 0) ++c;
      }
      return c;
    }

    bool mapped(uint32_t u) const {
      return mapping[u-1] != 0;
    }

    bool data_mapped(uint32_t v) const {
      for (const auto &m : mapping) {
        if (m == v) return true;
      }
      return false;
    }

    std::vector<uint32_t> mapping;
  };

  template <>
  struct partial_match<Graph::UNLABELLED>
  {
    partial_match(size_t n) : mapping(n, 0) {}
    partial_match(const partial_match &other) : mapping(other.mapping) {}
    partial_match(const partial_match &other, uint32_t qv, uint32_t dv) : mapping(other.mapping) {
      mapping[qv - 1] = dv;
    }

    uint32_t at(uint32_t u) const {
      return mapping[u - 1];
    }

    void map(uint32_t u, uint32_t v) {
      mapping[u-1] = v;
    }

    void map_label(uint32_t , uint32_t)
    {
      assert(false);
    }

    uint32_t label(uint32_t) const
    {
      assert(false);
      return 0;
    }

    void unmap(uint32_t u) {
      mapping[u-1] = 0;
    }

    size_t max_size() const {
      return mapping.size();
    }

    size_t size() const {
      size_t c = 0;
      for (const auto &v : mapping) {
        if (v > 0) ++c;
      }
      return c;
    }

    bool mapped(uint32_t u) const {
      return mapping[u-1] != 0;
    }

    bool data_mapped(uint32_t v) const {
      for (const auto &m : mapping) {
        if (m == v) return true;
      }
      return false;
    }

    std::vector<uint32_t> mapping;
  };



  template <bool HAV, StoppableOption stoppable, typename Func>
  struct Matcher
  {
    Matcher(const AnalyzedPattern &r, const DataGraph * const dg, uint8_t vgsi, std::vector<std::vector<uint32_t>> &cands, const Func &f)
      : rbi(r), vgsi(vgsi), vgs(rbi.vgs[vgsi]), gpb(dg), qo(rbi.qo_book[vgsi]), cands(cands), user_process(f) {}

    template <bool HAE, Graph::Labelling L>
    void nonCoreVertexMatching(const partial_match<L> &mvgs)
    {
      const SmallGraph &p = rbi.query_graph;
      const auto &qs_list = rbi.qs[vgsi];
      for (uint32_t qsi = 0; qsi < qs_list.size(); ++qsi) {
        const auto &qs = qs_list[qsi];
        partial_match<L> m(p.num_vertices());
        for (unsigned j = 0; j < qs.size(); ++j)
        {
          // update map to be consistent with query graph
          uint32_t qv = rbi.vmap[vgsi][qs[j]][qsi];
          m.map(qv, mvgs.at(qs[j]));
          if constexpr (L != Graph::UNLABELLED && L != Graph::LABELLED) m.map_label(qv, mvgs.label(qs[j]));
        }

        complete_match<L, HAE>(m);
      }
    }

    bool check_conditions(const std::vector<uint32_t> &mapping)
    {
      const auto &conditions = rbi.conditions;
      for (const auto &[l, g] : conditions)
      {
        if (mapping[l-1] >= mapping[g-1]) return false;
      }

      return true;
    }

    template <Graph::Labelling L>
    bool check_anti_vertices(const partial_match<L> &m)
    {
      const auto &p = rbi.query_graph;
      for (auto av : rbi.anti_vertices)
      {
        const std::vector<uint32_t> &anti_parents = rbi.query_graph.get_anti_neighbours(av);
        std::vector<uint32_t> &candidates = *(cands.end()-2);
        std::vector<uint32_t> &t = cands.back();
        const adjlist &fadj = gpb->get_adj(m.at(anti_parents.front()));
        candidates.assign(fadj.ptr, fadj.ptr + fadj.length);
        for (uint32_t i = 1; i < anti_parents.size(); ++i)
        {
          const uint32_t du = m.at(anti_parents[i]);
          const adjlist &duadj = gpb->get_adj(du);
          t.clear();
          t.swap(candidates);
          std::set_intersection(std::execution::unseq,
              t.cbegin(), t.cend(),
              duadj.ptr, duadj.ptr + duadj.length,
              std::back_inserter(candidates));
        }

        std::vector<uint32_t> reg_v = p.get_neighbours(anti_parents.front());
        for (uint32_t i = 1; i < anti_parents.size(); ++i)
        {
          const uint32_t pu = anti_parents[i];
          const auto &puadj = p.get_neighbours(pu);
          t.clear();
          t.swap(reg_v);
          std::set_intersection(std::execution::unseq,
              t.cbegin(), t.cend(),
              puadj.cbegin(), puadj.cend(),
              std::back_inserter(reg_v));
        }

        if (!candidates.empty())
        {
          if constexpr (L == Graph::UNLABELLED)
          {
            // if pattern is unlabelled, then we don't care about anti-vertex labels
            if (anti_parents.size() == p.num_vertices()) return false;

            // can never meet constraint
            if (candidates.size() > reg_v.size()) return false;

            // candidates.size() <= reg_v.size(), so this is cheap
            for (const uint32_t pv : reg_v)
            {
              if (!utils::search(candidates, m.at(pv))) return false;
            }
          }
          else
          {
            // if pattern is labelled or we're discovering labels, then an
            // unlabelled anti-vertex is the same as in the unlabelled pattern case
            if (p.label(av) == static_cast<uint32_t>(-1))
            {
              if (anti_parents.size() == p.num_vertices()) return false;

              // can never meet constraint
              if (candidates.size() > reg_v.size()) return false;

              // candidates.size() <= reg_v.size(), so this is cheap
              for (const uint32_t pv : reg_v)
              {
                if (!utils::search(candidates, m.at(pv))) return false;
              }
            }
            else
            {
              // otherwise, we just need no element with its label in the candidate set
              uint32_t undesirlable = p.label(av); // I think I'm funny

              // count number of pattern vertices with same label
              uint32_t c = std::count(std::execution::unseq,
                  p.get_labels().cbegin(), p.get_labels().cend(),
                  undesirlable);

              // don't count anti-vertices with the undesirlable
              for (uint32_t av : rbi.anti_vertices)
              {
                if (p.label(av) == undesirlable) c -= 1;
              }

              // don't count anti-parents with the undesirlable, as they are
              // not part of the candidates
              for (uint32_t pv : anti_parents)
              {
                if (p.label(pv) == undesirlable) c -= 1;
              }

              if (c < candidates.size())
              {
                // candidates size can be bigger, so long as there are only c undesirlables
                uint32_t check = 0;
                for (const uint32_t dv : candidates)
                {
                  if (gpb->label(dv) == undesirlable) check += 1;
                  if (check > c) return false;
                }
              }
            }
          }
        }
      }

      return true;
    }

    /**
     * Sibling groups share the same candidate set.
     * Order groups are totally ordered subsets of sibling groups.
     */
    template <Graph::Labelling L, bool HAE>
    void complete_match(partial_match<L> &m)
    {
      const uint8_t vgs_sz = vgs.num_vertices();
      const auto &p = rbi.query_graph;

      const std::vector<uint32_t> &candidate_idxs = rbi.candidate_idxs;
      const std::vector<std::vector<uint32_t>> &ordgs = rbi.order_groups;
      std::vector<std::pair<uint32_t, uint32_t>> views(ordgs.size());
      const std::vector<uint32_t> &sibgroups = rbi.sibling_groups;

      { // generate candidate sets

        std::vector<uint32_t> &t = cands.back();

        // generate candidate sets for sibling groups
        uint8_t idx = 0;
        for (const uint32_t g : sibgroups)
        {
          if constexpr (stoppable == STOPPABLE) pthread_testcancel();

          const auto &true_parents = p.get_neighbours(g);
          const auto &anti_parents = p.get_anti_neighbours(g);
          auto &candidates = cands[vgs_sz + idx];
          const adjlist &a = gpb->get_adj(m.at(true_parents.front()));
          uint32_t *start = a.ptr;
          uint32_t *end = a.ptr + a.length;

          // views for sibling groups
          const auto &lower = rbi.lower_bounds[g];
          const auto &upper = rbi.upper_bounds[g];
          if (!lower.empty())
          {
            uint32_t lower_bound = *std::max_element(lower.begin(), lower.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            start = std::upper_bound(start, end, m.at(lower_bound));
          }

          if (!upper.empty())
          {
            bool check_upper_bound = true;
            if constexpr (HAE)
            {
              check_upper_bound = rbi.check_sibg_bound[idx];
            }

            if (check_upper_bound)
            {
              uint32_t upper_bound = *std::min_element(upper.begin(), upper.end(),
                  [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
              end = std::lower_bound(start, end, m.at(upper_bound));
            }
          }

          candidates.assign(start, end);

          for (uint32_t i = 1; i < true_parents.size(); ++i)
          {
            if constexpr (stoppable == STOPPABLE) pthread_testcancel();
            const uint32_t du = m.at(true_parents[i]);
            const adjlist &duadj = gpb->get_adj(du);
            t.clear();
            t.swap(candidates);
            std::set_intersection(std::execution::unseq,
                t.cbegin(), t.cend(),
                duadj.ptr, duadj.ptr + duadj.length,
                std::back_inserter(candidates));
          }

          // can't guard this with HAE since vertex induced doesn't count as HAE
          {
            for (uint32_t i = 0; i < anti_parents.size(); ++i)
            {
              if constexpr (stoppable == STOPPABLE) pthread_testcancel();
              if (p.is_anti_vertex(anti_parents[i])) continue;
              const uint32_t du = m.at(anti_parents[i]);
              const adjlist &duadj = gpb->get_adj(du);
              t.clear();
              t.swap(candidates);
              std::set_difference(std::execution::unseq,
                  t.cbegin(), t.cend(),
                  duadj.ptr, duadj.ptr + duadj.length,
                  std::back_inserter(candidates));
            }
          }

          idx += 1;
        }
      }

      { // generate views
        if constexpr (stoppable == STOPPABLE) pthread_testcancel();
        uint8_t oidx = 0;
        for (const auto &g : ordgs)
        {
          // can't mutate this! it's shared by different order groups!
          const auto &candidates = cands[candidate_idxs[oidx]];
          const auto &lower = gpb->get_lower_bounds(g.front());
          const auto &upper = gpb->get_upper_bounds(g.back());
          uint32_t start = 0;
          uint32_t end = candidates.size();

          auto cand_end = candidates.end();

          if (!lower.empty())
          {
            uint32_t lower_bound = *std::max_element(lower.begin(), lower.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            start = std::distance(candidates.begin(),
                std::upper_bound(candidates.begin(), cand_end,
                  m.at(lower_bound)));
          }

          if (!upper.empty())
          {
            uint32_t upper_bound = *std::min_element(upper.begin(), upper.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            end = std::distance(candidates.begin(),
                std::lower_bound(candidates.begin() + start, cand_end,
                  m.at(upper_bound)));
          }

          views[oidx] = {start, end};

          oidx += 1;
        }
      }

      { // matching
        if constexpr (stoppable == STOPPABLE) pthread_testcancel();
        const auto &unmatched = rbi.ncore_vertices;
        if (unmatched.size() == 1)
        {
          // special case for one order group with one noncore vertex (very common)
          uint32_t qv = unmatched.front();
          auto &candidates = cands[vgs_sz];

          auto start = candidates.begin();
          auto end = candidates.end();
          if constexpr (L == Graph::LABELLED)
          {
            uint32_t l = p.label(qv);
            end = std::remove_if(std::execution::unseq,
                candidates.begin(),
                candidates.end(),
                [this, &l](uint32_t dv) { return gpb->label(dv) != l; });
          }

          auto range_start = start;
          auto range_end = end;

          std::vector<uint32_t> sorted_mapping(m.mapping);
          std::sort(sorted_mapping.begin(), sorted_mapping.end());
          for (uint32_t matched : sorted_mapping)
          {
            auto it = std::lower_bound(start, end, matched);

            // new range
            if (it != end && *it == matched) range_end = it;
            else continue;

            for (auto vit = range_start; vit < range_end; ++vit)
            {
              const uint32_t v = *vit;
              MAP_UNCHECKED(m, qv, v);
              PROCESS(m);
            }
            m.unmap(qv);
            range_start = range_end + 1;
          }

          for (auto vit = range_start; vit < end; ++vit)
          {
            const uint32_t v = *vit;
            MAP_UNCHECKED(m, qv, v);
            PROCESS(m);
          }
        }
        else if (ordgs.size() == 1)
        {
          // special case for one order group with multiple noncore vertices
          // (fairly common).  note that this will never be partially labelled:
          // unmatched.size() > 1, and partially labelled vertices get their
          // own order group because they have a unique label
          auto &candidates = cands[vgs_sz];
          auto [start, end] = views.front();

          if constexpr (L == Graph::LABELLED)
          {
            uint32_t l = p.label(unmatched.front());
            // remove labels not present in this sibling group
            auto e = std::remove_if(std::execution::unseq,
                candidates.begin(),
                candidates.end(),
                [this, &m, &l](uint32_t dv) { return gpb->label(dv) != l || m.data_mapped(dv); });
            end = std::distance(candidates.begin(), e);
          }

          if constexpr (HAE)
          {
            auto e = std::remove_if(std::execution::unseq,
                candidates.begin(),
                candidates.end(),
                [this, &m](uint32_t dv) { return m.data_mapped(dv); });
            end = std::distance(candidates.begin(), e);
          }

          int8_t sz = p.num_vertices() - vgs_sz;
          // last one is unused, but required for safety (when ascending from middle to last level)
          std::vector<uint32_t> cursors(sz);
          for (uint32_t i = 0; i < cursors.size(); ++i)
          {
            cursors[i] = start + i;
          }

          uint8_t idx = 0;
          while (true)
          {
            if constexpr (stoppable == STOPPABLE) pthread_testcancel();
            uint32_t qv = unmatched[idx];
            if (idx == sz-1)
            {
              for (uint32_t i = cursors[idx-1]; i < end; ++i)
              {
                uint32_t dv = candidates[i];
                MAP_UNCHECKED(m, qv, dv);
                PROCESS(m);
              }
              idx -= 1;
            }
            else
            {
              uint32_t &cursor = cursors[idx];
              if (cursor < end)
              {
                uint32_t dv = candidates[cursor++];
                MAP_UNCHECKED(m, qv, dv);
                cursors[idx+1] = cursor;
                idx += 1;
              }
              else if (idx == 0)
              {
                break;
              }
              else
              {
                cursor = cursors[idx-1] + 1;
                m.unmap(qv);
                idx -= 1;
              }
            }
          }
        }
        else if (unmatched.size() == 2 && sibgroups.size() == 2)
        {
          assert(unmatched[0] == sibgroups[0]);

          std::vector<uint32_t> &cands1 = cands[vgs_sz];
          std::vector<uint32_t> &cands2 = cands[vgs_sz + 1];

          if constexpr (L == Graph::LABELLED)
          {
            uint32_t l = p.label(sibgroups.front());
            // remove labels not present in this sibling group
            cands1.erase(std::remove_if(std::execution::unseq,
                cands1.begin(),
                cands1.end(),
                [this, &l](uint32_t dv) { return gpb->label(dv) != l; }),
                cands1.end());

            l = p.label(sibgroups.back());
            // remove labels not present in this sibling group
            cands2.erase(std::remove_if(std::execution::unseq,
                cands2.begin(),
                cands2.end(),
                [this, &l](uint32_t dv) { return gpb->label(dv) != l; }),
                cands2.end());
          }
          else if constexpr (L == Graph::PARTIALLY_LABELLED)
          {
            // if partially labelled, know two things:
            //  - the partially labelled vertex's candidates will not contain
            //    any of the core vertices
            //      (because those have known labels)
            //  - the intersection with the labelled vertex will be empty
            //      (because it has a known label)

            // to make things easier, have cands1 always point to candidates for labelled vertex,
            // and cands2 point to candidates for -1 labelled vertex
            uint32_t l = p.label(sibgroups.front());
            if (l == static_cast<uint32_t>(-1))
            {
              cands1 = cands[vgs_sz + 1];
              cands2 = cands[vgs_sz];
              l = p.label(sibgroups.back());
            }

            // remove labels not present in this sibling group
            cands1.erase(std::remove_if(std::execution::unseq,
                cands1.begin(),
                cands1.end(),
                [this, &l](uint32_t dv) { return gpb->label(dv) != l; }),
                cands1.end());

            // remove known labels
            cands2.erase(std::remove_if(std::execution::unseq,
                cands2.begin(),
                cands2.end(),
                [this](uint32_t dv) { return gpb->known_label(gpb->label(dv)); }),
                cands2.end());
          }

          std::vector<uint32_t> &sintersection = cands.back();
          sintersection.clear();

          auto rem = [&]() {
            // remove core vertices from the candidate sets
            // note that this messes up the sorting on cands
            std::vector<std::vector<uint32_t>::iterator> to_delete1;
            std::vector<std::vector<uint32_t>::iterator> to_delete2;
            std::vector<std::vector<uint32_t>::iterator> to_delete3;

            std::vector<uint32_t> sorted_mapping(m.mapping);
            std::sort(sorted_mapping.begin(), sorted_mapping.end(), std::greater<uint32_t>());
            for (uint32_t matched : sorted_mapping)
            {
              // swap to end and resize
              auto it = std::lower_bound(cands1.begin(), cands1.end(), matched);
              if (it != cands1.end() && (*it) == matched)
              {
                to_delete1.push_back(it);
              }

              if constexpr (L != Graph::PARTIALLY_LABELLED)
              {
                it = std::lower_bound(cands2.begin(), cands2.end(), matched);
                if (it != cands2.end() && (*it) == matched)
                {
                  to_delete2.push_back(it);
                }

                // remove from intersection as well
                it = std::lower_bound(cands.back().begin(), cands.back().end(), matched);
                if (it != cands.back().end() && (*it) == matched)
                {
                  to_delete3.push_back(it);
                }
              }
            }

            uint8_t i = 0;
            for (auto it : to_delete1)
            {
              std::iter_swap(it, cands1.end() - ++i);
            }
            std::vector<uint32_t>::iterator del_cur = cands1.end() - to_delete1.size();
            cands1.erase(del_cur, cands1.end());


            if constexpr (L != Graph::PARTIALLY_LABELLED)
            {
              i = 0;
              for (auto it : to_delete2)
              {
                std::iter_swap(it, cands2.end() - ++i);
              }

              i = 0;
              for (auto it : to_delete3)
              {
                std::iter_swap(it, cands.back().end() - ++i);
              }

              del_cur = cands2.end() - to_delete2.size();
              cands2.erase(del_cur, cands2.end());

              del_cur = cands.back().end() - to_delete3.size();
              cands.back().erase(del_cur, cands.back().end());
            }
          };

          if constexpr (L != Graph::PARTIALLY_LABELLED)
          {
            std::set_intersection(std::execution::unseq,
                cands1.cbegin(), cands1.cend(),
                cands2.cbegin(), cands2.cend(),
                std::back_inserter(sintersection));
          }

          if (sintersection.empty())
          {
            rem();
            for (uint32_t u : cands1)
            {
              MAP_UNCHECKED(m, unmatched[0], u);
              for (uint32_t v : cands2)
              {
                MAP_UNCHECKED(m, unmatched[1], v);
                PROCESS(m);
              }
            }
          }
          else
          {
            // this means we are not partially labelled!

            std::vector<uint32_t> t;

            t.swap(cands1);
            std::set_difference(std::execution::unseq,
                t.cbegin(), t.cend(),
                sintersection.cbegin(), sintersection.cend(),
                std::back_inserter(cands1));

            t.clear();

            t.swap(cands2);
            std::set_difference(std::execution::unseq,
                t.cbegin(), t.cend(),
                sintersection.cbegin(), sintersection.cend(),
                std::back_inserter(cands2));

            rem();

            for (uint32_t u : cands1)
            {
              MAP_UNCHECKED(m, unmatched[0], u);
              for (uint32_t v : cands2)
              {
                MAP_UNCHECKED(m, unmatched[1], v);
                PROCESS(m);
              }

              // sintersection mapped to second noncore
              for (uint32_t v : sintersection)
              {
                MAP_UNCHECKED(m, unmatched[1], v);
                PROCESS(m);
              }
            }

            // sintersection mapped to first noncore
            for (uint32_t u : sintersection)
            {
              MAP_UNCHECKED(m, unmatched[0], u);
              for (uint32_t v : cands2)
              {
                MAP_UNCHECKED(m, unmatched[1], v);
                PROCESS(m);
              }
            }

            // this part requires checking data_mapped, so unmap first
            m.unmap(unmatched[0]);
            m.unmap(unmatched[1]);

            // sintersection mapped to both
            for (uint32_t iu = 0; iu < sintersection.size(); ++iu)
            {
              uint32_t u = sintersection[iu];
              MAP_UNCHECKED(m, unmatched[0], u);

              for (uint32_t iv = iu + 1; iv < sintersection.size(); ++iv)
              {
                uint32_t v = sintersection[iv];
                MAP_UNCHECKED(m, unmatched[1], v);

                PROCESS(m);

                MAP_UNCHECKED(m, unmatched[0], v);
                MAP_UNCHECKED(m, unmatched[1], u);

                if constexpr (HAE)
                {
                  if (check_conditions(m.mapping))
                  {
                    PROCESS(m); // XXX check_conditions called twice
                  }
                }
                else
                {
                  PROCESS(m);
                }

                // remap u (always unchecked)
                MAP_UNCHECKED(m, unmatched[0], u);
              }
            }
          }
        }
        else
        {
          {
            int8_t sz = p.num_vertices() - vgs_sz;
            std::vector<uint32_t> cursors(sz);
            std::vector<uint32_t> oidxs(sz);
            std::vector<uint32_t> flat_ordgs(sz);
            uint8_t ii = 0;
            for (uint8_t oidx = 0; oidx < ordgs.size(); ++oidx)
            {
              for (uint32_t jj = 0; jj < ordgs[oidx].size(); ++jj)
              {
                cursors[ii] = views[oidx].first + jj;
                oidxs[ii] = oidx;
                // XXX: no reason not to do this in advance
                flat_ordgs[ii] = ordgs[oidx][jj];
                ii += 1;
              }
            }

            uint8_t idx = 0;
            while (true)
            {
              uint8_t oidx = oidxs[idx];
              const auto [start, end] = views[oidx];
              uint32_t qv = flat_ordgs[idx];
              uint32_t cursor = cursors[idx];
              const auto &candidates = cands[candidate_idxs[oidx]];

              auto cend = candidates.cbegin() + end;

              if (idx == sz-1)
              {

                auto range_start = candidates.cbegin() + cursor;
                auto range_end = cend;

                // Important to use m.mapping (instead of mvgs.mapping like in the unmatched.size() == 1 case)
                // since we may have mapped other things outside of this order group that may be in the candidate set.
                // Think of chains for example.
                std::vector<uint32_t> sorted_mapping(m.mapping);
                std::sort(std::execution::unseq, sorted_mapping.begin(), sorted_mapping.end());
                for (uint32_t matched : sorted_mapping)
                {
                  auto it = std::lower_bound(range_start, cend, matched);

                  // new range
                  if (it != cend && *it == matched)
                  {
                    range_end = it;
                    for (auto vit = range_start; vit < range_end; ++vit)
                    {
                      const uint32_t v = *vit;
                      CHECK_LABEL(qv, v);
                      MAP_UNCHECKED(m, qv, v);
                      m.map(qv, v);
                      PROCESS(m);
                    }
                    m.unmap(qv);
                    // skip what you're supposed to skip
                    range_start = range_end + 1;
                  }
                }

                for (auto vit = range_start; vit < cend; ++vit)
                {
                  const uint32_t v = *vit;
                  CHECK_LABEL(qv, v);
                  MAP_UNCHECKED(m, qv, v);
                  PROCESS(m);
                };

                cursors[idx] = start;
                m.unmap(qv);
                if (idx == 0) break;
                idx -= 1;
              }
              else
              {
                if (cursors[idx] < end)
                {
                  uint32_t dv = candidates[cursor];
                  cursors[idx] += 1;

                  CHECK_LABEL(qv, dv);

                  // this can't come out of the loop, since you don't know what you mapped in previous iterations
                  if (!m.data_mapped(dv))
                  {
                    MAP_UNCHECKED(m, qv, dv);
                    // in the same order group, can't have different labels or anything
                    // so we won't miss anything by setting cursors[idx+1] to cursor+1
                    if (oidxs[idx+1] == oidx) cursors[idx+1] = cursor + 1;
                    idx += 1;
                  }
                }
                else if (idx == 0)
                {
                  break;
                }
                else
                {
                  cursors[idx] = start;
                  m.unmap(qv);
                  idx -= 1;
                }
              }
            }
          }
        }
      }
    }

    template <Graph::Labelling L, bool HAE>
    void get_next_cand(std::vector<uint32_t> &xsection, unsigned idx, const partial_match<L> &mvgs)
    {
      // query vertices already matched, and connected from qo[key]
      std::vector<uint32_t> Ucon;
      const auto &true_edges = vgs.get_neighbours(qo[idx]);
      // all vertices in the range [levels, key) have been matched
      // want all such vertices that are adjacent to qo[key]
      for (const uint32_t m : true_edges) {
        if constexpr (stoppable == STOPPABLE) pthread_testcancel();
        for (uint32_t qi = 0; qi < idx; ++qi) {
          if (m == qo[qi]) {
            assert(mvgs.mapped(m));
            Ucon.push_back(m);
          }
        }
      }

      const adjlist &adj = gpb->get_adj(mvgs.at(Ucon[0]));

      uint32_t *start = adj.ptr;
      uint32_t *end = adj.ptr + adj.length;

      uint32_t lower_bound = 0;
      uint32_t upper_bound = static_cast<uint32_t>(-1); // largest possible
      for (uint32_t i = 0; i < idx; ++i)
      {
        if constexpr (stoppable == STOPPABLE) pthread_testcancel();
        if (qo[i] < qo[idx]) lower_bound = std::max(lower_bound, qo[i]);
        else if (qo[i] > qo[idx]) upper_bound = std::min(upper_bound, qo[i]);
      }

      // enforce total ordering on core vertices by only considering correct
      // range of candidate data vertices

      start = lower_bound > 0
        ? std::upper_bound(adj.ptr, adj.ptr + adj.length, mvgs.at(lower_bound))
        : start;

      end = upper_bound < (uint32_t)-1
        ? std::lower_bound(start, adj.ptr + adj.length, mvgs.at(upper_bound))
        : end;

      xsection.assign(start, end);

      for (unsigned i = 1; i < Ucon.size(); ++i)
      {
        if constexpr (stoppable == STOPPABLE) pthread_testcancel();

        std::vector<uint32_t> t;
        t.swap(xsection);

        const uint32_t n = Ucon[i];
        const uint32_t v = mvgs.at(n);
        const adjlist &adj = gpb->get_adj(v);

        std::set_intersection(std::execution::unseq,
            t.cbegin(), t.cend(),
            adj.ptr, adj.ptr + adj.length,
            std::back_inserter(xsection));
      }

      // remove bad labels
      if constexpr (L == Graph::LABELLED)
      {
        uint32_t label = vgs.label(qo[idx]);
        xsection.erase(std::remove_if(std::execution::unseq,
              xsection.begin(),
              xsection.end(),
              [this, &label](uint32_t dv) { return gpb->label(dv) != label; }),
            xsection.end());
      }


      // can't guard this with HAE since vertex induced doesn't count as HAE
      {
        std::vector<uint32_t> nUcon;
        const auto &not_edges = vgs.get_anti_neighbours(qo[idx]);
        // all vertices in the range [levels, key) have been matched
        // want all such vertices that are adjacent to qo[key]
        for (const uint32_t m : not_edges) {
          if constexpr (stoppable == STOPPABLE) pthread_testcancel();
          for (uint32_t qi = 0; qi < idx; ++qi) {
            if (m == qo[qi]) {
              assert(mvgs.mapped(m));
              nUcon.push_back(m);
            }
          }
        }

        for (const uint32_t n : nUcon)
        {
          if constexpr (stoppable == STOPPABLE) pthread_testcancel();

          std::vector<uint32_t> t;
          t.swap(xsection);

          const uint32_t v = mvgs.at(n);
          const adjlist &adj = gpb->get_adj(v);

          std::set_difference(std::execution::unseq,
              t.cbegin(), t.cend(),
              adj.ptr, adj.ptr + adj.length,
              std::back_inserter(xsection));
        }
      }
    }

    void discover_labels(partial_match<Graph::DISCOVER_LABELS> &m)
    {
      // XXX: 3-star only
      uint32_t centre = rbi.vmap[0][1][0] - 1;
      uint32_t left = (centre + 1) % 3;
      uint32_t right = (centre + 2) % 3;
      const adjlist &cands = gpb->get_adj(m.mapping[centre]);

      for (uint32_t i = 0; i < cands.length; ++i)
      {
        uint32_t u = cands.ptr[i];
        uint32_t ulab = gpb->label(u);
        for (uint32_t j = i + 1; j < cands.length; ++j)
        {
          uint32_t v = cands.ptr[j];
          uint32_t vlab = gpb->label(v);
          if (ulab <= vlab)
          {
            m.mapping[left] = u;
            m.labels[left] = ulab;

            m.mapping[right] = v;
            m.labels[right] = vlab;
          }
          else
          {
            m.mapping[left] = v;
            m.labels[left] = vlab;

            m.mapping[right] = u;
            m.labels[right] = ulab;
          }

          user_process({m.mapping, m.labels});
        }
      }
    }


    template <Graph::Labelling L, bool HAE>
    void star_map_into(uint32_t vertex_id)
    {
      partial_match<L> m(rbi.query_graph.num_vertices());
      if constexpr (L == Graph::UNLABELLED || L == Graph::LABELLED)
      {
        m.map(rbi.vmap[0][1][0], vertex_id);
      }
      else
      {
        m.map(rbi.vmap[0][1][0], vertex_id, gpb->label(vertex_id));
      }


      if constexpr (L == Graph::DISCOVER_LABELS)
      {
        discover_labels(m);
        return;
      }
      else
      {
        complete_match<L, HAE>(m);
      }
    }

    template <Graph::Labelling L, bool HAE = false>
    void map_into(uint32_t vertex_id)
    {
      const uint8_t vgs_size = vgs.num_vertices();
      if constexpr (L == Graph::LABELLED)
      {
        if (gpb->label(vertex_id) != vgs.label(qo[0])) return;
      }
      else if constexpr (L == Graph::PARTIALLY_LABELLED)
      {
        if ((vgs.label(qo[0]) == static_cast<uint32_t>(-1) && gpb->known_label(gpb->label(vertex_id)))
            || (vgs.label(qo[0]) != static_cast<uint32_t>(-1) && gpb->label(vertex_id) != vgs.label(qo[0])))
        {
          return;
        }
      }

      if (vgs_size == 1) return star_map_into<L, HAE>(vertex_id);

      partial_match<L> mvgs(vgs_size);
      if constexpr (L == Graph::UNLABELLED || L == Graph::LABELLED)
      {
        mvgs.map(qo[0], vertex_id);
      }
      else
      {
        mvgs.map(qo[0], vertex_id, gpb->label(vertex_id));
      }

      uint8_t idx = 1;
      uint32_t *cursors = new uint32_t[vgs_size-1](); // last cursor unused?

      get_next_cand<L, HAE>(cands[idx], idx, mvgs);

      // don't need to worry about duplicates (i.e. no !data_mapped() checks)
      // because the candidates we generate for core vertices are totally
      // ordered
      if constexpr (stoppable == STOPPABLE) pthread_testcancel();
      if (idx == vgs_size-1)
      {
        const auto &candidates = cands[idx];
        const uint32_t key = qo[idx];
        for (const uint32_t dv : candidates)
        {
          MAPR(mvgs, key, dv);

          nonCoreVertexMatching<HAE>(mvgs);
        }
      }
      else
      {
        while (idx > 0)
        {
          if constexpr (stoppable == STOPPABLE) pthread_testcancel();
          const auto &candidates = cands[idx];

          const uint32_t key = qo[idx];
          if (idx == vgs_size-1)
          {
            for (const uint32_t dv : candidates)
            {
              MAPR(mvgs, key, dv);
              nonCoreVertexMatching<HAE>(mvgs);
            }
            mvgs.unmap(key);
            idx -= 1;
          }
          else
          {
            uint32_t &cursor = cursors[idx];
            if (cursor < candidates.size())
            {
              uint32_t dv = candidates[cursor++];
              // XXX: do CHECK_LABEL optimization
              // prepare next level candidates
              MAPR(mvgs, key, dv);
              idx += 1;
              get_next_cand<L, HAE>(cands[idx], idx, mvgs);
            }
            else
            {
              cursor = 0;
              mvgs.unmap(key);
              idx -= 1;
            }
          }
        }
      }

      delete[] cursors;
    }

    const AnalyzedPattern &rbi;
    unsigned vgsi;
    const SmallGraph &vgs;
    const DataGraph * const gpb;
    std::vector<uint32_t> qo;
    std::vector<std::vector<uint32_t>> &cands;
    const Func &user_process;
  };


  template <bool HAV>
  struct Counter
  {
    Counter(const AnalyzedPattern &r, const DataGraph * const dg, uint8_t vgsi, std::vector<std::vector<uint32_t>> &cands)
      : rbi(r), vgsi(vgsi), vgs(rbi.vgs[vgsi]), gpb(dg), qo(rbi.qo_book[vgsi]), cands(cands)
    {
      assert(!r.has_anti_edges()); // Counter doesn't support it!
    }

    template <Graph::Labelling L>
    uint64_t nonCoreVertexMatching(const partial_match<L> &mvgs)
    {
      const SmallGraph &p = rbi.query_graph;
      const auto &qs_list = rbi.qs[vgsi];
      uint64_t c = 0;
      for (uint32_t qsi = 0; qsi < qs_list.size(); ++qsi)
      {
        const auto &qs = qs_list[qsi];
        partial_match<L> m(p.num_vertices());
        for (unsigned j = 0; j < qs.size(); ++j)
        {
          // update map to be consistent with query graph
          uint32_t qv = rbi.vmap[vgsi][qs[j]][qsi];
          m.map(qv, mvgs.at(qs[j]));
          if constexpr (L != Graph::UNLABELLED && L != Graph::LABELLED)
          {
            m.map_label(qv, mvgs.label(qs[j]));
          }
        }

        c += complete_match<L>(m);
      }

      return c;
    }

    template <Graph::Labelling L>
    bool check_anti_vertices(const partial_match<L> &m)
    {
      const auto &p = rbi.query_graph;
      for (auto av : rbi.anti_vertices)
      {
        const std::vector<uint32_t> &anti_parents = rbi.query_graph.get_anti_neighbours(av);
        std::vector<uint32_t> &candidates = *(cands.end()-2);
        std::vector<uint32_t> &t = cands.back();
        const adjlist &fadj = gpb->get_adj(m.at(anti_parents.front()));
        candidates.assign(fadj.ptr, fadj.ptr + fadj.length);
        for (uint32_t i = 1; i < anti_parents.size(); ++i)
        {
          const uint32_t du = m.at(anti_parents[i]);
          const adjlist &duadj = gpb->get_adj(du);
          t.clear();
          t.swap(candidates);
          std::set_intersection(std::execution::unseq,
              t.cbegin(), t.cend(),
              duadj.ptr, duadj.ptr + duadj.length,
              std::back_inserter(candidates));
        }

        std::vector<uint32_t> reg_v = p.get_neighbours(anti_parents.front());
        for (uint32_t i = 1; i < anti_parents.size(); ++i)
        {
          const uint32_t pu = anti_parents[i];
          const auto &puadj = p.get_neighbours(pu);
          t.clear();
          t.swap(reg_v);
          std::set_intersection(std::execution::unseq,
              t.cbegin(), t.cend(),
              puadj.cbegin(), puadj.cend(),
              std::back_inserter(reg_v));
        }

        if (!candidates.empty())
        {
          if constexpr (L == Graph::UNLABELLED)
          {
            // if pattern is unlabelled, then we don't care about anti-vertex labels
            if (anti_parents.size() == p.num_vertices()) return false;

            // can never meet constraint
            if (candidates.size() > reg_v.size()) return false;

            // candidates.size() <= reg_v.size(), so this is cheap
            for (const uint32_t pv : reg_v)
            {
              if (!utils::search(candidates, m.at(pv))) return false;
            }
          }
          else
          {
            // if pattern is labelled or we're discovering labels, then an
            // unlabelled anti-vertex is the same as in the unlabelled pattern case
            if (p.label(av) == static_cast<uint32_t>(-1))
            {
              if (anti_parents.size() == p.num_vertices()) return false;

              // can never meet constraint
              if (candidates.size() > reg_v.size()) return false;

              // candidates.size() <= reg_v.size(), so this is cheap
              for (const uint32_t pv : reg_v)
              {
                if (!utils::search(candidates, m.at(pv))) return false;
              }
            }
            else
            {
              // otherwise, we just need no element with its label in the candidate set
              uint32_t undesirlable = p.label(av); // I think I'm funny

              // count number of pattern vertices with same label
              uint32_t c = std::count(std::execution::unseq,
                  p.get_labels().cbegin(), p.get_labels().cend(),
                  undesirlable);

              // don't count anti-vertices with the undesirlable
              for (uint32_t av : rbi.anti_vertices)
              {
                if (p.label(av) == undesirlable) c -= 1;
              }

              // don't count anti-parents with the undesirlable, as they are
              // not part of the candidates
              for (uint32_t pv : anti_parents)
              {
                if (p.label(pv) == undesirlable) c -= 1;
              }

              if (c < candidates.size())
              {
                // candidates size can be bigger, so long as there are only c undesirlables
                uint32_t check = 0;
                for (const uint32_t dv : candidates)
                {
                  if (gpb->label(dv) == undesirlable) check += 1;
                  if (check > c) return false;
                }
              }
            }
          }
        }
      }

      return true;
    }

    /**
     * Sibling groups share the same candidate set.
     * Order groups are totally ordered subsets of sibling groups.
     */
    template <Graph::Labelling L>
    uint64_t complete_match(partial_match<L> &m)
    {
      const uint8_t vgs_sz = vgs.num_vertices();
      const auto &p = rbi.query_graph;

      const std::vector<uint32_t> &candidate_idxs = rbi.candidate_idxs;
      const std::vector<std::vector<uint32_t>> &ordgs = rbi.order_groups;
      std::vector<std::pair<uint32_t, uint32_t>> views(ordgs.size());
      const std::vector<uint32_t> &sibgroups = rbi.sibling_groups;

      { // generate candidate sets
        std::vector<uint32_t> &t = cands.back();

        // generate candidate sets for sibling groups
        uint8_t idx = 0;
        // sibgrouops must have some true parents
        for (const uint32_t g : sibgroups)
        {
          const auto &true_parents = p.get_neighbours(g);
          const auto &anti_parents = p.get_anti_neighbours(g);
          auto &candidates = cands[vgs_sz + idx];
          const adjlist &a = gpb->get_adj(m.at(true_parents.front()));
          uint32_t *start = a.ptr;
          uint32_t *end = a.ptr + a.length;

          // views for sibling groups
          const auto &lower = rbi.lower_bounds[g];
          const auto &upper = rbi.upper_bounds[g];
          if (!lower.empty()) {
            uint32_t lower_bound = *std::max_element(lower.begin(), lower.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            start = std::upper_bound(start, end, m.at(lower_bound));
          }

          if (!upper.empty()) {
            uint32_t upper_bound = *std::min_element(upper.begin(), upper.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            end = std::lower_bound(start, end, m.at(upper_bound));
          }

          candidates.assign(start, end);

          for (uint32_t i = 1; i < true_parents.size(); ++i)
          {
            const uint32_t du = m.at(true_parents[i]);
            const adjlist &duadj = gpb->get_adj(du);
            t.clear();
            t.swap(candidates);
            std::set_intersection(std::execution::unseq,
                t.cbegin(), t.cend(),
                duadj.ptr, duadj.ptr + duadj.length,
                std::back_inserter(candidates));
          }

          for (uint32_t i = 0; i < anti_parents.size(); ++i)
          {
            if (p.is_anti_vertex(anti_parents[i])) continue;
            const uint32_t du = m.at(anti_parents[i]);
            const adjlist &duadj = gpb->get_adj(du);
            t.clear();
            t.swap(candidates);
            std::set_difference(std::execution::unseq,
                t.cbegin(), t.cend(),
                duadj.ptr, duadj.ptr + duadj.length,
                std::back_inserter(candidates));
          }

          idx += 1;
        }
      }

      { // generate views
        uint8_t oidx = 0;
        for (const auto &g : ordgs)
        {
          const auto &candidates = cands[candidate_idxs[oidx]];
          const auto &lower = gpb->get_lower_bounds(g.front());
          const auto &upper = gpb->get_upper_bounds(g.back());
          uint32_t start = 0;
          uint32_t end = candidates.size();
          if (!lower.empty())
          {
            uint32_t lower_bound = *std::max_element(lower.begin(), lower.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            start = std::distance(candidates.cbegin(),
                std::upper_bound(candidates.cbegin(), candidates.cend(),
                  m.at(lower_bound)));
          }

          if (!upper.empty())
          {
            uint32_t upper_bound = *std::min_element(upper.begin(), upper.end(),
                [&m](uint32_t a, uint32_t b) { return m.at(a) < m.at(b); });
            end = std::distance(candidates.cbegin(),
                std::lower_bound(candidates.cbegin() + start, candidates.cend(),
                  m.at(upper_bound)));
          }

          views[oidx] = {start, end};

          oidx += 1;
        }
      }

      { // matching
        const auto &unmatched = rbi.ncore_vertices;

        if (unmatched.size() == 1)
        {
          // special case for one order group with one noncore vertex (very common)
          auto &candidates = cands[vgs_sz];

          if constexpr (HAV)
          {
            if constexpr (L == Graph::LABELLED)
            {
              uint32_t qv = unmatched.front();
              uint32_t label = p.label(qv);
              auto counter = [this, qv, label, &m](uint32_t dv)
              {
                if (gpb->label(dv) == label && !m.data_mapped(dv))
                {
                  m.map(qv, dv);
                  return check_anti_vertices(m);
                }
                return false;
              };
              return std::count_if(std::execution::unseq,
                  candidates.cbegin(),
                  candidates.cend(),
                  counter);
            }
            else if constexpr (L == Graph::UNLABELLED)
            {
              uint32_t qv = unmatched.front();
              auto counter = [this, qv, &m](uint32_t dv)
              {
                if (!m.data_mapped(dv))
                {
                  m.map(qv, dv);
                  return check_anti_vertices(m);
                }
                return false;
              };
              return std::count_if(std::execution::unseq,
                  candidates.cbegin(),
                  candidates.cend(),
                  counter);
            }
            else
            {
              assert(false); // count is only for labelled/unlabelled
            }
          }
          else if constexpr (L == Graph::UNLABELLED)
          {
            uint64_t n = candidates.size();
            for (uint32_t dv : m.mapping)
            {
              n -= std::binary_search(candidates.cbegin(), candidates.cend(), dv)
                ? 1 : 0;
            }

            return n;
          }
          else if constexpr (L == Graph::LABELLED)
          {
            uint32_t label = p.label(unmatched.front());
            uint64_t n = std::count_if(std::execution::unseq,
                candidates.cbegin(),
                candidates.cend(),
                [this, &label](uint32_t dv) { return gpb->label(dv) == label; });

            for (uint32_t dv : m.mapping)
            {
              n -= std::binary_search(candidates.cbegin(), candidates.cend(), dv)
                ? 1 : 0;
            }

            return n;
          }
          else
          {
            assert(false); // count is only for labelled/unlabelled
          }
        }
        else if (ordgs.size() == 1 && !HAV)
        {
          // vertices might be matched twice: consider [1-2][1-4][2-3][3-4][2-5][4-5]
          // it's like 4-cycle but with two non-core vertices.

          uint8_t sz = unmatched.size();
          // special case for one order group with multiple noncore vertices
          // (fairly common)
          const auto &candidates = cands[vgs_sz];
          uint32_t n = candidates.size();
          if constexpr (L == Graph::LABELLED)
          {
            uint32_t label = p.label(unmatched.front());
            // one ordg => one label
            // so count how many have same label
            n = std::count_if(std::execution::unseq,
                candidates.cbegin(),
                candidates.cend(),
                [this, &label](uint32_t dv) { return gpb->label(dv) == label; });
          }

          for (uint32_t dv : m.mapping)
          {
            if (dv != 0)
            {
              n -= std::binary_search(candidates.cbegin(), candidates.cend(), dv)
                ? 1 : 0;
            }
          }

          return binom(n, sz);
        }
        else if (unmatched.size() == 2 && sibgroups.size() == 2 && !HAV && L == Graph::UNLABELLED)
        {
          // TODO: only unlabelled!

          // Just inclusion-exclusion! (can be for ordgs.size() == sibgroups.size())
          // compute intersection X of the candidate sets

          const std::vector<uint32_t> &cands1 = cands[vgs_sz];
          const std::vector<uint32_t> &cands2 = cands[vgs_sz + 1];
          std::vector<uint32_t> &sintersection = cands.back();
          sintersection.clear();
          std::set_intersection(std::execution::unseq,
              cands1.cbegin(), cands1.cend(),
              cands2.cbegin(), cands2.cend(),
              std::back_inserter(sintersection));


          uint64_t c1 = cands1.size();
          uint64_t c2 = cands2.size();
          uint64_t xs = sintersection.size();

          for (uint32_t u : m.mapping)
          {
            c1 -= std::binary_search(cands1.cbegin(), cands1.cend(), u) ? 1 : 0;
            c2 -= std::binary_search(cands2.cbegin(), cands2.cend(), u) ? 1 : 0;
            xs -= std::binary_search(sintersection.cbegin(), sintersection.cend(), u) ? 1 : 0;
          }

          sintersection.clear();
          return (c1 * c2) - xs;
        }
        else
        {
          const std::vector<uint32_t> &candidate_idxs = rbi.candidate_idxs;
          int8_t sz = unmatched.size();

          // general matcher
          std::vector<uint32_t> ends(sz);
          std::vector<uint32_t> cursors(sz);
          std::vector<uint32_t> oidxs(sz);
          std::vector<uint32_t> flat_ordgs(sz);

          uint8_t ii = 0;
          for (uint8_t oidx = 0; oidx < ordgs.size(); ++oidx)
          {
            for (uint32_t jj = 0; jj < ordgs[oidx].size(); ++jj)
            {
              cursors[ii] = views[oidx].first + jj;
              oidxs[ii] = oidx;
              flat_ordgs[ii] = ordgs[oidx][jj];
              ii += 1;
            }
          }


          uint64_t c = 0;
          uint8_t idx = 0;
          while (true)
          {
            uint8_t oidx = oidxs[idx];
            auto [start, end] = views[oidx];
            const auto &candidates = cands[candidate_idxs[oidx]];
            uint32_t qv = flat_ordgs[idx];

            auto &cursor = cursors[idx];
            if (idx == sz-1)
            {
              for (uint32_t i = cursor; i < end; ++i)
              {
                uint32_t dv = candidates[i];
                CHECK_LABEL(qv, dv);
                if (!m.data_mapped(dv))
                {
                  if constexpr (HAV)
                  {
                    MAP_UNCHECKED(m, qv, dv);
                    if (!check_anti_vertices(m)) continue;
                  }
                  c += 1;
                }
              }
              cursor = start;
              m.unmap(qv);
              idx -= 1;
            }
            else
            {
              if (cursor < end)
              {
                uint32_t dv = candidates[cursor++];
                CHECK_LABEL(qv, dv);
                if (!m.data_mapped(dv))
                {
                  MAP_UNCHECKED(m, qv, dv);
                  // in the same order group, can't have different labels or anything
                  // so we won't miss anything by setting cursors[idx+1] to cursor
                  if (oidxs[idx+1] == oidx) cursors[idx+1] = cursor;
                  idx += 1;
                }
              }
              else if (idx == 0)
              {
                return c;
              }
              else
              {
                cursor = start;
                m.unmap(qv);
                idx -= 1;
              }
            }
          }
        }
      }
      return 0;
    }

    template <Graph::Labelling L>
    void get_next_cand(std::vector<uint32_t> &xsection, unsigned idx, const partial_match<L> &mvgs)
    {
      // query vertices already matched, and connected from qo[key]
      std::vector<uint32_t> Ucon;
      const auto &true_edges = vgs.get_neighbours(qo[idx]);
      // all vertices in the range [levels, key) have been matched
      // want all such vertices that are adjacent to qo[key]
      for (const uint32_t m : true_edges) {
        for (uint32_t qi = 0; qi < idx; ++qi) {
          if (m == qo[qi]) {
            assert(mvgs.mapped(m));
            Ucon.push_back(m);
          }
        }
      }

      const adjlist &adj = gpb->get_adj(mvgs.at(Ucon[0]));

      uint32_t *start = adj.ptr;
      uint32_t *end = adj.ptr + adj.length;

      uint32_t lower_bound = 0;
      uint32_t upper_bound = static_cast<uint32_t>(-1); // largest possible
      for (uint32_t i = 0; i < idx; ++i) {
        if (qo[i] < qo[idx]) lower_bound = std::max(lower_bound, qo[i]);
        else if (qo[i] > qo[idx]) upper_bound = std::min(upper_bound, qo[i]);
      }

      start = lower_bound > 0
        ? std::upper_bound(adj.ptr, adj.ptr + adj.length, mvgs.at(lower_bound))
        : start;

      end = upper_bound < (uint32_t)-1
        ? std::lower_bound(start, adj.ptr + adj.length, mvgs.at(upper_bound))
        : end;

      xsection.assign(start, end);

      for (unsigned i = 1; i < Ucon.size(); ++i) {
        std::vector<uint32_t> t;
        t.swap(xsection);

        const uint32_t n = Ucon[i];
        const uint32_t v = mvgs.at(n);
        const adjlist &adj = gpb->get_adj(v);

        std::set_intersection(std::execution::unseq,
            t.cbegin(), t.cend(),
            adj.ptr, adj.ptr + adj.length,
            std::back_inserter(xsection));
      }


      std::vector<uint32_t> nUcon;
      const auto &not_edges = vgs.get_anti_neighbours(qo[idx]);
      // all vertices in the range [levels, key) have been matched
      // want all such vertices that are adjacent to qo[key]
      for (const uint32_t m : not_edges) {
        for (uint32_t qi = 0; qi < idx; ++qi) {
          if (m == qo[qi]) {
            assert(mvgs.mapped(m));
            nUcon.push_back(m);
          }
        }
      }

      for (const uint32_t n : nUcon) {
        std::vector<uint32_t> t;
        t.swap(xsection);

        const uint32_t v = mvgs.at(n);
        const adjlist &adj = gpb->get_adj(v);

        std::set_difference(std::execution::unseq,
            t.cbegin(), t.cend(),
            adj.ptr, adj.ptr + adj.length,
            std::back_inserter(xsection));
      }
    }

    template <Graph::Labelling L>
    uint64_t star_map_into(uint32_t vertex_id)
    {
      partial_match<L> m(rbi.query_graph.num_vertices());
      if constexpr (L == Graph::UNLABELLED || L == Graph::LABELLED)
      {
        m.map(rbi.vmap[0][1][0], vertex_id);
      }
      else
      {
        m.map(rbi.vmap[0][1][0], vertex_id, gpb->label(vertex_id));
      }


      if constexpr (L == Graph::DISCOVER_LABELS)
      {
        assert(false);
        return 0;
      }
      else
      {
        return complete_match<L>(m);
      }
    }

    template <Graph::Labelling L>
    uint64_t map_into(uint32_t vertex_id)
    {
      const uint8_t vgs_size = vgs.num_vertices();
      if constexpr (L == Graph::LABELLED)
      {
        if (gpb->label(vertex_id) != vgs.label(qo[0])) return 0;
      }
      else if constexpr (L == Graph::PARTIALLY_LABELLED)
      {
        if ((vgs.label(qo[0]) == static_cast<uint32_t>(-1) && gpb->known_label(gpb->label(vertex_id)))
            || (vgs.label(qo[0]) != static_cast<uint32_t>(-1) && gpb->label(vertex_id) != vgs.label(qo[0])))
        {
          return 0;
        }
      }

      if (vgs_size == 1) return star_map_into<L>(vertex_id);

      partial_match<L> mvgs(vgs_size);
      if constexpr (L == Graph::UNLABELLED || L == Graph::LABELLED)
      {
        mvgs.map(qo[0], vertex_id);
      }
      else
      {
        mvgs.map(qo[0], vertex_id, gpb->label(vertex_id));
      }

      uint8_t idx = 1;
      uint32_t *cursors = new uint32_t[vgs_size-1](); // last cursor unused?

      get_next_cand<L>(cands[idx], idx, mvgs);

      uint64_t c = 0;
      if (idx == vgs_size-1)
      {
        const auto &candidates = cands[idx];
        const uint32_t key = qo[idx];
        for (const uint32_t dv : candidates)
        {
          MAPR(mvgs, key, dv);

          c += nonCoreVertexMatching(mvgs);
        }
      }
      else
      {
        while (idx > 0) {
          const auto &candidates = cands[idx];
          const uint32_t key = qo[idx];
          if (idx == vgs_size-1)
          {
            for (const uint32_t dv : candidates)
            {
              MAPR(mvgs, key, dv);
              c += nonCoreVertexMatching(mvgs);
            }
            mvgs.unmap(key);
            idx -= 1;
          }
          else
          {
            uint32_t &cursor = cursors[idx];
            if (cursor < candidates.size())
            {
              uint32_t dv = candidates[cursor++];
              if (!mvgs.data_mapped(dv))
              {
                // prepare next level candidates
                MAPR(mvgs, key, dv);
                idx += 1;
                get_next_cand<L>(cands[idx], idx, mvgs);
              }
            }
            else
            {
              cursor = 0;
              mvgs.unmap(key);
              idx -= 1;
            }
          }
        }
      }

      delete[] cursors;

      return c;
    }


    const AnalyzedPattern &rbi;
    unsigned vgsi;
    const SmallGraph &vgs;
    const DataGraph * const gpb;
    std::vector<uint32_t> qo;
    std::vector<std::vector<uint32_t>> &cands;
  };

  uint32_t num_mappings(SmallGraph &data_graph, const SmallGraph &p)
  {
    std::vector<std::vector<uint32_t>> cands(p.num_vertices()+1, std::vector<uint32_t>{});

    DataGraph dg(data_graph);
    dg.set_rbi(p);

    uint32_t num_vertices = dg.get_vertex_count();
    uint32_t vgs_count = dg.get_vgs_count();

    uint64_t count = 0;
    const auto process = [&count](const CompleteMatch &) -> void { count += 1; };
    for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
    {
      Matcher<false, UNSTOPPABLE, decltype(process)> m(dg.rbi, &dg, vgsi, cands, process);
      for (uint32_t v = 1; v <= num_vertices; ++v)
      {
        m.map_into<Graph::UNLABELLED>(v);
      }
    }
    return count;
  }
}

#endif
