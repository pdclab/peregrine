#include <string>
#include "Graph.hh"

auto from_string(const std::string &str)
{
  Peregrine::SmallGraph p;
  {
    unsigned c = 0;
    while (c < str.size())
    {
      if (str[c] == '[')
      {
        unsigned edge_start = c;
        while (str[c] != ']' && c < str.size()) c += 1;
        // edges are in format [u-v] where each is a one digit integer
        assert(c-edge_start == 4);

        uint32_t u = str[edge_start+1] - '0';
        uint32_t v = str[edge_start+3] - '0';

        p.add_edge(u, v);
      }
      else if (str[c] == '(')
      {
        unsigned edge_start = c;
        while (str[c] != ')' && c < str.size()) c += 1;
        // edges are in format (u~v) where each is a one digit integer
        assert(c-edge_start == 4);

        uint32_t u = str[edge_start+1] - '0';
        uint32_t v = str[edge_start+3] - '0';

        p.add_anti_edge(u, v);
      }

      c += 1;
    }
  }

  return p;
}

auto from_string_labelled(const std::string &str)
{
  Peregrine::SmallGraph p;
  {
    unsigned c = 0;
    while (c < str.size())
    {
      if (str[c] == '[')
      {
        unsigned edge_start = c;
        while (str[c] != ']' && c < str.size()) c += 1;
        // edges are in format [u,lu-v,lv] where each is a one digit integer
        assert(c-edge_start == 8);

        uint32_t u = str[edge_start+1] - '0';
        uint32_t lu = str[edge_start+3] - '0';
        uint32_t v = str[edge_start+5] - '0';
        uint32_t lv = str[edge_start+7] - '0';

        p.add_edge(u, v);
        p.set_label(u, lu);
        p.set_label(v, lv);
      }
      else if (str[c] == '(')
      {
        unsigned edge_start = c;
        while (str[c] != ')' && c < str.size()) c += 1;
        // edges are in format (u,lu~v,lv) where each is a one digit integer
        assert(c-edge_start == 8);

        uint32_t u = str[edge_start+1] - '0';
        uint32_t lu = str[edge_start+3] - '0';
        uint32_t v = str[edge_start+5] - '0';
        uint32_t lv = str[edge_start+7] - '0';

        p.add_anti_edge(u, v);
        p.set_label(u, lu);
        p.set_label(v, lv);
      }
      c += 1;
    }
  }

  return p;
}

