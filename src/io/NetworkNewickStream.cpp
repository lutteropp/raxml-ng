#include "file_io.hpp"

using namespace std;

char * newick_name_cb(const pll_unetwork_node_t * node)
{
  return node->label ? strdup(node->label) : strdup("");
}

char * newick_print_cb(const pll_unetwork_node_t * node)
{
  // that's ugly, but cannot find a better solution so far...
  const unsigned int precision = logger().precision(LogElement::brlen);

  char * newick;

  if (node_is_reticulation(node))
  {
    if (asprintf(&newick, "%s:%.*lf",
             node->label ? node->label : "" , precision, node->length) < 0)
      return NULL;
  }
  else
  {
    // TODO: implement this...
	return NULL;
  }

  return newick;
}

std::string to_newick_string(const Network& network)
{
  char * newick_str = pll_unetwork_export_newick(&network.pll_unetwork_root(), NULL);
  std::string result = newick_str;
  free(newick_str);
  return result;
}

NetworkNewickStream& operator<<(NetworkNewickStream& stream, const pll_unetwork_node_t& root)
{
  auto print_cb = stream.brlens() ? newick_print_cb : newick_name_cb;
  char * newick_str = pll_unetwork_export_newick(&root, print_cb);
  if (newick_str)
  {
    stream << newick_str << std::endl;
    free(newick_str);
  }
  else
  {
    assert(pll_errno);
    libpll_check_error("Failed to generate Newick");
  }
  return stream;
}

NetworkNewickStream& operator<<(NetworkNewickStream& stream, const pll_unetwork_t& network)
{
  stream << *network.vroot;
  return stream;
}

NetworkNewickStream& operator<<(NetworkNewickStream& stream, const Network& network)
{
  stream << network.pll_unetwork_root();
  return stream;
}

NetworkNewickStream& operator>>(NetworkNewickStream& stream, Network& network)
{
  string newick_str;

  std::getline(stream, newick_str, ';');

  // discard any trailing spaces, newlines etc.
  stream >> std::ws;

  if (!newick_str.empty())
  {
    newick_str += ";";

    pll_unetwork_t * unetwork = pll_unetwork_parse_newick_string(newick_str.c_str());

    libpll_check_error("ERROR reading network file");

    assert(unetwork);

    network = Network(*unetwork);

    pll_unetwork_destroy(unetwork, nullptr);
  }

  return stream;
}

/*NetworkNewickStream& operator<<(NetworkNewickStream& stream, const SupportNetwork& network)
{
  stream << network.pll_unetwork_root();
  return stream;
}*/


