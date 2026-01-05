/*
 * Author: qopgx (real name to be disclosed)
 * */
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using natural = std::uint32_t;
using ixt = natural;
using label = natural;
using degree = natural;

using vertex_neighbours = std::vector<label>;
using degree_bucket = std::vector<label>;

/*
 * Bucket sort for Chiba-Nishizeki groups by degrees instead of ranges.
 * So this is the special case of bucket sort for ranges with only one integer
 * value. All buckets are thus automatically sorted.
 *
 * How we specifically implement bucket sort in this algorithm is by first
 * statically sorting the elements into buckets, but then dynamically fetching
 * the ordering by tracking the highest nonempty degree, popping constantly
 * from the back O(1), which also means that we don't need to keep around a
 * running index for the current nonempty degree bucket.
 * */
class Buckets {
  natural n;
  degree_bucket *buckets;

  void _initialize_buckets(degree *degrees) {
    for (ixt i = 0; i < n; i++) {
      buckets[degrees[i]].push_back(i);
    }
  };

public:
  degree highest_degree_nonempty;

  Buckets(natural n, degree *degrees)
      : n(n), buckets(new degree_bucket[n]()), highest_degree_nonempty(n - 1) {
    _initialize_buckets(degrees);
  }
  ~Buckets() { delete[] buckets; }

  bool pop_next_vx_in_order(label &result) {
    degree_bucket *bucket;
    label l;
    while (highest_degree_nonempty >= 0) {
      if (!(bucket = &buckets[highest_degree_nonempty])->empty()) {
        result = bucket->back();
        bucket->pop_back();
        return true;
      } else
        highest_degree_nonempty--;
    }
    return false;
  }
};
/*
 * The input of the program is structured as follows:
 *
 * |V(G)| IS1INDEXED?
 * i j
 *  .
 *  .
 *  .
 *
 * Where V(G) is the number of vertices in the graph, ISI1INDEXED?
 * is either 0 or 1, representing whether the graph is 1- or 0-indexed,
 * and the following rows contain the labels of the vertices of the edges
 * of the graph. Input terminates on an empty line.
 *
 *                          PLEASE NOTE:
 *
 * The program does not contain out of bounds checks for invalid INPUT values,
 * errors in case wrong indexation has been specified, other error checking etc.
 *
 * The program is made to work as-is for valid input values.
 *
 * An example in zachary_karate_graph.txt is given.
 * */
int main() {
  // The number of vertices
  natural n;
  natural indexation;
  std::cin >> n >> indexation;
  std::cin.ignore();

  // Because of O(d(v_i)) removal
  vertex_neighbours *adjacency = new vertex_neighbours[n]();

  /*
   * Chiba-Nishizeki guarantees O(d(v_i)) removal of elements,
   * this can only work if we store the indices of the removed vertex
   * in all of the other adjacency lists, meaning that when iterating
   * a removed vertex adjacency list, we use the current index 0 <= i < deg(v)
   * to find mirrors[rem_label][i], which must be the index of rem_label in the
   * adjacency list of the vertex with label adjacency[rem_label][i].
   *
   * adjacency[removed_label][ix_neighbour_in_rem] = label_neighbour
   * mirrors[removed_label][ix_neighbour_in_rem] = ix_rem_in_neighbour
   *
   * Not storing mirrors would require re-discovery of ix_rem_in_neighbour
   * through iterative search through neighbour adjacency list,
   * worst-case O(deg(neighbour)).
   *
   * This would make the program slower for high-degree nodes.
   * */
  std::vector<ixt> *mirrors = new std::vector<ixt>[n]();

  // Degrees
  degree *degrees = new degree[n]();

  // Start loading graph
  std::string edge;

  /*
   * Just as specified in the original paper the degrees of vertices are
   * computed in O(m) time, where m is the number of edges.
   *
   * We do this by incrementing the vertex degrees for each incoming edge.
   *
   *    'Clearly the degrees of vertices can be computed in O(m) time.'
   * */
  while (std::getline(std::cin, edge)) {
    std::stringstream ss(edge);
    label l1, l2;
    ss >> l1 >> l2;

    if (indexation) {
      l1--;
      l2--;
    }

    ixt ix_l1 = adjacency[l1].size();
    ixt ix_l2 = adjacency[l2].size();

    adjacency[l1].push_back(l2);
    adjacency[l2].push_back(l1);

    mirrors[l1].push_back(ix_l2);
    mirrors[l2].push_back(ix_l1);

    // Handshaking lemma
    ++degrees[l1];
    ++degrees[l2];
  }

  /* Chiba-Nishizeki specifies bucket sort as an average-case O(n)
   * sorting algorithm for creating the required vertex order by degree.
   *
   * Specifically:
   *
   *   'sort the vertices v_1, v_2, ..., v_n of G in such a way that d(v_1) >=
   *    d(v_2) >= ... >= d(v_n)'
   * */

  Buckets buckets(n, degrees);

  bool *mark = new bool[n]();

  natural triangle_count = 0;

  /*
   * Loop indices are equivalent to the bounds 1 <= i <= n-2 in the original
   * algorithm.
   * */
  for (ixt i = 0; i < n - 2; i++) {
    label vx;

    /*
     * ORDER BASED VERTEX SELECTION
     *
     * Select max d(v1) in G from the leftover vertices in G
     *
     * This step also automatically handles vertex removal
     * */

    buckets.pop_next_vx_in_order(vx);

    /*
     * MARKING
     *
     * We mark all neighbours of vx to create the set A,
     * we compare marks for all neighbours of all vxs in A.
     * (These are the B sets).
     *
     * The intersection of A and B gives a triangle.
     * */
    vertex_neighbours &vx_adjacent_edges = adjacency[vx];

    for (ixt ixa = 0; ixa < vx_adjacent_edges.size(); ixa++) {
      mark[vx_adjacent_edges[ixa]] = true;
    }

    /*
     * TRIANGLE COUNTING
     *
     * If any neighbour of a neighbour vx is mark[nnvx] == true
     * then we have found a triangle.
     * */
    for (ixt ixa = 0; ixa < vx_adjacent_edges.size(); ixa++) {
      label nvx = vx_adjacent_edges[ixa];
      vertex_neighbours &nvx_adjacent_edges = adjacency[nvx];

      for (ixt ixan = 0; ixan < nvx_adjacent_edges.size(); ixan++) {
        label nnvx = nvx_adjacent_edges[ixan];

        if (mark[nnvx]) {
          std::cout << "(" << vx << "," << nvx << "," << nnvx << "), ";
          triangle_count++;
        }
      }

      mark[nvx] = false;
    }

    /*
     * VERTEX DELETION
     *
     * Finally remove vx from the graph.
     *
     * Removal complexity of vx adjacency records is O(d(vx))
     * due to the storage of mirror indices as described earlier.
     *
     * O(1) removal from some nvx_adjacency is guaranteed by swapping
     * (O(1)) the element into the back and pop_back (O(1)).
     *
     * Furthermore, besides respective adjacency records, the entire adjacency
     * list of vx should be cleared, this is possible O(1) due to the fact
     * that in CPP vector.clear() is O(1) for trivially destructed types like
     * scalars.
     *
     * Such vx will not appear in our adjacency lists further.
     * */
    for (ixt ixa = 0; ixa < vx_adjacent_edges.size(); ixa++) {
      label &nvx = vx_adjacent_edges[ixa];
      vertex_neighbours *nvx_adjacency = &adjacency[nvx];

      ixt mirror_ix = mirrors[vx][ixa];
      label nvx_last_neighbour = nvx_adjacency->back();

      std::swap(nvx_adjacency->at(mirror_ix), nvx_adjacency->back());

      /**
       * Once we remove vx with swap and pop we additionally have to correct
       * where mirror_ix of the of nvx points to in mirrors[nvx_last_neighbour],
       * due to the fact that index changed.
       *
       * The mirror of the last neighbour of nvx points to where nvx is in
       * the adjacency list of the last neighbour of nvx. At precisely the
       * same index in mirrors, we find where the last element of nvx is in
       * the adjacency list of nvx. Before correction, that integer value is
       * showing nvx_adjacency->size() - 1. This is wrong and precisely what
       * we are correcting.
       * */
      mirrors[nvx_last_neighbour][mirrors[nvx][nvx_adjacency->size() - 1]] =
          mirror_ix;

      mirrors[nvx][mirror_ix] = mirrors[nvx][nvx_adjacency->size() - 1];

      nvx_adjacency->pop_back();
    }
    vx_adjacent_edges.clear();
  }

  std::cout << std::endl << "Num K3 / C3: " << triangle_count << std::endl;

  // Memory cleanup
  delete[] adjacency;
  delete[] mirrors;
  delete[] degrees;
  delete[] mark;
}
