# Chiba-Nishizeki Triangle Search

This is an implementation of the original Chiba-Nishizeki triangle search algorithm introduced in the 1985 paper [ARBORICITY AND SUBGRAPH LISTING ALGORITHMS](https://www.cs.cornell.edu/courses/cs6241/2019sp/readings/Chiba-1985-arboricity.pdf).

It attempts to follow the pseudocode logic of the original algorithm, but it necessarily includes additional logic and data structures which enable the stated complexity. The ordering process is abstracted into a dedicated class.

## Motivation

The program was written as a personal exercise and for the educational value of having a sample implementation of the original algorithm.

## Running

Once the binary is compiled using CMake, you should first create an input `.txt` file in the format (a sample):

```
10 1
1 2
8 3
# AND SO ON...
```

The number 10 is the number of vertices of the graph, the second number tells the program whether the vertices are 0 or 1 indexed (labeled). Finally, all of the edge labels follow.

Please note the program will not error for invalid inputs (your machine will).

Run with:

```
./triangles < input.txt
```

## Testing

A sample test file `zachary_karate_graph.txt` is provided, from the popular [Zachary karate club](http://konect.cc/networks/ucidata-zachary/) test set, a description:

> This is the well-known and much-used Zachary karate club network. The data was collected from the members of a university karate club by Wayne Zachary in 1977. Each node represents a member of the club, and each edge represents a tie between two members of the club. The network is undirected. An often discussed problem using this dataset is to find the two groups of people into which the karate club split after an argument between two teachers.

Importantly, the **number of vertices is 34** and the **number of triangles is 45**.
