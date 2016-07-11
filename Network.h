// Header file
#ifndef NETWORK_H 
#define NETWORK_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Node;
class Edge;
class Graph;

class Edge{
	private:
		Node* origin;
		Node* Destn;
		double EdgeWeight;

	public:
		Edge(); //Constructor
		Edge(Node* Org, Node* Dest, double weight); // Overload constructor
		double GetEdgeWeight() const; // Getting Edge weight
		Node* GetOrgin();
		Node* GetDestn();
		~Edge(); //Destructor
};

class Node{
	private:
		string NodeName;
		double NodeWeight;
		vector<Edge> Edges;
		
	public:
		Node();
		void addEdge(Node* Neighbour, double weight);
		void setName(string);
		void setWeight(double);
		void printEdges() const;
		vector<Edge> GetEdges();
		double GetNodeWeight() const;
		string getName() const;
		~Node();
};

class Graph{
	private:
		vector<Node*> graph;
		
	public:
		Graph();
		void insertNode(Node* newNode);
		void printGraph() const;		
};

#define MAX(a,b) (a>b?a:b)

#endif
