#include <math.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iomanip>
#include "Network.h"
////////////////////////////////////////////////////////////////////////////////////////////////
Edge::Edge(){
	origin = NULL; Destn = NULL; EdgeWeight =0;}
	
Edge::Edge(Node* Org, Node* Dest, double weight){
	origin = Org; Destn = Dest; EdgeWeight = weight;}

double Edge::GetEdgeWeight() const{return EdgeWeight;}

Node* Edge::GetOrgin() {return origin;}

Node* Edge::GetDestn() {return Destn;}

Edge::~Edge(){}
////////////////////////////////////////////////////////////////////////////////////////////////
Node::Node(){NodeName =' '; NodeWeight =0; Edges = vector<Edge>();}

void Node::addEdge(Node* Neighbour, double weight){
	Edge newEdge(this, Neighbour, weight);
	Edges.push_back(newEdge);}

void Node::printEdges() const {
	cout<<NodeName<<" : ";
	 for (unsigned int i = 0; i < Edges.size(); i++){
            Edge e = Edges.at(i);
            cout << e.GetDestn()->getName()<<" - " << e.GetEdgeWeight() << endl;}
    cout << endl;}

void Node::setName(string name) {NodeName =name;}

void Node::setWeight(double weight) {NodeWeight = weight;}

vector<Edge> Node::GetEdges() {return Edges;}

double Node::GetNodeWeight() const {return NodeWeight;}

string Node::getName() const {return NodeName;}

Node::~Node(){}
//////////////////////////////////////////////////////////////////////////////////////////////////
Graph::Graph(){graph = vector<Node*>();}

void Graph::insertNode(Node* newNode){
	graph.push_back(newNode);}

void Graph::printGraph() const{
	for (unsigned int i=0; i<graph.size(); i++)
	graph.at(i)->printEdges();}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getCorrCoeff (int sampleSize, int orgNodeIndex, int destNodeIndex, double** geneParam);

int main()
{
	string geneFile = "METABRIC_Exp_common.txt", PPIfile = "iRefIndex_PPI_Edited.txt", CNVfile= "METABRIC_CNV_common_transpose.txt"; //Used R to get data pre-processing
	int no_Of_Gene = 13038, sampleSize = 997;
	double mean =0.0, sd =0.0, w, cutoff = 0.0;
	int i, k=0, maxItrNum;
	int no_Of_Proteins;
	
	cout<<"Provide maximum number of iteration ";
	cin>>maxItrNum;
	cout<<"Provide correlation cutoff ";
	cin>>cutoff;
	
	// Create variables to store geneName, geneExpression, mean and standard deviation over Sample size	
	string* geneNameList = new string[no_Of_Gene];
	double **geneParam = new double*[no_Of_Gene];
	for (i = 0; i < no_Of_Gene; i++)
		geneParam[i] = new double[sampleSize + 2];
	
	// Read gene expression file
	ifstream geneRead;
	geneRead.open(geneFile.c_str());
	cout<<"Reading header of "<<geneFile<<endl;
	string s;
	for (i = 0; i<sampleSize+1; i++) // Passing over the header (Sample names)
		geneRead>>s;
	//Calculate mean and standard deviation for gene expression
	while(!geneRead.eof()){
		geneRead>>geneNameList[k];
		mean = 0.0;
		sd = 0.0;
		for(i =0; i<sampleSize; i++){
			geneRead>>geneParam[k][i];
			mean +=geneParam[k][i];
		}
		mean/= sampleSize;
		for(i = 0; i<sampleSize; i++)
			sd += pow(geneParam[k][i]-mean,2);
		sd = sqrt(sd/sampleSize);
		geneParam[k][sampleSize] = mean;
		geneParam[k][sampleSize+1] = sd;
		k++;
	}
	
	cout<<"Mean and standard deviation calculated."<<endl;
	geneRead.close();
	cout<<"Closing file "<<geneFile<<endl;
	
	no_Of_Proteins = k;
	Node* Vertex = new Node[no_Of_Proteins];
	cout<< no_Of_Proteins << " nodes are created" <<endl;
	for (i = 0; i<no_Of_Proteins; i++){
		Vertex[i].setName(geneNameList[i]);
		Vertex[i].setWeight(0.0);	
		}
	
	bool **EdgeCheck = new bool*[no_Of_Proteins];
	for (i = 0; i<no_Of_Proteins; i++) EdgeCheck[i] = new bool[no_Of_Proteins];
	for (i = 0; i<no_Of_Proteins; i++){
		for (k=0; k<no_Of_Proteins; k++)
		EdgeCheck[i][k]= false;}
	
	ifstream proteinRead;
	proteinRead.open(PPIfile.c_str());
	cout<<"Reading Protein Protein interaction file"<<endl;
	cout<<"Connecting Vertices with edges "<<endl;
	string x,y;
	int orgNodeIndex, destNodeIndex;
	int totalEdges=0;
	bool flag1, flag2;
	double coVariance, corrCoeff;
	while(!proteinRead.eof()){
		proteinRead>>x;
		proteinRead>>y;
		if(x.compare(y)!=0){
			flag1 = false; flag2 = false;
			// Finding edge between two Nodes(Genes)
			for(i=0; i<no_Of_Proteins; i++){
                if(geneNameList[i]==x && flag1==false) {orgNodeIndex=i; flag1=true;}
                if(geneNameList[i]==y && flag2==false) {destNodeIndex=i; flag2=true;}
                if(flag1==true && flag2==true) break;
            }
            // Calculating weight of the edge using Pearson Correlation coefficient
			if(flag1==true && flag2==true && EdgeCheck[orgNodeIndex][destNodeIndex]==false){
				corrCoeff = getCorrCoeff(sampleSize, orgNodeIndex, destNodeIndex, geneParam);
			    // Assgining edge to the graph
			    if( Vertex[orgNodeIndex].getName() != Vertex[destNodeIndex].getName()){ //Just in case if gene names in expression table are same
                	Vertex[orgNodeIndex].addEdge( &(Vertex[destNodeIndex]), corrCoeff );
                	Vertex[destNodeIndex].addEdge( &(Vertex[orgNodeIndex]), corrCoeff );
                	totalEdges++;
				}
				EdgeCheck[orgNodeIndex][destNodeIndex]=true;
				EdgeCheck[destNodeIndex][orgNodeIndex]=true;
			}
		}
	}
	
	double degreeOfGraph = double(totalEdges)/double(no_Of_Proteins);
	cout<< "Degree of Network is " << degreeOfGraph << endl;
	cout<< "Total edges in the network are : " << totalEdges <<" edges " << endl;
	proteinRead.close();
	cout<< PPIfile <<" has been closed" <<endl; 
	
    // Creating graph
    Graph PPI_network;
    for(i = 0; i< no_Of_Proteins; i++)
    PPI_network.insertNode(&Vertex[i]);
    cout<< "PPI unweighted graph has been created " << endl;
   
   // The following code is to find out trans-associated genes
   // The corresponding equations can be found in referred paper in section 4.2: equation 2 and 3
   /////////////////////////////////////////////////////////////////////////////////////////////////// Assigning weights to Vertices
    ifstream CNV_read;
	CNV_read.open(CNVfile.c_str());
	cout<<"Reading " <<CNVfile << endl;
	// Just to make sure that Genes are common in CNV and expression file
	s.clear();
	CNV_read>>s;
	for (i=0; i<no_Of_Proteins;i++){
		CNV_read>>s;
		if(geneNameList[i] != s){
			cout<<"Gene names not matched" <<endl;
			return 1;}
		}
	cout<<"Gene names matched" <<endl;
	s.clear();
	
	string a,b, fileName;
	double lambda =0.0; //Penalizing factor
	int edgeSign, CurNodeSign, destnNodeSign; //sign of w(gx, gy), gx and gy
	double *NghbrContri = new double[no_Of_Proteins];
	vector<Edge> NeighborSet;
	double* PositiveContri = new double[no_Of_Proteins];
	double* NegativeContri = new double[no_Of_Proteins];
	for (i = 0; i< no_Of_Proteins; i++){
		PositiveContri[i] =0.0;
		NegativeContri[i] =0.0;}
	
	ofstream weight_out;
	k=1; //Indicating patient index
	while(!CNV_read.eof())
	{
		CNV_read>>s;
		cout<<"Read Patient #"<<k<<"\t"<<s<<endl;
		for (i=0; i<no_Of_Proteins; i++)
		{
			CNV_read>>w;
			Vertex[i].setWeight(w); // Each node is assigned some initial weight
		}
		a=static_cast<ostringstream*>(&(ostringstream()<<k))->str();
		for(int itr=0; itr<maxItrNum; itr++)
		{
			for(i=0; i<no_Of_Proteins; i++){
				lambda =0.0;
				NghbrContri[i] = 0.0;
				NeighborSet = Vertex[i].GetEdges();
				if(Vertex[i].GetNodeWeight() > 0.0) CurNodeSign =1;
				else if(Vertex[i].GetNodeWeight() < 0.0) CurNodeSign = -1;
				else  CurNodeSign = 0;
				for (int j =0; j<NeighborSet.size(); j++){
					
					if(NeighborSet[j].GetEdgeWeight() > 0.0) edgeSign = 1;
                    else if(NeighborSet[j].GetEdgeWeight() < 0.0) edgeSign = -1;
                    else edgeSign = 0;
                    
                    if(NeighborSet[j].GetDestn()->GetNodeWeight() > 0.0) destnNodeSign = 1;
                    else if(NeighborSet[j].GetDestn()->GetNodeWeight() < 0.0) destnNodeSign = -1;
                    else destnNodeSign = 0;
                    
					if( fabs(NeighborSet[j].GetEdgeWeight()) >= cutoff && CurNodeSign != 0 && edgeSign != 0 && destnNodeSign != 0 && (edgeSign * destnNodeSign == CurNodeSign))
						{   //Calculating numerator equation 2
						NghbrContri[i] += ( NeighborSet[j].GetDestn()->GetNodeWeight() * (double)edgeSign * pow(fabs(NeighborSet[j].GetEdgeWeight()),itr+1) ); 
						lambda = lambda + 1.0;
						}
					}
				if(lambda>0.0){
					lambda = MAX(0.0, (double)(degreeOfGraph - NeighborSet.size()) );
					NghbrContri[i] /= (lambda + (double)(NeighborSet.size()));
				}
				else NghbrContri[i] = 0.0;
			}
			// Update node weight for all genes after each iteration
			for(i=0; i<no_Of_Proteins; i++){
				w = Vertex[i].GetNodeWeight() + NghbrContri[i];
				Vertex[i].setWeight(w);
				if(NghbrContri[i] > 0.0) PositiveContri[i] += NghbrContri[i];
				else if (NghbrContri[i] < 0.0) NegativeContri[i] += NghbrContri[i];
			}
			// Get the data out for each patients
			b=static_cast<ostringstream*>(&(ostringstream()<<i+1))->str();
			fileName="Patient_"+a+"_Iteration_"+b+".txt";
			if(itr = maxItrNum-1)
			{
				weight_out.open(string("METABRIC_" + fileName).c_str());
                weight_out<<"Gene Weight\n";
                for(i=0; i<no_Of_Proteins; i++){
                    if(i<no_Of_Proteins-1) weight_out<<Vertex[i].getName()<<" "<<std::setprecision(4)<<Vertex[i].GetNodeWeight()<<endl;	//4-decimal precision.
                    else weight_out<<Vertex[i].getName()<<" "<<std::setprecision(4)<<Vertex[i].GetNodeWeight();
                }
                weight_out.close();
			}
		}
		s.clear();
		k++;
	}

    CNV_read.close();
    cout<< "Data is ready for k-means clustering" <<endl;
	return 0;
}

double getCorrCoeff (int sampleSize, int orgNodeIndex, int destNodeIndex, double** geneParam){
	double coVariance = 0.0, corrCoeff;
	for(int i = 0; i<sampleSize; i++){
		coVariance += (geneParam[orgNodeIndex][i]-geneParam[orgNodeIndex][sampleSize]) * (geneParam[destNodeIndex][i]-geneParam[destNodeIndex][sampleSize]);
		coVariance /= sampleSize;
	}
	corrCoeff = coVariance / (geneParam[orgNodeIndex][sampleSize+1]*geneParam[destNodeIndex][sampleSize+1]);
	return corrCoeff;
}



