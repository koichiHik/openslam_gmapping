#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
//#include <gsl/gsl_blas.h>

#include <gmapping/utils/stat.h>
#include <gmapping/gridfastslam/gridslamprocessor.h>
#include <gmapping/gridfastslam/particle.h>

namespace GMapping {

using namespace std;

TNodeVector GridSlamProcessor::getTrajectories() const{
  TNodeVector v;
  TNodeMultimap parentCache;
  TNodeDeque border;
	
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
		TNode* node=it->node;
		while(node){
			node->flag=false;
			node=node->parent;
		}
	}
	
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
		TNode* newnode=new TNode(* (it->node) );
		
		v.push_back(newnode);
		assert(newnode->childs==0);
		if (newnode->parent){
			parentCache.insert(make_pair(newnode->parent, newnode));
			if (! newnode->parent->flag){
				newnode->parent->flag=true;
				border.push_back(newnode->parent);
			}
		}
	}
	
	while (! border.empty()){
		const TNode* node=border.front();
		border.pop_front();
		if (! node)
			continue;
			
		TNode* newnode=new TNode(*node);
		node->flag=false;
		
		//update the parent of all of the referring childs 
		pair<TNodeMultimap::iterator, TNodeMultimap::iterator> p=parentCache.equal_range(node);
		double childs=0;
		for (TNodeMultimap::iterator it=p.first; it!=p.second; it++){
			assert(it->second->parent==it->first);
			(it->second)->parent=newnode;
			childs++;
		}
		parentCache.erase(p.first, p.second);
		assert(childs==newnode->childs);
		
		//unmark the node
		if ( node->parent ){
			parentCache.insert(make_pair(node->parent, newnode));
			if(! node->parent->flag){
				border.push_back(node->parent);
				node->parent->flag=true;
			}	
		}
		//insert the parent in the cache
	}
	for (unsigned int i=0; i<v.size(); i++){
		TNode* node= v[i];
		while (node){
			node=node->parent;
		}
	}	
	
	return v;

}

void GridSlamProcessor::updateTreeWeights(bool weightsAlreadyNormalized){

  if (!weightsAlreadyNormalized) {
    normalize();
  }
  resetTree();
  propagateWeights();
}

void GridSlamProcessor::resetTree(){

	// Iterate for all particles in the vector.
	for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++){
		TNode* n = it->node;

		// Iterate all nodes in trajectory owned by THIS particle.
		while (n){
			n->accWeight = 0;
			n->visitCounter = 0;
			n = n->parent;
		}
	}
}

double propagateWeight(TNode* n, double weight){
	if (!n)
		return weight;
	double w=0;
	n->visitCounter++;
	n->accWeight+=weight;
	if (n->visitCounter==n->childs){
		w=propagateWeight(n->parent,n->accWeight);
	}
	assert(n->visitCounter<=n->childs);
	return w;
}

double GridSlamProcessor::propagateWeights(){

  // Don't calls this function directly, use updateTreeWeights(..) !
  // All nodes must be resetted to zero and weights normalized.
  
	// The accumulated weight of the root.
	double lastNodeWeight = 0;

	// Sum of the weights in the leafs.
	double aw = 0;                   

	std::vector<double>::iterator w = m_weights.begin();

	for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++){
		double weight = *w;
		aw += weight;
		TNode *n = it->node;
		n->accWeight = weight;
		lastNodeWeight += propagateWeight(n->parent,n->accWeight);
		w++;
	}
	
	if (fabs(aw - 1.0) > 0.0001 || fabs(lastNodeWeight - 1.0) > 0.0001) {
	  cerr << "ERROR: ";
	  cerr << "root->accWeight=" << lastNodeWeight << "    sum_leaf_weights=" << aw << endl;
	  assert(0);         
	}
	return lastNodeWeight;
}

};