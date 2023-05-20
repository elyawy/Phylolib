#ifndef ___SEQUENCE_CONTAINER_TREE_MAP
#define ___SEQUENCE_CONTAINER_TREE_MAP
#include "definitions.h"
#include "tree.h"
#include "treeIt.h"
#include "sequenceContainer.h"

// Explanations by Tal Pupko 2/6/2017.
// This class connects the tree and the sequence container. Each node in the tree has an ID and each sequence in the sequence container has an ID.
// The class provides means to find the sequence ID that is associated with a specific node ID.
// It stores this information in a vector _V so that _V[nodeID] = the sequence ID that is associated with the sequence in that nodeID.
// _V of internal nodes is set to -1.

// CLASS UTILITY FUNCTIONS (static functions, no need to instantiate a class to use them).
// The class also provides some utility functions, for example to check that there's an agreement between the sequence container names 
// and the names of the node sequences.
// In the function checkThatNamesInTreeAreSameAsNamesInSequenceContainer, bLeavesOnly = true indicates that agreement is only searched
// between the names of the leaves of the tree and the sequence container.
// If it is false, that the function issue an error also if the sequence container does not include sequences of internal nodes;
// Usually, only leaf sequences are strored in the sequence container.
// But sometimes, also sequences of internal nodes are stored (e.g., when reconstructing ancestral sequences, or when simulating sequences).
// This utility function thus also handles the case in which internal nodes in the tree are asociated with sequences in the sequence container.

// A DESCRIPTION OF THE FUNCTION < intersectNamesInTreeAndSequenceContainer > IS MISSING. TO DO
// A DESCRIPTION OF THE FUNCTION < getLeavesSequences > IS MISSING. TO DO

class seqContainerTreeMap {
public:
	explicit seqContainerTreeMap(const sequenceContainer& sc, const tree& et);
	seqContainerTreeMap() {};
	void fillSeqContainerTreeMap(const sequenceContainer& sc, const tree& et);

	int seqIdOfNodeI(const int nodeID) const {return _V[nodeID];}

	static void checkThatNamesInTreeAreSameAsNamesInSequenceContainer(const tree& et, const sequenceContainer & sc, bool bLeavesOnly = true);
	static void intersectNamesInTreeAndSequenceContainer(tree& et, sequenceContainer & sc, bool bLeavesOnly);
	static void getLeavesSequences(const sequenceContainer& sc, const tree& tr, sequenceContainer& sc_leaves);
private:
	vector<int> _V;// _V[i] is the sequenceId of node I.
};

#endif
