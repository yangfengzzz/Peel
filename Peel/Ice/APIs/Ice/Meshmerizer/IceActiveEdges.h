///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code to create active edges.
 *	\file		IceActiveEdges.h
 *	\author		Pierre Terdiman
 *	\date		January, 17, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEACTIVEEDGES_H
#define ICEACTIVEEDGES_H

#ifdef SUPPORT_SEPARATE_ACTIVE_EDGES
class MESHMERIZER_API ActiveEdges {
public:
    ActiveEdges();
    ~ActiveEdges();

    bool Compute(const EdgeList& edges, const IndexedTriangle* faces, const Point* verts, float epsilon = 0.001f);

    // Data access
    inline_ const bool* GetActiveEdges() const { return mActiveEdges; }
    inline_ bool GetActiveEdge(udword edge_index) const { return mActiveEdges[edge_index]; }

private:
    bool* mActiveEdges;  //!< mNbEdges bools marking active edges
};
#endif

#endif  // ICEACTIVEEDGES_H
