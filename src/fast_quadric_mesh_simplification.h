//////////////////////////////////////////////////////////////////////////
//
// OpenMesh based fast quadric mesh simplification
//
// (C) by grassofsky in 2022
//
// License : MIT
// http://opensource.org/licenses/MIT
// 
// https://github.com/grassofsky/OpenMesh-Fast-Quadric-Simplification
//
//////////////////////////////////////////////////////////////////////////

#pragma once 

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES_WAS_DEFINED
#define _USE_MATH_DEFINES
#endif

#include <iostream>

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Utils/Property.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

#include "symetric_matrix.h"

#ifdef _USE_MATH_DEFINES_WAS_DEFINED
#undef _USE_MATH_DEFINES
#endif

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class MeshSimplification
{
public:
    MeshSimplification(MyMesh& mesh);
    ~MeshSimplification();

    void SimplifyVertexTo(size_t uiRemainedVertexNum, double dAgressiveness = 7);
    void SimplifyFaceTo(size_t uiRemainedFaceNum, double dAgreesiveness = 7);

private:
    void Initialize();
    double CalculateError(OpenMesh::EdgeHandle edge, MyMesh::Point& ptResult);
    // v*Q*v
    double VertexError(const SymetricMatrix& q, double x, double y, double z);
    /// \brief When v0 collapsed to ptTarget, test if there will be some flipped triangles around v0
    bool IsFlipped(OpenMesh::VertexHandle v0, OpenMesh::VertexHandle v1, const MyMesh::Point& ptTarget);
    void UpdateFaceNormal(OpenMesh::VertexHandle v0);
    void UpdateEdgePropertyAroundV(OpenMesh::VertexHandle v0);

private:
    MyMesh& m_mesh;

    OpenMesh::VPropHandleT<SymetricMatrix>  m_vpQuadrics;
    OpenMesh::EPropHandleT<double> m_epError;
    OpenMesh::EPropHandleT<MyMesh::Point> m_epTargetPoints;
    OpenMesh::EPropHandleT<bool> m_epDirty;
};
