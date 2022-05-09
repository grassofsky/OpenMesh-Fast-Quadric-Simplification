#include "fast_quadric_mesh_simplification.h"

Simplification::Simplification(MyMesh& mesh) : m_mesh(mesh)
{
    m_mesh.add_property(m_vpQuadrics);
    m_mesh.add_property(m_epError);
    m_mesh.add_property(m_fpError);
    m_mesh.add_property(m_epTargetPoints);
    m_mesh.add_property(m_fpDirty);
}

Simplification::~Simplification()
{
    m_mesh.remove_property(m_vpQuadrics);
    m_mesh.remove_property(m_epError);
    m_mesh.remove_property(m_fpError);
    m_mesh.remove_property(m_epTargetPoints);
    m_mesh.remove_property(m_fpDirty);
}

void Simplification::SimplifyVertexTo(int iRemainedVertexNum, double dAgressiveness/*=7*/)
{
    int nCurCollapses = 0;
    int nCollapses = m_mesh.n_vertices() - iRemainedVertexNum;

    if (!m_mesh.has_vertex_status()) m_mesh.request_vertex_status();
    if (!m_mesh.has_face_status())   m_mesh.request_face_status();
    if (!m_mesh.has_edge_status())   m_mesh.request_edge_status();
    if (!m_mesh.has_face_normals())
    {
        m_mesh.request_face_normals();
    }
    m_mesh.update_face_normals();

    Initialize();

    for (int iteration = 0; iteration < 100; iteration++)
    {
        if (nCurCollapses >= nCollapses) break;
        
        // Magic equation
        double dThreshold = 0.000000001 * pow(double(iteration + 3), dAgressiveness);

        for (auto fIt = m_mesh.faces_begin(); fIt != m_mesh.faces_end(); ++fIt)
        {
            m_mesh.property(m_fpDirty, *fIt) = false;
        }

        for (auto fIt = m_mesh.faces_begin(); fIt != m_mesh.faces_end(); ++fIt)
        {
            if (m_mesh.property(m_fpError, *fIt) > dThreshold) continue;
            if (!m_mesh.is_valid_handle(*fIt)) continue;
            if (m_mesh.status(*fIt).deleted()) continue;
            if (m_mesh.property(m_fpDirty, *fIt)) continue;

            for (auto hIt = m_mesh.fh_iter(*fIt); hIt.is_valid(); ++hIt)
            {
                auto eh = m_mesh.edge_handle(hIt);
                if (m_mesh.property(m_epError, eh) < dThreshold)
                {
                    OpenMesh::HalfedgeHandle h0 = *hIt;
                    OpenMesh::VertexHandle v0 = m_mesh.from_vertex_handle(h0);
                    OpenMesh::VertexHandle v1 = m_mesh.to_vertex_handle(h0);

                    if (m_mesh.is_boundary(v0) != m_mesh.is_boundary(v1))
                    {
                        continue;
                    }

                    MyMesh::Point ptTarget = m_mesh.property(m_epTargetPoints, eh);
                    if (IsFlipped(v0, v1, ptTarget)) continue;
                    if (IsFlipped(v1, v0, ptTarget)) continue;

                    auto h1 = m_mesh.opposite_halfedge_handle(h0);
                    if (m_mesh.is_collapse_ok(h0))
                    {
                        m_mesh.collapse(h0);
                        m_mesh.point(v1) = ptTarget;
                    }
                    else if (m_mesh.is_collapse_ok(h1))
                    {
                        m_mesh.collapse(h1);
                        m_mesh.point(v0) = ptTarget;
                    }
                    else
                    {
                        continue;
                    }
                    nCurCollapses++;

                    // Update related face normal
                    UpdateFaceNormal(v1);
                    UpdateFaceNormal(v0);
                    auto v0Quadric = m_mesh.property(m_vpQuadrics, v0);
                    auto v1Quadric = m_mesh.property(m_vpQuadrics, v1);
                    m_mesh.property(m_vpQuadrics, v1) += v0Quadric;
                    m_mesh.property(m_vpQuadrics, v0) += v1Quadric;
                    UpdateEdgePropertyAroundV(v1);
                    UpdateEdgePropertyAroundV(v0);
                    break;
                }
                if (nCurCollapses >= nCollapses) break;
            }
        }
    }

    m_mesh.garbage_collection();

    if (m_mesh.has_vertex_status()) m_mesh.release_vertex_status();
    if (m_mesh.has_face_status())   m_mesh.release_face_status();
    if (m_mesh.has_edge_status())   m_mesh.release_edge_status();
}

void Simplification::SimplifyFaceTo(int iRemainedFaceNum, double dAgreesiveness /*= 7*/)
{
    // because in every decimate:
    // remove one vertex, will remove three edges, and remove two triangles
    int iRemovedFaceNum = m_mesh.n_faces() - iRemainedFaceNum;
    int iRemovedVertexNum = iRemovedFaceNum / 2;
    int iRemainedVertexNum = m_mesh.n_vertices() - iRemovedVertexNum;
    SimplifyVertexTo(iRemainedVertexNum, dAgreesiveness);
}

void Simplification::Initialize()
{
    if (!m_vpQuadrics.is_valid())
    {
        m_mesh.add_property(m_vpQuadrics);
    }

    for (auto vIt = m_mesh.vertices_begin(); vIt != m_mesh.vertices_end(); ++vIt)
    {
        m_mesh.property(m_vpQuadrics, *vIt).Clear();
    }

    for (auto fIt = m_mesh.faces_begin(); fIt != m_mesh.faces_end(); ++fIt)
    {
        m_mesh.property(m_fpError, *fIt) = FLT_MAX;
    }

    MyMesh::FaceVertexIter fvIt;
    MyMesh::VertexHandle vh0, vh1, vh2;
    for (auto fIt = m_mesh.faces_begin(); fIt != m_mesh.faces_end(); ++fIt)
    {
        fvIt = m_mesh.fv_iter(*fIt);
        vh0 = *fvIt; ++fvIt;
        vh1 = *fvIt; ++fvIt;
        vh2 = *fvIt;

        OpenMesh::Vec3d v0 = OpenMesh::vector_cast<OpenMesh::Vec3d>(m_mesh.point(vh0));
        OpenMesh::Vec3d v1 = OpenMesh::vector_cast<OpenMesh::Vec3d>(m_mesh.point(vh1));
        OpenMesh::Vec3d v2 = OpenMesh::vector_cast<OpenMesh::Vec3d>(m_mesh.point(vh2));

        OpenMesh::Vec3d n = (v1 - v0) % (v2 - v0);
        double area = n.norm();
        if (area > FLT_MIN)
        {
            n /= area;
        }

        const double a = n[0];
        const double b = n[1];
        const double c = n[2];
        const double d = -(OpenMesh::vector_cast<OpenMesh::Vec3d>(m_mesh.point(vh0)) | n);

        SymetricMatrix q(a, b, c, d);

        m_mesh.property(m_vpQuadrics, vh0) += q;
        m_mesh.property(m_vpQuadrics, vh1) += q;
        m_mesh.property(m_vpQuadrics, vh2) += q;
    }

    double dError = 0;
    MyMesh::Point ptResult;
    for (auto eIt = m_mesh.edges_sbegin(); eIt != m_mesh.edges_end(); ++eIt)
    {
        dError = CalculateError(*eIt, ptResult);
        m_mesh.property(m_epError, *eIt) = dError;
        m_mesh.property(m_epTargetPoints, *eIt) = ptResult;
        auto fh0 = m_mesh.face_handle(m_mesh.halfedge_handle(*eIt, 0));
        auto fh1 = m_mesh.face_handle(m_mesh.halfedge_handle(*eIt, 1));

        if (m_mesh.is_valid_handle(fh0))
        {
            m_mesh.property(m_fpError, fh0) = std::min(m_mesh.property(m_fpError, fh0), dError);
        }
        if (m_mesh.is_valid_handle(fh1))
        {
            m_mesh.property(m_fpError, fh1) = std::min(m_mesh.property(m_fpError, fh1), dError);
        }
    }
}

double Simplification::CalculateError(OpenMesh::EdgeHandle edge, MyMesh::Point& ptResult)
{

    OpenMesh::HalfedgeHandle h0 = m_mesh.halfedge_handle(edge, 0);
    OpenMesh::VertexHandle v0 = m_mesh.from_vertex_handle(h0);
    OpenMesh::VertexHandle v1 = m_mesh.to_vertex_handle(h0);
    SymetricMatrix q = m_mesh.property(m_vpQuadrics, v0) + m_mesh.property(m_vpQuadrics, v1);
    double dError = 0;
    double det = q.Det(0, 1, 2, 1, 4, 5, 2, 5, 7);
    if (det != 0 && !m_mesh.is_boundary(edge))
    {
        ptResult[0] = -1 / det * (q.Det(1, 2, 3, 4, 5, 6, 5, 7, 8));	// vx = A41/det(q_delta)
        ptResult[1] = 1 / det * (q.Det(0, 2, 3, 1, 5, 6, 2, 7, 8));	    // vy = A42/det(q_delta)
        ptResult[2] = -1 / det * (q.Det(0, 1, 3, 1, 4, 6, 2, 5, 8));	// vz = A43/det(q_delta)

        dError = VertexError(q, ptResult[0], ptResult[1], ptResult[2]);
    }
    else
    {
        // find if v0, or v1, or midpoint
        MyMesh::Point pt0 = m_mesh.point(v0);
        MyMesh::Point pt1 = m_mesh.point(v1);
        MyMesh::Point pt2 = (pt0 + pt1) / 2;
        double error0 = VertexError(q, pt0[0], pt0[1], pt0[2]);
        double error1 = VertexError(q, pt1[0], pt1[1], pt1[2]);
        double error2 = VertexError(q, pt2[0], pt2[1], pt2[2]);
        dError = std::min(error0, std::min(error1, error2));
        if (error0 == dError)
        {
            ptResult = pt0;
        }
        if (error1 == dError)
        {
            ptResult = pt1;
        }
        if (error2 == dError)
        {
            ptResult = pt2;
        }
    }

    return dError;
}

double Simplification::VertexError(const SymetricMatrix& q, double x, double y, double z)
{
    return   q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x + q[4] * y * y
        + 2 * q[5] * y * z + 2 * q[6] * y + q[7] * z * z + 2 * q[8] * z + q[9];
}

// 基本判断逻辑：
// 不考虑即将被删除的三角形
// 如果变化后的三角形，有近似相同的边，则返回true
// 如果变化后的三角形之间的法向量之间的夹角过大，则返回true
bool Simplification::IsFlipped(OpenMesh::VertexHandle v0, OpenMesh::VertexHandle v1, const MyMesh::Point& ptTarget)
{
    for (auto fIt = m_mesh.vf_iter(v0); fIt.is_valid(); ++fIt)
    {
        if (m_mesh.status(*fIt).deleted())
        {
            continue;
        }

        auto vIt = m_mesh.fv_iter(*fIt);
        OpenMesh::VertexHandle fv[3];
        fv[0] = *vIt; ++vIt;
        fv[1] = *vIt; ++vIt;
        fv[2] = *vIt;

        // ignore the face will be deleted
        if (fv[0] == v1 || fv[1] == v1 || fv[2] == v1)
        {
            continue;
        }

        int idxV0 = 0;
        for (int i = 0; i < 3; ++i)
        {
            if (fv[i] == v0)
            {
                idxV0 = i;
            }
        }

        MyMesh::Point pt0 = m_mesh.point(v0);
        MyMesh::Point pt1 = m_mesh.point(fv[(idxV0 + 1) % 3]);
        MyMesh::Point pt2 = m_mesh.point(fv[(idxV0 + 2) % 3]);

        MyMesh::Point dir1 = pt1 - ptTarget;
        dir1.normalize();
        MyMesh::Point dir2 = pt2 - ptTarget;
        dir2.normalize();

        if (fabs(dir1.dot(dir2)) > 0.999) return true;

        MyMesh::Point normold = m_mesh.normal(*fIt);
        MyMesh::Point norm = dir1.cross(dir2);
        norm.normalize();

        if (normold.dot(norm) < 0.2) return true;
    }
    return false;
}

void Simplification::UpdateFaceNormal(OpenMesh::VertexHandle v0)
{
    if (!m_mesh.is_valid_handle(v0) || m_mesh.is_isolated(v0))
    {
        return;
    }
    for (auto fIt = m_mesh.vf_iter(v0); fIt.is_valid(); ++fIt)
    {
        if (!m_mesh.status(*fIt).deleted())
        {
            m_mesh.set_normal(*fIt, m_mesh.calc_face_normal(*fIt));
            m_mesh.property(m_fpDirty, *fIt) = true;
        }
    }
}

void Simplification::UpdateEdgePropertyAroundV(OpenMesh::VertexHandle v0)
{
    if (!m_mesh.is_valid_handle(v0) || m_mesh.is_isolated(v0))
    {
        return;
    }
    double dError = 0;
    MyMesh::Point ptResult;
    for (auto hIt = m_mesh.voh_iter(v0); hIt.is_valid(); ++hIt)
    {
        OpenMesh::EdgeHandle eh = m_mesh.edge_handle(*hIt);

        if (!m_mesh.is_valid_handle(*hIt) || m_mesh.status(eh).deleted())
        {
            continue;
        }

        dError = CalculateError(eh, ptResult);
        m_mesh.property(m_epError, eh) = dError;
        m_mesh.property(m_epTargetPoints, eh) = ptResult;
    }
}
