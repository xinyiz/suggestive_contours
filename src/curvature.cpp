#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <math.h>
#include "curvature.h"

using namespace OpenMesh;
using namespace Eigen;
using namespace std;

float computeFaceArea(const Mesh &mesh, const Mesh::FaceHandle &fh){

  Vec3f point[3];
  Mesh::ConstFaceVertexIter cfv_it;
  cfv_it = mesh.cfv_iter(fh);
  point[0] = mesh.point(*cfv_it);
  point[1] = mesh.point(*(++cfv_it));
  point[2] = mesh.point(*(++cfv_it));
  Vec3f e1 = point[1] - point[0];
  Vec3f e2 = point[2] - point[0];
  Vector3d E1(e1[0], e1[1], e1[2]);
  Vector3d E2(e2[0], e2[1], e2[2]);

  return 0.5*E1.cross(E2).norm();
}
void computeCurvature(Mesh& mesh, OpenMesh::VPropHandleT<CurvatureInfo>& curvature) {
  std::vector<Vector3d> Ts; 
  std::vector<float> ks; 
  std::vector<float> ws; 
  float total_w;

  for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {

    // WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
    const Mesh::VertexHandle vih = (*v_it);
    const Vec3f normal = mesh.normal(vih);
    const Vec3f p = mesh.point(vih);

    const Vector3d  N_vi = Vector3d(normal[0], normal[1], normal[2]).normalized();
    const Vector3d vi = Vector3d(p[0],p[1],p[2]);

    Ts.clear(); 
    ks.clear(); 
    ws.clear(); 
    total_w = 0;

    // Iterate over all outgoing halfedges...
    for(Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(vih); voh_it; ++voh_it) {

      const Mesh::HalfedgeHandle vi_vjh = (*voh_it);
      const Mesh::VertexHandle vjh = mesh.to_vertex_handle(vi_vjh);
      const Vec3f vjp = mesh.point(vjh);

      // Calculate T_ij
      const Vector3d vj = Vector3d(vjp[0], vjp[1], vjp[2]);
     
      // vi - vj
      const Vector3d delta = (vi - vj);
      
      const Matrix3d I = Matrix<double,3,3>::Identity();
      Vector3d T_ij = ((I - N_vi*N_vi.transpose())*delta).normalized();
      Ts.push_back(T_ij);

      // Calculate k_ij
      double k_ij = -2*N_vi.transpose()*delta;
      k_ij = k_ij/(delta.norm()*delta.norm());
      ks.push_back(k_ij);

      // Calculate w_ij
      const Mesh::FaceHandle f1 = mesh.face_handle(vi_vjh);
      const Mesh::FaceHandle f2 = mesh.opposite_face_handle(vi_vjh);
      float weight = computeFaceArea(mesh,f1) + computeFaceArea(mesh,f2);
      ws.push_back(weight);
      total_w+=weight;

    }

    // Normalize weights

    for(unsigned int i; i < ws.size(); i++){
      ws[i] = ws[i]/total_w;
    }

    // Calculate M
    Matrix3d M = Matrix<double,3,3>::Zero();
    for(unsigned int i = 0; i < ws.size(); i++){
      Matrix3d TT = Ts[i]*Ts[i].transpose();
      M = M + ws[i]*ks[i]*TT;
    }
    SelfAdjointEigenSolver<Matrix3d> es(M);
    Matrix3d evecs = es.eigenvectors();

    double evals0 = abs(es.eigenvalues()[0]);
    double evals1 = abs(es.eigenvalues()[1]);
    double evals2 = abs(es.eigenvalues()[2]);
    unsigned i1;
    unsigned i2;
    unsigned min_i;
    if( (evals0 < evals1) && (evals0 < evals2)){
      i1 = 1;
      i2 = 2;
      min_i = 0;
    } else if ((evals1 < evals2) && (evals1 < evals0)) {
      i1 = 0;
      i2 = 2;
      min_i = 1;
    } else {
      i1 = 0;
      i2 = 1;
      min_i = 2;
    }
    // In the end you need to fill in this struct
    CurvatureInfo info;
    info.curvatures[0] = es.eigenvalues()[i1];
    info.curvatures[1] = es.eigenvalues()[i2];
    Vector3d ev1 = es.eigenvectors().col(i1).normalized();
    Vector3d ev2 = es.eigenvectors().col(i2).normalized();
    info.directions[0] = Vec3f(ev1[0], ev1[1], ev1[2]);
    info.directions[1] = Vec3f(ev2[0], ev2[1], ev2[2]);

    mesh.property(curvature, vih) = info;
    // -------------------------------------------------------------------------------------------------------------
  }
}

void computeViewCurvature(Mesh& mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo>& curvature,
  OpenMesh::VPropHandleT<double>& viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative) {
  // WRITE CODE HERE TO COMPUTE CURVATURE IN THE VIEW PROJECTION PROJECTED ON THE TANGENT PLANE ------------------
  // Compute vector to viewer and project onto tangent plane, then use components in principal directions to find curvature
  // -------------------------------------------------------------------------------------------------------------
  for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {

    const Mesh::VertexHandle vih = (*v_it);
    const Vec3f normal = mesh.normal(vih);
    const Vec3f p = mesh.point(vih);

    const Vector3d  N_vi = Vector3d(normal[0], normal[1], normal[2]).normalized();
    const Vector3d vi = Vector3d(p[0],p[1],p[2]);

    CurvatureInfo info = mesh.property(curvature, vih);
    float K1 = info.curvatures[0];
    float K2 = info.curvatures[1];
    Vector3d T1 = Vector3d(info.directions[0][0],info.directions[0][1],info.directions[0][2]);
    Vector3d T2 = Vector3d(info.directions[1][0],info.directions[1][1],info.directions[1][2]);

    // Project view to tangent plane
    Vector3d viewVec = Vector3d(camPos[0],camPos[1],camPos[2]) - vi;
    Vector3d viewTanDir = (viewVec - N_vi.dot(viewVec)*viewVec).normalized();
    float cosTheta2 = pow(viewTanDir.dot(T1),2.0);
    float sinTheta2 = 1.0 - cosTheta2;
    float c_view = K1*cosTheta2 + K2*sinTheta2;


    mesh.property(viewCurvature, vih) = c_view;
    // -------------------------------------------------------------------------------------------------------------
  }

  // We'll use the finite elements piecewise hat method to find per-face gradients of the view curvature
  // CS 348a doesn't cover how to differentiate functions on a mesh (Take CS 468!) so we provide code here

  for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
    const Mesh::FaceHandle fh = (*f_it);

    Vec3f p[3];
    double c[3];
    Mesh::ConstFaceVertexIter fv_it = mesh.cfv_begin(fh);
    for (int i = 0; i < 3 && fv_it != mesh.cfv_end(fh); ++i, ++fv_it) {
      const Mesh::VertexHandle n_vh = (*fv_it);
      p[i] = mesh.point(n_vh);
      c[i] = mesh.property(viewCurvature, n_vh);
    }

    const Vec3f n = mesh.normal(fh);
    double area = mesh.calc_sector_area(mesh.halfedge_handle(fh));

    mesh.property(viewCurvatureDerivative, fh) =
      cross(n, (p[0] - p[2])) * (c[1] - c[0]) / (2 * area) +
      cross(n, (p[1] - p[0])) * (c[2] - c[0]) / (2 * area);
  }
}
