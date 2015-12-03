#include "mesh_features.h"
#include <iostream>
#include <math.h>
using namespace OpenMesh;
using namespace Eigen;

Vector3d computeFaceNormal(Mesh &mesh, const OpenMesh::FaceHandle &fh){

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

  return E1.cross(E2).normalized();
}
bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {

  FaceHandle fh1 = mesh.face_handle(mesh.halfedge_handle(e,0));
  FaceHandle fh2 = mesh.face_handle(mesh.halfedge_handle(e,1));
  
  //First face normal
  const Vec3f N1 = mesh.normal(fh1);
  const Vec3f N2 = mesh.normal(fh2);
  Vector3d n1 = Vector3d(N1[0], N1[1], N1[2]).normalized();
  Vector3d n2 = Vector3d(N2[0], N2[1], N2[2]).normalized();

  VertexHandle vih = mesh.to_vertex_handle(mesh.halfedge_handle(e,0));
  Vec3f p = mesh.point(vih);
  Vector3d point(p[0], p[1], p[2]);

  Vector3d viewVec = Vector3d(cameraPos[0], cameraPos[1], cameraPos[2]) - point;

  if((n1.dot(viewVec))*(n2.dot(viewVec)) < 0){
    return true;
  }
  return false;
  // CHECK IF e IS A SILHOUETTE HERE -----------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {

  FaceHandle fh1 = mesh.face_handle(mesh.halfedge_handle(e,0));
  FaceHandle fh2 = mesh.face_handle(mesh.halfedge_handle(e,1));
  
  //First face normal
  const Vec3f N1 = mesh.normal(fh1);
  const Vec3f N2 = mesh.normal(fh2);
  Vector3d n1 = Vector3d(N1[0], N1[1], N1[2]).normalized();
  Vector3d n2 = Vector3d(N2[0], N2[1], N2[2]).normalized();

  if(n1.dot(n2) < 0.5){
    return true;
  }
  return false;
  // CHECK IF e IS SHARP HERE ------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
  return mesh.is_boundary(e) || isSilhouette(mesh, e, cameraPos) || isSharpEdge(mesh, e);
}
