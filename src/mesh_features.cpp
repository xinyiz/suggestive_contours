#include "mesh_features.h"
#include <iostream>
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
  Vector3d n1 = computeFaceNormal(mesh, fh1);
  Vector3d n2 = computeFaceNormal(mesh, fh2);
  Vector3d viewVec = Vector3d(cameraPos[0], cameraPos[1], cameraPos[2]);

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
  Vector3d n1 = computeFaceNormal(mesh, fh1);
  Vector3d n2 = computeFaceNormal(mesh, fh2);

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
