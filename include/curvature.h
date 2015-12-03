#ifndef CURVATURE_H
#define CURVATURE_H

#include "mesh_definitions.h"

struct CurvatureInfo {
  OpenMesh::Vec3f directions[2];
  double curvatures[2];
};

float computeFaceArea(const Mesh &mesh, const Mesh::FaceHandle &fh);
float computeCurvature(Mesh& mesh, OpenMesh::VPropHandleT<CurvatureInfo>& curvature, float max_curvature);
void computeViewCurvature(Mesh& mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo>& curvature,
  OpenMesh::VPropHandleT<double>& viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative);

#endif
