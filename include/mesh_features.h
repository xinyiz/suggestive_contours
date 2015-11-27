#ifndef MESH_FEATURES_H
#define MESH_FEATURES_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "mesh_definitions.h"

Eigen::Vector3d computeFaceNormal(Mesh &mesh, const OpenMesh::FaceHandle &fh);
bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, OpenMesh::Vec3f cameraPos);
bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e);
bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, OpenMesh::Vec3f cameraPos);

#endif
