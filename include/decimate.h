#define DECIMATE_HH

//== INCLUDES =================================================================

#include "mesh_definitions.h"

void simplify(Mesh& _mesh, const float _percentage, const std::string _output_filename);


//== CLASS/STRUCT DEFINITION =========================================================

/** /struct VertexPriority
Stores a vertex handle, minimum priority, the corresponding halfedge handle, and version.
The priority is valid when the version is the latest one.
**/

struct VertexPriority {
  VertexPriority(const Mesh::VertexHandle _vh, const Mesh::HalfedgeHandle _heh,
    const double _priority, const int _version)
    : vh_(_vh), heh_(_heh), priority_(_priority), version_(_version) { }

  VertexPriority(const VertexPriority& _other)
    : vh_(_other.vh_), heh_(_other.heh_), priority_(_other.priority_), version_(_other.version_) { }

  Mesh::VertexHandle vh_;
  Mesh::HalfedgeHandle heh_;
  double priority_;
  int version_;
};

struct VertexPriorityCompare {
  bool operator()(const VertexPriority& _vp_0, const VertexPriority& _vp_1) const {
    return _vp_0.priority_ > _vp_1.priority_;
  }
};
