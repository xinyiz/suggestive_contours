#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/Core>
#include <Eigen/QR>

#include "decimate.h"
#include <iostream>
#include <limits>
#include <queue>


using namespace OpenMesh;
using namespace Eigen;

// NOTE: We do not use pointers in the queue to avoid any memory leak issues.
// When it is guaranteed that everyone can use c+11 functionality,
// we can consider to use std::unique_ptr<> or std::shared_ptr<> later.
typedef std::priority_queue<VertexPriority, std::vector<VertexPriority>,
  VertexPriorityCompare> VertexPriorityQueue;
VPropHandleT<Eigen::Matrix4d> vprop_quadric;
VPropHandleT<int> vprop_latest_version;


// Mesh property accessors
Eigen::Matrix4d& vertex_quadric(Mesh& _mesh, const Mesh::VertexHandle _vh) {
  return _mesh.property(vprop_quadric, _vh);
}

// The vertex version is used for tracking the latest out-going halfedge with
// the minimum priority
int& vertex_latest_version(Mesh& _mesh, const Mesh::VertexHandle _vh) {
  return _mesh.property(vprop_latest_version, _vh);
}


// Functions
void intialize(Mesh& _mesh);
double compute_priority(Mesh& _mesh, const Mesh::HalfedgeHandle _heh);
bool is_collapse_valid(Mesh& _mesh, const Mesh::HalfedgeHandle _heh);
bool is_vertex_priority_valid(Mesh& _mesh, const VertexPriority& _vp);
void enqueue_vertex(Mesh& _mesh, VertexPriorityQueue& _queue, const Mesh::VertexHandle _vh);
void decimate(Mesh& _mesh, const unsigned int _target_num_vertices);


void simplify(Mesh& _mesh, const float _percentage, const std::string _output_filename) {
  // Add required properties
  _mesh.request_vertex_status();
  _mesh.request_edge_status();
  _mesh.request_face_status();
  _mesh.request_face_normals();
  _mesh.add_property(vprop_quadric);
  _mesh.add_property(vprop_latest_version);

  // Compute normals and quadrics
  intialize(_mesh);

  // Decimate
  decimate(_mesh, (int)(_percentage * _mesh.n_vertices()));
  std::cout << "Simplified to #vertices: " << _mesh.n_vertices() << std::endl;

  _mesh.remove_property(vprop_quadric);
  _mesh.remove_property(vprop_latest_version);

  // Write to file
  IO::Options opt;
  std::cout << "Writing to file '" << _output_filename << "'... ";
  if (!IO::write_mesh(_mesh, _output_filename, opt)) {
    std::cout << "Failed!" << std::endl;;
  }
  std::cout << "Done." << std::endl;;
}

void intialize(Mesh& _mesh) {
  // Compute face normals
  _mesh.update_face_normals();

  for (Mesh::ConstVertexIter v_it = _mesh.vertices_begin();
    v_it != _mesh.vertices_end(); ++v_it) {
    const Mesh::VertexHandle vh = (*v_it);
    vertex_quadric(_mesh, vh).setZero();
    vertex_latest_version(_mesh, vh) = 0;

    // INSERT CODE HERE FOR PART 1-------------------------------------------------------------------------------
    // Calculate vertex quadrics from incident triangles
    // ----------------------------------------------------------------------------------------------------------
    for (Mesh::ConstVertexFaceIter f_it = _mesh.cvf_begin(vh); f_it != _mesh.cvf_end(vh); ++f_it) {
      const auto fh = *f_it;
      Eigen::Vector4d fcoeff;
      // calculate face normal
      const auto normal = _mesh.calc_face_normal(fh).normalized();

      // calculate face plane coefficients
      const auto point = _mesh.point(vh);
      auto d = dot(point, normal);
      fcoeff << normal[0], normal[1], normal[2], -d;

      // add face quadric contribution to the vertex quadric
      vertex_quadric(_mesh, vh) += fcoeff * fcoeff.transpose();
    }
  }

  std::cout << "Finished initialization." << std::endl;
}

// calculate optimal position for a given pair of vertices and associated quadric
Vector4d optimal_pos(const Mesh& _mesh, const Matrix4d& Q, const VertexHandle& v1, const VertexHandle& v2) {
  Vector4d ret;
  Matrix4d s_mat = Q;
  s_mat.row(3) << 0., 0., 0., 1.;
  Eigen::ColPivHouseholderQR<Matrix4d> qr;
  qr.setThreshold(1.1);
  qr.compute(s_mat);
  // is s_mat full rank?
  if (qr.rank() == 4) {
    Vector4d rhs;
    rhs << 0., 0., 0., 1.;
    ret = qr.solve(rhs);
    ret /= ret[3]; // normalize to w = 1
    std::cout << "Hello" << std::endl;
  } else {
    // handle degenerate cases: naive implementation is to return mid point instead
    const auto l_pt = _mesh.point(v1);
    const auto m_pt = 0.5 * (_mesh.point(v1) + _mesh.point(v2));
    const auto r_pt = _mesh.point(v2);
    ret << r_pt[0], r_pt[1], r_pt[2], 1.0;
  }
  return ret;
}

double compute_priority(Mesh& _mesh, const Mesh::HalfedgeHandle _heh) {
  double priority = 0.0;

  // INSERT CODE HERE FOR PART 2---------------------------------------------------------------------------------
  // Return priority: The smaller the better
  // Use quadrics to estimate approximation error
  // -------------------------------------------------------------------------------------------------------------
  const auto vh_f = _mesh.from_vertex_handle(_heh);
  const auto vh_t = _mesh.to_vertex_handle(_heh);
  // calculate poor man's quadric
  const auto Q = vertex_quadric(_mesh, vh_f) + vertex_quadric(_mesh, vh_t);
  // calculate optimal position to reduce the quadric error
  const auto new_pos = optimal_pos(_mesh, Q, vh_f, vh_t);
  // calculate associated error
  priority = new_pos.transpose() * Q * new_pos;

  return priority;
}

bool is_collapse_valid(Mesh& _mesh, const Mesh::HalfedgeHandle _heh) {
  const Mesh::VertexHandle from_vh = _mesh.from_vertex_handle(_heh);
  const Mesh::VertexHandle to_vh = _mesh.to_vertex_handle(_heh);

  // Collect faces
  const Mesh::FaceHandle fh_0 = _mesh.face_handle(_heh);
  const Mesh::FaceHandle fh_1 = _mesh.face_handle(_mesh.opposite_halfedge_handle(_heh));

  // Backup point positions
  const Mesh::Point from_p = _mesh.point(from_vh);
  const Mesh::Point to_p = _mesh.point(to_vh);

  // Topological test
  if (!_mesh.is_collapse_ok(_heh))
    return false;

  // Test boundary
  if (_mesh.is_boundary(from_vh) && !_mesh.is_boundary(to_vh))
    return false;

  // Test for normal flipping
  for (Mesh::ConstVertexFaceIter vf_it = _mesh.cvf_begin(from_vh); vf_it != _mesh.cvf_end(from_vh); ++vf_it) {
    const Mesh::FaceHandle n_fh = (*vf_it);
    if (fh_0 == n_fh || fh_1 == n_fh) continue;
    const Mesh::Normal n_before = _mesh.normal(n_fh).normalized();

    Vec3f nf_p[3];
    Mesh::ConstFaceVertexIter n_fv_it = _mesh.cfv_begin(n_fh);
    for (int i = 0; n_fv_it != _mesh.cfv_end(n_fh) && i < 3; ++n_fv_it, ++i) {
      const Mesh::VertexHandle nn_vh = (*n_fv_it);
      nf_p[i] = _mesh.point(nn_vh);

      // Replace 'from' point to 'to' point.
      if (nf_p[0] == from_p) nf_p[0] = to_p;
    }

    const Mesh::Normal cross_prod = cross(nf_p[1] - nf_p[0], nf_p[2] - nf_p[0]);

    if (std::abs(cross_prod.norm()) > 1.0E-8) {
      const Mesh::Normal n_after = cross_prod.normalized();

      // Consider the triangle is flipped if the normal angle is changed more than 45 degrees
      const Mesh::Scalar cos_pi_over_4 = 1 / sqrt(2.0);
      if (dot(n_before, n_after) < cos_pi_over_4)
        return false;
    }
  }

  // Collapse passed all tests
  return true;
}

bool is_vertex_priority_valid(Mesh& _mesh, const VertexPriority& _vp) {
  // The halfedge priority is valid only when its version is equal to the
  // 'from' vertex version.
  return (_vp.version_ == vertex_latest_version(_mesh, _vp.vh_));
}

void enqueue_vertex(Mesh& _mesh, VertexPriorityQueue& _queue,
  const Mesh::VertexHandle _vh) {
  double min_priority = std::numeric_limits<double>::max();
  Mesh::HalfedgeHandle min_heh;

  // Find the minimum priority out-going halfedge
  for (Mesh::ConstVertexOHalfedgeIter vh_it = _mesh.cvoh_begin(_vh); vh_it != _mesh.cvoh_end(_vh); ++vh_it) {
    if (is_collapse_valid(_mesh, *vh_it)) {
      const double priority = compute_priority(_mesh, *vh_it);
      if (priority < min_priority) {
        min_priority = priority;
        min_heh = (*vh_it);
      }
    }
  }

  // Update queue
  if (min_priority < std::numeric_limits<double>::max()) {
    // Increase the vertex version and use the updated version for the halfedge
    int& version = vertex_latest_version(_mesh, _vh);
    ++version;
    _queue.emplace(_vh, min_heh, min_priority, version);
  }
}

void decimate(Mesh& _mesh, const unsigned int _target_num_vertices) {
  std::cout << "Starting decimation... ";

  // Build priority queue
  VertexPriorityQueue queue;
  for (Mesh::ConstVertexIter v_it = _mesh.vertices_begin();
    v_it != _mesh.vertices_end(); ++v_it) {
    const Mesh::VertexHandle vh = (*v_it);
    enqueue_vertex(_mesh, queue, vh);
  }

  int num_vertices = _mesh.n_vertices();

  // INSERT CODE HERE FOR PART 3-----------------------------------------------------------------------------------
  // Decimate using priority queue:
  //
  uint decimate_count = 0;
  while (num_vertices > _target_num_vertices && !queue.empty()) {
    // 1) Take first element of queue
    const auto vp = queue.top();
    queue.pop();

    //  - Check whether the vertex priority is valid using 'is_vertex_priority_valid()'
    if (!is_vertex_priority_valid(_mesh, vp)) { continue;}

    // 2) Collapse this halfedge
    const auto vheh = vp.heh_;
    //  - Check whether the halfedge collapse is valid using 'is_collapse_valid()'
    if (!is_collapse_valid(_mesh, vheh)) { continue;}

    // get new vertex location
    const auto vh_f = _mesh.from_vertex_handle(vheh);
    const auto vh_t = _mesh.to_vertex_handle(vheh);
    const auto Q = vertex_quadric(_mesh, vh_f) + vertex_quadric(_mesh, vh_t);
    const auto new_pos = optimal_pos(_mesh, Q, vh_f, vh_t);

    // update relations to be deleted for the to_vertex
    _mesh.collapse(vheh); //this collapses the from vertex to the to vertex.

    // change point position
    // _mesh.set_point(vh_t, Vec3f(new_pos[0], new_pos[1], new_pos[2]));

    // 3) Update queue
    // enqueue all neighboring vertices including this one
    enqueue_vertex(_mesh, queue, vh_t);
    for (Mesh::ConstVertexVertexIter vv_it = _mesh.vv_begin(vh_t); vv_it != _mesh.vv_end(vh_t); ++vv_it) {
      enqueue_vertex(_mesh, queue, *vv_it);
    }

    // Delete the items marked to be deleted
    _mesh.garbage_collection();

    // update num vertices
    num_vertices = _mesh.n_vertices();
  }
  // --------------------------------------------------------------------------------------------------------------

  std::cout << "Done." << std::endl;
}
