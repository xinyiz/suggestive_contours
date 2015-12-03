#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <GLUT/glut.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
#include "mesh_features.h"
#include "image_generation.h"
#include "decimate.h"

using namespace std;
using namespace OpenMesh;
using namespace Eigen;

VPropHandleT<double> viewCurvature;
FPropHandleT<Vec3f> viewCurvatureDerivative;
VPropHandleT<CurvatureInfo> curvature;
Mesh mesh;


bool leftDown = false, rightDown = false, middleDown = false;
int lastPos[2];
float cameraPos[4] = { 0, 0, 4, 1 };
Vec3f up, pan;
int windowWidth = 640, windowHeight = 480;
bool showSurface = true, showAxes = false, showCurvature = false, showNormals = false, renderFace = false;

float specular[] = { 1.0, 1.0, 1.0, 1.0 };
float shininess[] = { 50.0 };

void renderSuggestiveContours(Vec3f actualCamPos) { // use this camera position to account for panning etc.
  glBegin(GL_LINES);
  glColor3f(0.8, 0.8, 0.5);
  glLineWidth(2.0f);
  for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
    const Mesh::FaceHandle fh = (*f_it);

    Vec3f p[3];
    double c[3];
    unsigned int zero_crossing[3] = {0};
    Vec3f vcd = mesh.property(viewCurvatureDerivative, fh);
    Vector3d viewCurveDeriv = Vector3d(vcd[0], vcd[1], vcd[2]);
    const Vec3f n = mesh.normal(fh);
    Vector3d N = Vector3d(n[0], n[1], n[2]).normalized();


    Mesh::ConstFaceVertexIter fv_it = mesh.cfv_begin(fh);
    for (int i = 0; i < 3 && fv_it != mesh.cfv_end(fh); ++i, ++fv_it) {
      const Mesh::VertexHandle n_vh = (*fv_it);
      p[i] = mesh.point(n_vh);
      c[i] = mesh.property(viewCurvature, n_vh);
    }
     
    // check if there is zero crossing
    if(c[0]*c[1] < 0){
      zero_crossing[0] = 1;
    }
    if (c[1]*c[2] < 0){
      zero_crossing[1] = 1;
    }
    if (c[2]*c[0] < 0){
      zero_crossing[2] = 1;
    }
      // check viewCurvatureDerivative
    if((zero_crossing[0] + zero_crossing[1] + zero_crossing[2]) == 2){
      const Vector3d vi1 = Vector3d(p[0][0],p[0][1],p[0][2]);
      const Vector3d vi2 = Vector3d(p[1][0],p[1][1],p[1][2]);
      const Vector3d vi3 = Vector3d(p[2][0],p[2][1],p[2][2]);
      const Vector3d vi = (vi1 + vi2 + vi3)/3.0;
      Vector3d viewVec = Vector3d(actualCamPos[0],actualCamPos[1],actualCamPos[2])- vi;
      float DK = (viewVec.normalized()).dot(viewCurveDeriv);
      float DK_thresh = DK/viewVec.norm();
      float thetad = acos(N.dot(viewVec.normalized()));
      if(DK_thresh > 0.02 && thetad > 0.2){
        if(zero_crossing[0]){
          glVertex3f((p[0][0]+p[1][0])/2.0, 
                     (p[0][1]+p[1][1])/2.0,
                     (p[0][2]+p[1][2])/2.0);
        }
        if(zero_crossing[1]){
          glVertex3f((p[1][0]+p[2][0])/2.0, 
                     (p[1][1]+p[2][1])/2.0,
                     (p[1][2]+p[2][2])/2.0);
        }
        if(zero_crossing[2]){
          glVertex3f((p[0][0]+p[2][0])/2.0, 
                     (p[0][1]+p[2][1])/2.0,
                     (p[0][2]+p[2][2])/2.0);
        }
      }
    }

  }
  glEnd();

  // RENDER SUGGESTIVE CONTOURS HERE -----------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------
}

/* Compute Surface Color */
Vec3f calculateColor(CurvatureInfo c_info) {
	float c_curv1 = c_info.curvatures[0];
	float c_curv2 = c_info.curvatures[1];
	//Normalize the values so that they are between 0 and 1
	c_curv1 = (c_curv1 + 0.05f) / 0.1f;
	c_curv2 = (c_curv2 + 0.05f) / 0.1f;
	float c_curv;
	//Find the max of the two curvatures
	if (c_curv1 > c_curv2) c_curv = c_curv1;
	else c_curv = c_curv2;
	Vec3f color;
	//Choose color based on curvature (red = valley, blue = mountain)
	if (c_curv >= 0.75f) color = Vec3f(c_curv, 0, 0);
	if (c_curv < 0.75f && c_curv >= 0.5f) color = Vec3f(c_curv, 1 - c_curv, 0);
	if (c_curv < 0.5f && c_curv >= 0.25f) color = Vec3f(0, c_curv, 1 - c_curv);
	if (c_curv < 0.25f) color = Vec3f(0, 0, 1 - c_curv);
	return color;
}

void renderMesh() {
  if (!showSurface) glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // render regardless to remove hidden lines

  glDisable(GL_LIGHTING);
  glLightfv(GL_LIGHT0, GL_POSITION, cameraPos);

  glDepthRange(0.001, 1);
  glEnable(GL_NORMALIZE);

  // WRITE CODE HERE TO RENDER THE TRIANGLES OF THE MESH ---------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------
// show edges
  for (Mesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it){
    OpenMesh::Vec3f point[3];
	OpenMesh::Vec3f normal[3];
    Mesh::ConstFaceVertexIter cfv_it;
    cfv_it = mesh.cfv_iter(*f_it);
	const auto it0 = *cfv_it;
	const auto it1 = *(++cfv_it);
	const auto it2 = *(++cfv_it);
    point[0] = mesh.point(it0); normal[0] = mesh.normal(it0);
    point[1] = mesh.point(it1); normal[1] = mesh.normal(it1);
    point[2] = mesh.point(it2); normal[2] = mesh.normal(it2);

	/**************************Colors**************************/
    CurvatureInfo c_info1 = mesh.property(curvature, it0);
	CurvatureInfo c_info2 = mesh.property(curvature, it1);
	CurvatureInfo c_info3 = mesh.property(curvature, it2);
	Vec3f color1 = calculateColor(c_info1);
	Vec3f color2 = calculateColor(c_info2);
	Vec3f color3 = calculateColor(c_info3);
	/**********************************************************/

    glBegin(GL_TRIANGLES);
	glColor3f(color1[0], color1[1], color1[2]); 
	glNormal3f(normal[0][0], normal[0][1], normal[0][2]); 
	glVertex3f(point[0][0], point[0][1], point[0][2]);

	glColor3f(color3[0], color3[1], color3[2]); 
	glNormal3f(normal[2][0], normal[2][1], normal[2][2]); 
	glVertex3f(point[2][0], point[2][1], point[2][2]);

	glColor3f(color2[0], color2[1], color2[2]); 
	glNormal3f(normal[1][0], normal[1][1], normal[1][2]); 
	glVertex3f(point[1][0], point[1][1], point[1][2]);
    glEnd();
  }

  if (!showSurface) glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

  // show edges
  if (renderFace) {
    glDisable(GL_LIGHTING);
    glLineWidth(1.0f);
    glColor4f(0.0, 0.0, 0.0, 1.0);
    for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++ e_it) {
      OpenMesh::Vec3f point[2];
      const auto hh = mesh.halfedge_handle(*e_it, 1);
      point[0] = mesh.point(mesh.from_vertex_handle(hh));
      point[1] = mesh.point(mesh.to_vertex_handle(hh));
      glBegin(GL_LINES);
      glVertex3f(point[0][0], point[0][1], point[0][2]);
      glVertex3f(point[1][0], point[1][1], point[1][2]);
      glEnd();
    }
  }

  glDisable(GL_LIGHTING);
  glDepthRange(0, 0.999);

  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);
  renderSuggestiveContours(actualCamPos);

  // We'll be nice and provide you with code to render feature edges below
  glBegin(GL_LINES);
  glColor3f(0, 0, 0);
  glLineWidth(2.0f);
  for (Mesh::ConstEdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
    const Mesh::EdgeHandle eh = (*e_it);
    if (isFeatureEdge(mesh, eh, actualCamPos)) {
      const Mesh::HalfedgeHandle heh_0 = mesh.halfedge_handle(eh, 0);
      const Mesh::HalfedgeHandle heh_1 = mesh.halfedge_handle(eh, 1);
      const Vec3f source(mesh.point(mesh.from_vertex_handle(heh_0)));
      const Vec3f target(mesh.point(mesh.from_vertex_handle(heh_1)));
      glVertex3f(source[0], source[1], source[2]);
      glVertex3f(target[0], target[1], target[2]);
    }
  }
  glEnd();

  if (showCurvature) {
    // WRITE CODE HERE TO RENDER THE PRINCIPAL DIRECTIONS YOU COMPUTED ---------------------------------------------
    // -------------------------------------------------------------------------------------------------------------
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
      CurvatureInfo info = mesh.property(curvature,*v_it);
      const Mesh::VertexHandle vh = (*v_it);
      const Vec3f p = mesh.point(vh);
      const Vec3f dd1 = info.directions[0];
      const Vec3f dd2 = info.directions[1];
      const Vec3f d1 = p + dd1 * 0.03;
      const Vec3f d2 = p + dd2 * 0.03;
      glVertex3f(p[0], p[1], p[2]);
      glVertex3f(d1[0], d1[1], d1[2]);
      glVertex3f(p[0], p[1], p[2]);
      glVertex3f(d2[0], d2[1], d2[2]);
    }
    glEnd();
  }

  if (showNormals) {
    glBegin(GL_LINES);
    glColor3f(0, 1, 0);
    for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
      const Mesh::VertexHandle vh = (*v_it);
      const Vec3f n = mesh.normal(vh);
      const Vec3f p = mesh.point(vh);
      const Vec3f d = p + n * 0.01;
      glVertex3f(p[0], p[1], p[2]);
      glVertex3f(d[0], d[1], d[2]);
    }
    glEnd();
  }

  glDepthRange(0, 1);
}

void display() {
  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
  glEnable(GL_LIGHT0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0, 0, windowWidth, windowHeight);

  float ratio = (float)windowWidth / (float)windowHeight;
  gluPerspective(50, ratio, 1, 1000); // 50 degree vertical viewing angle, zNear = 1, zFar = 1000

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2],
    pan[0], pan[1], pan[2], up[0], up[1], up[2]);

  // Draw mesh
  renderMesh();

  // Draw axes
  if (showAxes) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glLineWidth(1);
    glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(1, 0, 0); // x axis
    glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 1, 0); // y axis
    glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 1); // z axis
    glEnd(/*GL_LINES*/);
  }

  glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) leftDown = (state == GLUT_DOWN);
  else if (button == GLUT_RIGHT_BUTTON) rightDown = (state == GLUT_DOWN);
  else if (button == GLUT_MIDDLE_BUTTON) middleDown = (state == GLUT_DOWN);

  lastPos[0] = x;
  lastPos[1] = y;
}

void mouseMoved(int x, int y) {
  const float speed = 30.0f;

  int dx = x - lastPos[0];
  int dy = y - lastPos[1];

  Vec3f curCamera(cameraPos[0], cameraPos[1], cameraPos[2]);
  Vec3f curCameraNormalized = curCamera.normalized();
  Vec3f right = up % curCameraNormalized;

  if (middleDown || (leftDown && rightDown)) {
    pan += -speed * (float)((float)dx / (float)windowWidth) * right +
      speed * (float)((float)dy / (float)windowHeight) * up;
  }
  else if (leftDown) {
    // Assume here that up vector is (0,1,0)
    Vec3f newPos = curCamera - speed * (float)((float)dx / (float)windowWidth) * right +
      speed * (float)((float)dy / (float)windowHeight) * up;
    newPos = newPos.normalized() * curCamera.length();

    up = up - (up | newPos) * newPos / newPos.sqrnorm();
    up.normalize();

    for (int i = 0; i < 3; i++) cameraPos[i] = newPos[i];
  }
  else if (rightDown) {
    for (int i = 0; i < 3; i++) cameraPos[i] *= pow(1.1, dy * 0.1);
  }


  lastPos[0] = x;
  lastPos[1] = y;

  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);
  computeViewCurvature(mesh, actualCamPos, curvature, viewCurvature, viewCurvatureDerivative);

  glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);

  if (key == 's' || key == 'S') showSurface = !showSurface;
  else if (key == 'a' || key == 'A') showAxes = !showAxes;
  else if (key == 'c' || key == 'C') showCurvature = !showCurvature;
  else if (key == 'n' || key == 'N') showNormals = !showNormals;
  else if (key == 'r' || key == 'R') renderFace = !renderFace;
  else if (key == 'd' || key == 'D') {
    float percentage = 1.0f;
    while (percentage <= 0.0f || percentage >= 1.0f) {
      std::cout << "Type percentage of vertices: ";
      std::cin >> percentage;
    }
    simplify(mesh, percentage, "output.off");
  }
  else if (key == 'w' || key == 'W') {
    writeImage(mesh, windowWidth, windowHeight, "renderedImage.svg", actualCamPos);
  }
  else if (key == 'q' || key == 'Q') exit(0);
  glutPostRedisplay();
}

void reshape(int width, int height) {
  windowWidth = width;
  windowHeight = height;
  glutPostRedisplay();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " mesh_filename\n";
    exit(0);
  }

  IO::Options opt;
  opt += IO::Options::VertexNormal;
  opt += IO::Options::FaceNormal;

  mesh.request_face_normals();
  mesh.request_vertex_normals();

  cout << "Reading from file " << argv[1] << "...\n";
  if (!IO::read_mesh(mesh, argv[1], opt)) {
    cout << "Read failed.\n";
    exit(0);
  }

  cout << "Mesh stats:\n";
  cout << '\t' << mesh.n_vertices() << " vertices.\n";
  cout << '\t' << mesh.n_edges() << " edges.\n";
  cout << '\t' << mesh.n_faces() << " faces.\n";

  mesh.update_normals();

  mesh.add_property(viewCurvature);
  mesh.add_property(viewCurvatureDerivative);
  mesh.add_property(curvature);

  // Move center of mass to origin
  Vec3f center(0, 0, 0);
  for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    center += mesh.point(*v_it);
  center /= mesh.n_vertices();

  for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    mesh.point(*v_it) -= center;

  // Fit in the unit sphere
  float maxLength = 0;
  for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    maxLength = max(maxLength, mesh.point(*v_it).length());

  if (maxLength > 0) {
    for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
      mesh.point(*v_it) /= maxLength;
  }

  computeCurvature(mesh, curvature);

  up = Vec3f(0, 1, 0);
  pan = Vec3f(0, 0, 0);

  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);
  computeViewCurvature(mesh, actualCamPos, curvature, viewCurvature, viewCurvatureDerivative);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(windowWidth, windowHeight);
  glutCreateWindow(argv[0]);

  glutDisplayFunc(display);
  glutMotionFunc(mouseMoved);
  glutMouseFunc(mouse);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);

  glutMainLoop();

  return 0;
}
