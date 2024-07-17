#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 

  CMU462::Matrix3x3 m[4];
  for (size_t i = 0; i < 4; i++)
	  m[i] = CMU462::Matrix3x3::identity();

  m[0][2] = CMU462::Vector3D(-centerX, -centerY, 1);
  m[1][0][0] = 1 / vspan;
  m[1][1][1] = 1 / vspan;
  m[2][2] = CMU462::Vector3D(1, 1, 1);
  m[3][0][0] = 0.5;
  m[3][1][1] = 0.5;
  set_svg_2_norm(m[3] * m[2] * m[1] * m[0]);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
