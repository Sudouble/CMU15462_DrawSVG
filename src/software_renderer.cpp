#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {
    

  // set top level transformation
  transformation = screen_to_buffer * svg_2_screen;
  
  std::fill(this->sample_buffer.begin(), this->sample_buffer.end(), 0);

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  screen_to_buffer = CMU462::Matrix3x3::identity();
  screen_to_buffer(0, 0) = sample_rate;
  screen_to_buffer(1, 1) = sample_rate;

  this->h = this->sample_rate * this->target_h;
  this->w = this->sample_rate * this->target_w;
  this->sample_buffer.resize(4 * (this->h * this->w));
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  this->h = this->sample_rate * this->target_h;
  this->w = this->sample_rate * this->target_w;
  this->sample_buffer.resize(4 * (this->h * this->w));
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack
  transformation = transformation * element->transform;
  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }
  transformation = transformation * element->transform.inv();

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates
void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= w ) return;
  if ( sy < 0 || sy >= h ) return;

  // fill sample - NOT doing alpha blending!
  for (int ki = (1 - int(sample_rate)) / 2; ki <= int(sample_rate) / 2; ki++) {
      for (int kj = (1 - int(sample_rate)) / 2; kj <= int(sample_rate) / 2; kj++)
      {
          int x = clamp<int>(sx + ki, 0, w - 1);
          int y = clamp<int>(sy + kj, 0, h - 1);
          sample_buffer[4 * (x + y * w) + 0] = (uint8_t)(color.r * 255);
          sample_buffer[4 * (x + y * w) + 1] = (uint8_t)(color.g * 255);
          sample_buffer[4 * (x + y * w) + 2] = (uint8_t)(color.b * 255);
          sample_buffer[4 * (x + y * w) + 3] = (uint8_t)(color.a * 255);
      }
  }
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

    // Task 2: 
    bool isIncreasingByX = true;
    float dx = x1 - x0;
    float dy = y1 - y0;
    if (x0 < x1
        && (((dy / dx) <= 1 && (dy / dx) >= 0)   //octant 1 clock wise
            || ((dy / dx) <= 0 && (dy / dx) >= -1)  //octant 2 clock wise
            )
        )
    {
        isIncreasingByX = true;
    }
    else if (x1 < x0
        && (((dy / dx) <= 1 && (dy / dx) >= 0)   //octant 5 clock wise
            || ((dy / dx) <= 0 && (dy / dx) >= -1)  //octant 6 clock wise
            )
        )
    {
        isIncreasingByX = true;
        swap(x0, x1);
        swap(y0, y1);
    }
    else if (y0 < y1
        && ((dy / dx) > 1       //octant 4 clock wise
            || (dy / dx) < -1   //octant 3 clock wise
            )
        )
    {
        isIncreasingByX = false;
        swap(x0, y0);
        swap(x1, y1);
    }
    else if (y1 < y0
        && ((dy / dx) > 1       //octant 7 clock wise
            || (dy / dx) < -1   //octant 8 clock wise
            )
        )
    {
        isIncreasingByX = false;
        swap(x0, y0);
        swap(x1, y1);

        swap(x0, x1);
        swap(y0, y1);
    }

    dx = x1 - x0;
    dy = y1 - y0;
    float eps = 0;
    int sy = int(floor(y0));

    for (int x = int(floor(x0)); x < int(floor(x1)); ++x)
    {
        if (isIncreasingByX == false)
        {
            swap(x, sy);
        }
        rasterize_point(x, sy, color);
        if (isIncreasingByX == false)
        {
            swap(x, sy);
        }

        if ((dx != 0 && dy/dx >= 0) || dx == 0)
        {
            eps += dy;
            if ((eps * 2) >= dx)
            {
                sy++;
                eps -= dx;
            }
        }
        else
        {
            eps = eps + dy;
            if ((eps*2) <= dx)
            {
                sy = floor(sy - 1);
                eps += dx;
            }
        }
    }
}

bool IsPointInTriangle(float px, float py,
    float x0, float y0,
    float x1, float y1,
    float x2, float y2)

{
    // https://blog.csdn.net/weixin_39872717/article/details/77368234
    // https://blog.csdn.net/weixin_44120025/article/details/123830197
    float AB_x = (x1 - x0); float AB_y = (y1 - y0);
    float BC_x = (x2 - x1); float BC_y = (y2 - y1);
    float CA_x = (x0 - x2); float CA_y = (y0 - y2);

    float AP_x = (px - x0); float AP_y = py - y0;
    float BP_x = px - x1; float BP_y = py - y1;
    float CP_x = px - x2; float CP_y = py - y2;

    bool isCounterClock = ((AB_x * BC_y - BC_x * AB_y) > 0);
    //std::cout << "isCounterClock: " << isCounterClock << std::endl;
    float AB_AP = AB_x * AP_y - AB_y * AP_x;
    float BC_BP = BC_x * BP_y - BC_y * BP_x;
    float CA_CP = CA_x * CP_y - CA_y * CP_x;
    //std::cout << AB_AP << " " << BC_BP << " " << CA_CP << std::endl;
    if ((isCounterClock && AB_AP > 0 && BC_BP > 0 && CA_CP > 0)
        || (!isCounterClock && AB_AP < 0 && BC_BP < 0 && CA_CP < 0)
        || fabs(AB_AP * BC_BP * CA_CP) < 0.0001)
    {
        return true;
    }
    return false;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
    // Task 3: 
    size_t left = floor(min({ x0,x1,x2 }));
    size_t right = ceil(std::max({ x0, x1, x2 }));
    size_t bottom = floor(min({ y0,y1,y2 }));
    size_t top = ceil(std::max({ y0,y1,y2 }));

    for (float x = left; x <= right; x++) {
        for (float y = bottom; y <= top; y++) {
            if (IsPointInTriangle(x, y, x0, y0, x1, y1, x2, y2))
            {
                rasterize_point(x, y, color);
            }
        }
    }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {
  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
    // fill in the nearest pixel
    for (size_t i = 0; i < target_h; i++) {
        for (size_t j = 0; j < target_w; j++) {
            CMU462::Color color(0, 0, 0, 0);
            for (size_t ki = 0; ki < sample_rate; ki++) {
                for (size_t kj = 0; kj < sample_rate; kj++) {
                    size_t idx = (i * sample_rate + ki) * target_w * sample_rate * 4 + (j * sample_rate + kj) * 4;
                    float r = sample_buffer[idx + 0];
                    float g = sample_buffer[idx + 1];
                    float b = sample_buffer[idx + 2];
                    float a = sample_buffer[idx + 3];
                    color += CMU462::Color(r, g, b, a);
                }
            }
            color *= 1.0f / (sample_rate * sample_rate);
            size_t idx = 4 * (i * target_w + j);
            if (color.a != 0) {
                color.r *= 1.0 / color.a * 255.0f;
                color.g *= 1.0 / color.a * 255.0f;
                color.b *= 1.0 / color.a * 255.0f;
            }

            render_target[idx + 0] = color.r;
            render_target[idx + 1] = color.g;
            render_target[idx + 2] = color.b;
            render_target[idx + 3] = color.a;
        }
    }

}


} // namespace CMU462
