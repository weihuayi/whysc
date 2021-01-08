#ifndef Cell_type_h
#define Cell_type_h

namespace WHYSC {

namespace Mesh {

// 这里的定义来自VTK
enum Cell_type 
{
  // Linear cells
  EMPTY_CELL       = 0,
  VERTEX           = 1,
  POLY_VERTEX      = 2,
  LINE             = 3,
  POLY_LINE        = 4,
  TRIANGLE         = 5,
  TRIANGLE_STRIP   = 6,
  POLYGON          = 7,
  PIXEL            = 8,
  QUAD             = 9,
  TETRA            = 10,
  VOXEL            = 11,
  HEXAHEDRON       = 12,
  WEDGE            = 13,
  PYRAMID          = 14,
  PENTAGONAL_PRISM = 15,
  HEXAGONAL_PRISM  = 16,
  POLYGON_FACE     = 17,

  // Quadratic, isoparametric cells
  QUADRATIC_EDGE                   = 21,
  QUADRATIC_TRIANGLE               = 22,
  QUADRATIC_QUAD                   = 23,
  QUADRATIC_POLYGON                = 36,
  QUADRATIC_TETRA                  = 24,
  QUADRATIC_HEXAHEDRON             = 25,
  QUADRATIC_WEDGE                  = 26,
  QUADRATIC_PYRAMID                = 27,
  BIQUADRATIC_QUAD                 = 28,
  TRIQUADRATIC_HEXAHEDRON          = 29,
  QUADRATIC_LINEAR_QUAD            = 30,
  QUADRATIC_LINEAR_WEDGE           = 31,
  BIQUADRATIC_QUADRATIC_WEDGE      = 32,
  BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33,
  BIQUADRATIC_TRIANGLE             = 34,

  // Cubic, isoparametric cell
  CUBIC_LINE                       = 35,

  // Special class of cells formed by convex group of points
  CONVEX_POINT_SET = 41,

  // Polyhedron cell (consisting of polygonal faces)
  POLYHEDRON = 42,

  // Higher order cells in parametric form
  PARAMETRIC_CURVE        = 51,
  PARAMETRIC_SURFACE      = 52,
  PARAMETRIC_TRI_SURFACE  = 53,
  PARAMETRIC_QUAD_SURFACE = 54,
  PARAMETRIC_TETRA_REGION = 55,
  PARAMETRIC_HEX_REGION   = 56,

  // Higher order cells
  HIGHER_ORDER_EDGE        = 60,
  HIGHER_ORDER_TRIANGLE    = 61,
  HIGHER_ORDER_QUAD        = 62,
  HIGHER_ORDER_POLYGON     = 63,
  HIGHER_ORDER_TETRAHEDRON = 64,
  HIGHER_ORDER_WEDGE       = 65,
  HIGHER_ORDER_PYRAMID     = 66,
  HIGHER_ORDER_HEXAHEDRON  = 67,

  // Arbitrary order Lagrange elements (formulated separated from generic higher order cells)
  LAGRANGE_CURVE           = 68,
  LAGRANGE_TRIANGLE        = 69,
  LAGRANGE_QUADRILATERAL   = 70,
  LAGRANGE_TETRAHEDRON     = 71,
  LAGRANGE_HEXAHEDRON      = 72,
  LAGRANGE_WEDGE           = 73,
  LAGRANGE_PYRAMID         = 74,
  LAGRANGE_POLYGON         = 75,
  LAGRANGE_POLYHEDRON      = 76,

  NUMBER_OF_CELL_TYPES
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Cell_type_h
