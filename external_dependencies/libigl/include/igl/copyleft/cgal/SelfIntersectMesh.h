// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_SELFINTERSECTMESH_H
#define IGL_COPYLEFT_CGAL_SELFINTERSECTMESH_H

#include "CGAL_includes.hpp"
#include "RemeshSelfIntersectionsParam.h"

#include <Eigen/Dense>
#include <list>
#include <map>
#include <vector>

//#define IGL_SELFINTERSECTMESH_DEBUG
#ifndef IGL_FIRST_HIT_EXCEPTION
#define IGL_FIRST_HIT_EXCEPTION 10
#endif

// The easiest way to keep track of everything is to use a class

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Kernel is a CGAL kernel like:
      //     CGAL::Exact_predicates_inexact_constructions_kernel
      // or 
      //     CGAL::Exact_predicates_exact_constructions_kernel
    
      template <
        typename Kernel,
        typename DerivedV,
        typename DerivedF,
        typename DerivedVV,
        typename DerivedFF,
        typename DerivedIF,
        typename DerivedJ,
        typename DerivedIM>
      class SelfIntersectMesh
      {
        typedef 
          SelfIntersectMesh<
          Kernel,
          DerivedV,
          DerivedF,
          DerivedVV,
          DerivedFF,
          DerivedIF,
          DerivedJ,
          DerivedIM> Self;
        public:
          // 3D Primitives
          typedef CGAL::Point_3<Kernel>    Point_3;
          typedef CGAL::Segment_3<Kernel>  Segment_3; 
          typedef CGAL::Triangle_3<Kernel> Triangle_3; 
          typedef CGAL::Plane_3<Kernel>    Plane_3;
          typedef CGAL::Tetrahedron_3<Kernel> Tetrahedron_3; 
          //typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3; 
          //typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3; 
          // 2D Primitives
          typedef CGAL::Point_2<Kernel>    Point_2;
          typedef CGAL::Segment_2<Kernel>  Segment_2; 
          typedef CGAL::Triangle_2<Kernel> Triangle_2; 
          // 2D Constrained Delaunay Triangulation types
          typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
          typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTFB_2;
          typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
          typedef CGAL::Exact_intersections_tag Itag;
          typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> 
            CDT_2;
          typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;
          // Axis-align boxes for all-pairs self-intersection detection
          typedef std::vector<Triangle_3> Triangles;
          typedef typename Triangles::iterator TrianglesIterator;
          typedef typename Triangles::const_iterator TrianglesConstIterator;
          typedef 
            CGAL::Box_intersection_d::Box_with_handle_d<double,3,TrianglesIterator> 
            Box;
    
          // Input mesh
          const Eigen::PlainObjectBase<DerivedV> & V;
          const Eigen::PlainObjectBase<DerivedF> & F;
          // Number of self-intersecting triangle pairs
          typedef typename DerivedF::Index Index;
          Index count;
          typedef std::vector<CGAL::Object> ObjectList;
          // Using a vector here makes this **not** output sensitive
          Triangles T;
          typedef std::vector<Index> IndexList;
          IndexList lIF;
          // #F-long list of faces with intersections mapping to the order in
          // which they were first found
          std::map<Index,std::pair<Index,ObjectList> > offending;
          // Make a short name for the edge map's key
          typedef std::pair<Index,Index> EMK;
          // Make a short name for the type stored at each edge, the edge map's
          // value
          typedef std::vector<Index> EMV;
          // Make a short name for the edge map
          typedef std::map<EMK,EMV> EdgeMap;
          // Maps edges of offending faces to all incident offending faces
          EdgeMap edge2faces;
        public:
          RemeshSelfIntersectionsParam params;
        public:
          // Constructs (VV,FF) a new mesh with self-intersections of (V,F)
          // subdivided
          //
          // See also: remesh_self_intersections.h
          inline SelfIntersectMesh(
              const Eigen::PlainObjectBase<DerivedV> & V,
              const Eigen::PlainObjectBase<DerivedF> & F,
              const RemeshSelfIntersectionsParam & params,
              Eigen::PlainObjectBase<DerivedVV> & VV,
              Eigen::PlainObjectBase<DerivedFF> & FF,
              Eigen::PlainObjectBase<DerivedIF> & IF,
              Eigen::PlainObjectBase<DerivedJ> & J,
              Eigen::PlainObjectBase<DerivedIM> & IM);
        private:
          // Helper function to mark a face as offensive
          //
          // Inputs:
          //   f  index of face in F
          inline void mark_offensive(const Index f);
          // Helper function to count intersections between faces
          //
          // Input:
          //   fa  index of face A in F
          //   fb  index of face B in F
          inline void count_intersection( const Index fa, const Index fb);
          // Helper function for box_intersect. Intersect two triangles A and B,
          // append the intersection object (point,segment,triangle) to a running
          // list for A and B
          //
          // Inputs:
          //   A  triangle in 3D
          //   B  triangle in 3D
          //   fa  index of A in F (and key into offending)
          //   fb  index of A in F (and key into offending)
          // Returns true only if A intersects B
          //
          inline bool intersect(
              const Triangle_3 & A, 
              const Triangle_3 & B, 
              const Index fa,
              const Index fb);
          // Helper function for box_intersect. In the case where A and B have
          // already been identified to share a vertex, then we only want to add
          // possible segment intersections. Assumes truly duplicate triangles are
          // not given as input
          //
          // Inputs:
          //   A  triangle in 3D
          //   B  triangle in 3D
          //   fa  index of A in F (and key into offending)
          //   fb  index of B in F (and key into offending)
          //   va  index of shared vertex in A (and key into offending)
          //   vb  index of shared vertex in B (and key into offending)
          //// Returns object of intersection (should be Segment or point)
          //   Returns true if intersection (besides shared point)
          //
          inline bool single_shared_vertex(
              const Triangle_3 & A,
              const Triangle_3 & B,
              const Index fa,
              const Index fb,
              const Index va,
              const Index vb);
          // Helper handling one direction
          inline bool single_shared_vertex(
              const Triangle_3 & A,
              const Triangle_3 & B,
              const Index fa,
              const Index fb,
              const Index va);
          // Helper function for box_intersect. In the case where A and B have
          // already been identified to share two vertices, then we only want to add
          // a possible coplanar (Triangle) intersection. Assumes truly degenerate
          // facets are not givin as input.
          inline bool double_shared_vertex(
              const Triangle_3 & A,
              const Triangle_3 & B,
              const Index fa,
              const Index fb,
              const std::vector<std::pair<Index,Index> > shared);
    
        public:
          // Callback function called during box self intersections test. Means
          // boxes a and b intersect. This method then checks if the triangles in
          // each box intersect and if so, then processes the intersections
          //
          // Inputs:
          //   a  box containing a triangle
          //   b  box containing a triangle
          inline void box_intersect(const Box& a, const Box& b);
        private:
          // Compute 2D delaunay triangulation of a given 3d triangle and a list of
          // intersection objects (points,segments,triangles). CGAL uses an affine
          // projection rather than an isometric projection, so we're not
          // guaranteed that the 2D delaunay triangulation here will be a delaunay
          // triangulation in 3D.
          //
          // Inputs:
          //   A  triangle in 3D
          //   A_objects_3  updated list of intersection objects for A
          // Outputs:
          //   cdt  Contrained delaunay triangulation in projected 2D plane
        public:
          // Getters:
          //const IndexList& get_lIF() const{ return lIF;}
          static inline void box_intersect_static(
            SelfIntersectMesh * SIM, 
            const Box &a, 
            const Box &b);
      };
    }
  }
}

// Implementation

#include "mesh_to_cgal_triangle_list.h"
#include "remesh_intersections.h"
#include "remesh_intersections.h"

#include "../../REDRUM.h"
#include "../../get_seconds.h"
#include "../../C_STR.h"


#include <functional>
#include <algorithm>
#include <exception>
#include <cassert>
#include <iostream>

// References:
// http://minregret.googlecode.com/svn/trunk/skyline/src/extern/CGAL-3.3.1/examples/Polyhedron/polyhedron_self_intersection.cpp
// http://www.cgal.org/Manual/3.9/examples/Boolean_set_operations_2/do_intersect.cpp

// Q: Should we be using CGAL::Polyhedron_3?
// A: No! Input is just a list of unoriented triangles. Polyhedron_3 requires
// a 2-manifold.
// A: But! It seems we could use CGAL::Triangulation_3. Though it won't be easy
// to take advantage of functions like insert_in_facet because we want to
// constrain segments. Hmmm. Actualy Triangulation_3 doesn't look right...

//static void box_intersect(SelfIntersectMesh * SIM,const Box & A, const Box & B)
//{
//  return SIM->box_intersect(A,B);
//}


// CGAL's box_self_intersection_d uses C-style function callbacks without
// userdata. This is a leapfrog method for calling a member function. It should
// be bound as if the prototype was:
//   static void box_intersect(const Box &a, const Box &b)
// using boost:
//  boost::function<void(const Box &a,const Box &b)> cb
//    = boost::bind(&::box_intersect, this, _1,_2);
//   
template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline void igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::box_intersect_static(
  Self * SIM, 
  const typename Self::Box &a, 
  const typename Self::Box &b)
{
  SIM->box_intersect(a,b);
}

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::SelfIntersectMesh(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const RemeshSelfIntersectionsParam & params,
  Eigen::PlainObjectBase<DerivedVV> & VV,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedIF> & IF,
  Eigen::PlainObjectBase<DerivedJ> & J,
  Eigen::PlainObjectBase<DerivedIM> & IM):
  V(V),
  F(F),
  count(0),
  T(),
  lIF(),
  offending(),
  edge2faces(),
  params(params)
{
  using namespace std;
  using namespace Eigen;

#ifdef IGL_SELFINTERSECTMESH_DEBUG
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  tictoc();
#endif

  // Compute and process self intersections
  mesh_to_cgal_triangle_list(V,F,T);
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  cout<<"mesh_to_cgal_triangle_list: "<<tictoc()<<endl;
#endif
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Box_intersection_d/Chapter_main.html#Section_63.5 
  // Create the corresponding vector of bounding boxes
  std::vector<Box> boxes;
  boxes.reserve(T.size());
  for ( 
    TrianglesIterator tit = T.begin(); 
    tit != T.end(); 
    ++tit)
  {
    boxes.push_back(Box(tit->bbox(), tit));
  }
  // Leapfrog callback
  std::function<void(const Box &a,const Box &b)> cb = 
    std::bind(&box_intersect_static, this, 
      // Explicitly use std namespace to avoid confusion with boost (who puts
      // _1 etc. in global namespace)
      std::placeholders::_1,
      std::placeholders::_2);
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  cout<<"boxes and bind: "<<tictoc()<<endl;
#endif
  // Run the self intersection algorithm with all defaults
  try{
    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(),cb);
  }catch(int e)
  {
    // Rethrow if not IGL_FIRST_HIT_EXCEPTION
    if(e != IGL_FIRST_HIT_EXCEPTION)
    {
      throw e;
    }
    // Otherwise just fall through
  }
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  cout<<"box_self_intersection_d: "<<tictoc()<<endl;
#endif

  // Convert lIF to Eigen matrix
  assert(lIF.size()%2 == 0);
  IF.resize(lIF.size()/2,2);
  {
    Index i=0;
    for(
      typename IndexList::const_iterator ifit = lIF.begin();
      ifit!=lIF.end();
      )
    {
      IF(i,0) = (*ifit);
      ifit++; 
      IF(i,1) = (*ifit);
      ifit++;
      i++;
    }
  }
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  cout<<"IF: "<<tictoc()<<endl;
#endif

  if(params.detect_only)
  {
    return;
  }

  remesh_intersections(V,F,T,offending,edge2faces,VV,FF,J,IM);

  // Q: Does this give the same result as TETGEN?
  // A: For the cow and beast, yes.

  // Q: Is tetgen faster than this CGAL implementation?
  // A: Well, yes. But Tetgen is only solving the detection (predicates)
  // problem. This is also appending the intersection objects (construction).
  // But maybe tetgen is still faster for the detection part. For example, this
  // CGAL implementation on the beast takes 98 seconds but tetgen detection
  // takes 14 seconds


}


template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline void igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::mark_offensive(const Index f)
{
  using namespace std;
  lIF.push_back(f);
  if(offending.count(f) == 0)
  {
    // first time marking, initialize with new id and empty list
    const Index id = offending.size();
    offending[f] = {id,{}};
    for(Index e = 0; e<3;e++)
    {
      // append face to edge's list
      Index i = F(f,(e+1)%3) < F(f,(e+2)%3) ? F(f,(e+1)%3) : F(f,(e+2)%3);
      Index j = F(f,(e+1)%3) < F(f,(e+2)%3) ? F(f,(e+2)%3) : F(f,(e+1)%3);
      edge2faces[EMK(i,j)].push_back(f);
    }
  }
}

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline void igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::count_intersection(
  const Index fa,
  const Index fb)
{
  mark_offensive(fa);
  mark_offensive(fb);
  this->count++;
  // We found the first intersection
  if(params.first_only && this->count >= 1)
  {
    throw IGL_FIRST_HIT_EXCEPTION;
  }
}

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline bool igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::intersect(
  const Triangle_3 & A, 
  const Triangle_3 & B, 
  const Index fa,
  const Index fb)
{
  // Determine whether there is an intersection
  if(!CGAL::do_intersect(A,B))
  {
    return false;
  }
  count_intersection(fa,fb);
  if(!params.detect_only)
  {
    // Construct intersection
    CGAL::Object result = CGAL::intersection(A,B);
    offending[fa].second.push_back(result);
    offending[fb].second.push_back(result);
  }
  return true;
}

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline bool igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::single_shared_vertex(
  const Triangle_3 & A,
  const Triangle_3 & B,
  const Index fa,
  const Index fb,
  const Index va,
  const Index vb)
{
  ////using namespace std;
  //CGAL::Object result = CGAL::intersection(A,B);
  //if(CGAL::object_cast<Segment_3 >(&result))
  //{
  //  // Append to each triangle's running list
  //  F_objects[fa].push_back(result);
  //  F_objects[fb].push_back(result);
  //  count_intersection(fa,fb);
  //}else
  //{
  //  // Then intersection must be at point
  //  // And point must be at shared vertex
  //  assert(CGAL::object_cast<Point_3>(&result));
  //}
  if(single_shared_vertex(A,B,fa,fb,va))
  {
    return true;
  }
  return single_shared_vertex(B,A,fb,fa,vb);
}

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline bool igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::single_shared_vertex(
  const Triangle_3 & A,
  const Triangle_3 & B,
  const Index fa,
  const Index fb,
  const Index va)
{
  // This was not a good idea. It will not handle coplanar triangles well.
  using namespace std;
  Segment_3 sa(
    A.vertex((va+1)%3),
    A.vertex((va+2)%3));

  if(CGAL::do_intersect(sa,B))
  {
    // can't put count_intersection(fa,fb) here since we use intersect below
    // and then it will be counted twice.
    if(params.detect_only)
    {
      count_intersection(fa,fb);
      return true;
    }
    CGAL::Object result = CGAL::intersection(sa,B);
    if(const Point_3 * p = CGAL::object_cast<Point_3 >(&result))
    {
      // Single intersection --> segment from shared point to intersection
      CGAL::Object seg = CGAL::make_object(Segment_3(
        A.vertex(va),
        *p));
      count_intersection(fa,fb);
      offending[fa].second.push_back(seg);
      offending[fb].second.push_back(seg);
      return true;
    }else if(CGAL::object_cast<Segment_3 >(&result))
    {
      //cerr<<REDRUM("Coplanar at: "<<fa<<" & "<<fb<<" (single shared).")<<endl;
      // Must be coplanar
      // WRONG:
      //// Segment intersection --> triangle from shared point to intersection
      //CGAL::Object tri = CGAL::make_object(Triangle_3(
      //  A.vertex(va),
      //  s->vertex(0),
      //  s->vertex(1)));
      //F_objects[fa].push_back(tri);
      //F_objects[fb].push_back(tri);
      //count_intersection(fa,fb);
      // Need to do full test. Intersection could be a general poly.
      bool test = intersect(A,B,fa,fb);
      ((void)test);
      assert(test && "intersect should agree with do_intersect");
      return true;
    }else
    {
      cerr<<REDRUM("Segment ??? triangle neither point nor segment?")<<endl;
      assert(false);
    }
  }

  return false;
}


template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline bool igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::double_shared_vertex(
  const Triangle_3 & A,
  const Triangle_3 & B,
  const Index fa,
  const Index fb,
  const std::vector<std::pair<Index,Index> > shared)
{
  using namespace std;

  // must be co-planar
  if(
    A.supporting_plane() != B.supporting_plane() &&
    A.supporting_plane() != B.supporting_plane().opposite())
  {
    return false;
  }
  // Since A and B are non-degenerate the intersection must be a polygon
  // (triangle). Either
  //   - the vertex of A (B) opposite the shared edge of lies on B (A), or
  //   - an edge of A intersects and edge of B without sharing a vertex



  // Determine if the vertex opposite edge (a0,a1) in triangle A lies in
  // (intersects) triangle B
  const auto & opposite_point_inside = [](
    const Triangle_3 & A, const Index a0, const Index a1, const Triangle_3 & B) 
    -> bool
  {
    // get opposite index
    Index a2 = -1;
    for(int c = 0;c<3;c++)
    {
      if(c != a0 && c != a1)
      {
        a2 = c;
        break;
      }
    }
    assert(a2 != -1);
    bool ret = CGAL::do_intersect(A.vertex(a2),B);
    //cout<<"opposite_point_inside: "<<ret<<endl;
    return ret;
  };

  // Determine if edge opposite vertex va in triangle A intersects edge
  // opposite vertex vb in triangle B.
  const auto & opposite_edges_intersect = [](
    const Triangle_3 & A, const Index va,
    const Triangle_3 & B, const Index vb) -> bool
  {
    Segment_3 sa( A.vertex((va+1)%3), A.vertex((va+2)%3));
    Segment_3 sb( B.vertex((vb+1)%3), B.vertex((vb+2)%3));
    //cout<<sa<<endl;
    //cout<<sb<<endl;
    bool ret = CGAL::do_intersect(sa,sb);
    //cout<<"opposite_edges_intersect: "<<ret<<endl;
    return ret;
  };


  if( 
    !opposite_point_inside(A,shared[0].first,shared[1].first,B) &&
    !opposite_point_inside(B,shared[0].second,shared[1].second,A) &&
    !opposite_edges_intersect(A,shared[0].first,B,shared[1].second) && 
    !opposite_edges_intersect(A,shared[1].first,B,shared[0].second))
  {
    return false;
  }

  // there is an intersection indeed
  count_intersection(fa,fb);
  if(params.detect_only)
  {
    return true;
  }
  // Construct intersection
  try
  {
    // This can fail for Epick but not Epeck
    CGAL::Object result = CGAL::intersection(A,B);
    if(!result.empty())
    {
      if(CGAL::object_cast<Segment_3 >(&result))
      {
        // not coplanar
        assert(false && 
          "Co-planar non-degenerate triangles should intersect over triangle");
        return false;
      } else if(CGAL::object_cast<Point_3 >(&result))
      {
        // this "shouldn't" happen but does for inexact
        assert(false && 
          "Co-planar non-degenerate triangles should intersect over triangle");
        return false;
      } else
      {
        // Triangle object
        offending[fa].second.push_back(result);
        offending[fb].second.push_back(result);
        //cerr<<REDRUM("Coplanar at: "<<fa<<" & "<<fb<<" (double shared).")<<endl;
        return true;
      }
    }else
    {
      // CGAL::intersection is disagreeing with do_intersect
      assert(false && "CGAL::intersection should agree with predicate tests");
      return false;
    }
  }catch(...)
  {
    // This catches some cgal assertion:
    //     CGAL error: assertion violation!
    //     Expression : is_finite(d)
    //     File       : /opt/local/include/CGAL/GMP/Gmpq_type.h
    //     Line       : 132
    //     Explanation: 
    // But only if NDEBUG is not defined, otherwise there's an uncaught
    // "Floating point exception: 8" SIGFPE
    return false;
  }
  // No intersection.
  return false;
}

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedF,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedIF,
  typename DerivedJ,
  typename DerivedIM>
inline void igl::copyleft::cgal::SelfIntersectMesh<
  Kernel,
  DerivedV,
  DerivedF,
  DerivedVV,
  DerivedFF,
  DerivedIF,
  DerivedJ,
  DerivedIM>::box_intersect(
  const Box& a, 
  const Box& b)
{
  using namespace std;
  // Could we write this as a static function of:
  //
  // F.row(fa)
  // F.row(fb)
  // A
  // B

  // index in F and T
  Index fa = a.handle()-T.begin();
  Index fb = b.handle()-T.begin();
  const Triangle_3 & A = *a.handle();
  const Triangle_3 & B = *b.handle();
  // I'm not going to deal with degenerate triangles, though at some point we
  // should
  assert(!a.handle()->is_degenerate());
  assert(!b.handle()->is_degenerate());
  // Number of combinatorially shared vertices
  Index comb_shared_vertices = 0;
  // Number of geometrically shared vertices (*not* including combinatorially
  // shared)
  Index geo_shared_vertices = 0;
  // Keep track of shared vertex indices
  std::vector<std::pair<Index,Index> > shared;
  Index ea,eb;
  for(ea=0;ea<3;ea++)
  {
    for(eb=0;eb<3;eb++)
    {
      if(F(fa,ea) == F(fb,eb))
      {
        comb_shared_vertices++;
        shared.emplace_back(ea,eb);
      }else if(A.vertex(ea) == B.vertex(eb))
      {
        geo_shared_vertices++;
        shared.emplace_back(ea,eb);
      }
    }
  }
  const Index total_shared_vertices = comb_shared_vertices + geo_shared_vertices;
  if(comb_shared_vertices== 3)
  {
    assert(shared.size() == 3);
    //// Combinatorially duplicate face, these should be removed by preprocessing
    //cerr<<REDRUM("Facets "<<fa<<" and "<<fb<<" are combinatorial duplicates")<<endl;
    goto done;
  }
  if(total_shared_vertices== 3)
  {
    assert(shared.size() == 3);
    //// Geometrically duplicate face, these should be removed by preprocessing
    //cerr<<REDRUM("Facets "<<fa<<" and "<<fb<<" are geometrical duplicates")<<endl;
    goto done;
  }
  //// SPECIAL CASES ARE BROKEN FOR COPLANAR TRIANGLES
  //if(total_shared_vertices > 0)
  //{
  //  bool coplanar = 
  //    CGAL::coplanar(A.vertex(0),A.vertex(1),A.vertex(2),B.vertex(0)) &&
  //    CGAL::coplanar(A.vertex(0),A.vertex(1),A.vertex(2),B.vertex(1)) &&
  //    CGAL::coplanar(A.vertex(0),A.vertex(1),A.vertex(2),B.vertex(2));
  //  if(coplanar)
  //  {
  //    cerr<<MAGENTAGIN("Facets "<<fa<<" and "<<fb<<
  //      " are coplanar and share vertices")<<endl;
  //    goto full;
  //  }
  //}

  if(total_shared_vertices == 2)
  {
    assert(shared.size() == 2);
    // Q: What about coplanar?
    //
    // o    o
    // |\  /|
    // | \/ |
    // | /\ |
    // |/  \|
    // o----o
    double_shared_vertex(A,B,fa,fb,shared);

    goto done;
  }
  assert(total_shared_vertices<=1);
  if(total_shared_vertices==1)
  {
//#ifndef NDEBUG
//    CGAL::Object result =
//#endif
    single_shared_vertex(A,B,fa,fb,shared[0].first,shared[0].second);
//#ifndef NDEBUG
//    if(!CGAL::object_cast<Segment_3 >(&result))
//    {
//      const Point_3 * p = CGAL::object_cast<Point_3 >(&result);
//      assert(p);
//      for(int ea=0;ea<3;ea++)
//      {
//        for(int eb=0;eb<3;eb++)
//        {
//          if(F(fa,ea) == F(fb,eb))
//          {
//            assert(*p==A.vertex(ea));
//            assert(*p==B.vertex(eb));
//          }
//        }
//      }
//    }
//#endif
  }else
  {
//full:
    // No geometrically shared vertices, do general intersect
    intersect(*a.handle(),*b.handle(),fa,fb);
  }
done:
  return;
}

#endif

