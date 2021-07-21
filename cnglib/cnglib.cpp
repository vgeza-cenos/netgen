/**************************************************************************/
/* File:   nglib.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   7. May. 2000                                                   */
/**************************************************************************/

/*
  
  Interface to the netgen meshing kernel
  
*/
#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>
#include <../visualization/soldata.hpp>

#ifdef OCCGEOMETRY
#include <occgeom.hpp>
namespace netgen
{
   DLL_HEADER extern OCCParameters occparam;
} // namespace netgen
#endif // OCCGEOMETRY

#include <nginterface.h>


namespace netgen {
   extern void MeshFromSpline2D (SplineGeometry2d & geometry,
                                 shared_ptr<Mesh> & mesh, 
                                 MeshingParameters & mp);
   extern MeshingParameters mparam;
   DLL_HEADER extern STLParameters stlparam;
}



#ifdef PARALLEL
#include <mpi.h>

namespace netgen
{
  // int id = 0, ntasks = 1;
  MPI_Comm mesh_comm;
}
#endif


namespace cnglib {
#include "cnglib.h"
}

using namespace netgen;

// constants and types:

namespace cnglib
{
  inline void NOOP_Deleter(void *) { ; }

  
   // initialize, deconstruct Netgen library:
   //copy of Ng_Init
   DLL_HEADER void CNg_Init ()
   {
      mycout = &cout;
      myerr = &cerr;
      // netgen::testout->SetOutStream (new ofstream ("test.out"));
      // testout = new ofstream ("test.out");
   }


   // Clean-up functions before ending usage of cnglib
   //copy of Ng_Exit
   DLL_HEADER void CNg_Exit ()
   {
      ;
   }


   // Create a new netgen mesh object
   //copy of Ng_NewMesh
   DLL_HEADER CNg_Mesh * CNg_NewMesh ()
   {
      Mesh * mesh = new Mesh;  
      //mesh->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
      return (CNg_Mesh*)(void*)mesh;
   }

   // Delete an existing netgen mesh object
   //copy of Ng_DeleteMesh
   DLL_HEADER void CNg_DeleteMesh (CNg_Mesh * mesh)
   {
      if(mesh != NULL)
      {
         // Delete the Mesh structures
         ((Mesh*)mesh)->DeleteMesh();

         // Now delete the Mesh class itself
         delete (Mesh*)mesh;

         // Set the Ng_Mesh pointer to NULL
         mesh = NULL;
      }
   }

   // Copy nodes and elements from origin_mesh into destination_mesh
   DLL_HEADER void CNg_MergeMesh(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh, int index)
   {
       Mesh* origin_mesh = (Mesh*)orig_mesh;
       Mesh* destination_mesh = (Mesh*)dest_mesh;

       //Geometries
       const OCCGeometry& origin_geom = *dynamic_pointer_cast<OCCGeometry>(origin_mesh->GetGeometry());
       const OCCGeometry& destination_geom = *dynamic_pointer_cast<OCCGeometry> (destination_mesh->GetGeometry());
       double eps = 1e-6 * destination_geom.GetBoundingBox().Diam();
			
       //if origin is edge
       int edge_index = 0;
       int solid_index = 0;
       int origin_mesh_dimension = 0;

       if (origin_geom.emap.Extent() == 1 && origin_geom.fmap.Extent() == 0 && origin_geom.somap.Extent() == 0)
       {
           if (index == 0)
               edge_index = destination_geom.emap.FindIndex(origin_geom.emap.FindKey(1));
           else
               edge_index = index;
           origin_mesh_dimension = 1;
       }
       else if (origin_geom.fmap.Extent() == 1 && origin_geom.somap.Extent() == 0) // origin is face
       {
           if (index == 0)
               destination_mesh->AddFaceDescriptor(FaceDescriptor(destination_mesh->GetNFD() + 1, 1, 0, 1));
           else
               destination_mesh->AddFaceDescriptor(FaceDescriptor(index, 1, 0, 1));
           origin_mesh_dimension = 2;
       }
       else if (origin_geom.somap.Extent() == 1)
       {
           if (index == 0)
               solid_index = destination_geom.somap.FindIndex(origin_geom.somap.FindKey(1));
           else
               solid_index = index;
           origin_mesh_dimension = 3;
       }
       else // not supported
       {
           std::string message = "CNg_MergeMesh() does supports only single edge or single face or single solid as origin geometry. Passed Edges: " + std::to_string(origin_geom.emap.Extent())
               + ", Faces: " + std::to_string(origin_geom.fmap.Extent()) + ", Solids: " + std::to_string(origin_geom.somap.Extent()) + ".";
           std::cout << message << std::endl;
           throw(message);
       }

      // map to hold node ids.
      // Example:
	  // node 1 will correspond to node 10, if there are already 9 nodes in destination mesh
      // node 3 will correspond to node 7, if they are equal (say, edge node)
      std::map<PointIndex, PointIndex> node_ids;

      int nNodes = origin_mesh->GetNP();
      for (PointIndex opi: origin_mesh->Points().Range())
      {
          const Point3d& p = (*origin_mesh)[opi];
          bool exists = false;
          for (PointIndex dpi : destination_mesh->Points().Range())
          {
              if (Dist2((*destination_mesh)[dpi], p) < eps * eps)
              {
                  exists = true;
                  node_ids[opi] = dpi;
                  break;
              }
          }
          
          if (!exists)
          {
              destination_mesh->AddPoint(p);
              node_ids[opi] = PointIndex(destination_mesh->GetNP());
          }
      }

      int nEdgeEl = origin_mesh->GetNSeg();
      int nFaceEl = origin_mesh->GetNSE();
      int nVolEl = origin_mesh->GetNE();

      //edge elements
	  if (origin_mesh_dimension == 1)
      {
		  for (int i = 1; i <= nEdgeEl; i++)
		  {
			  Segment seg = origin_mesh->LineSegment(i);

			  seg[0] = node_ids[seg[0]];
			  seg[1] = node_ids[seg[1]];
			  seg.epgeominfo[0].edgenr = edge_index;
			  seg.epgeominfo[1].edgenr = edge_index;
			  seg.edgenr = destination_mesh->GetNSeg() + 1;
			  seg.cd2i = -1;

			  destination_mesh->AddSegment(seg);
		  }
	  }

      //surface elements
      if (origin_mesh_dimension == 2)
      {
          for (int i = 1; i <= nFaceEl; i++)
          {
              //TODO check element type!!!
              Element2d el = origin_mesh->SurfaceElement(i);
              el.PNum(1) = node_ids[el.PNum(1)];
              el.PNum(2) = node_ids[el.PNum(2)];
              el.PNum(3) = node_ids[el.PNum(3)];
              if (el.GetType() == QUAD)
                  el.PNum(4) = node_ids[el.PNum(4)];

              el.SetIndex(destination_mesh->GetNFD());
              destination_mesh->AddSurfaceElement(el);
          }
      }

      //volume elements
      if (origin_mesh_dimension == 3)
      {
          for (int i = 1; i <= nVolEl; i++)
          {
              Element el = origin_mesh->VolumeElement(i);
              for (int n = 1; n <= el.GetNP(); n++)
                el.PNum(n) = node_ids[el.PNum(n)];
              el.SetIndex(solid_index);
              destination_mesh->AddVolumeElement(el);
          }
      }
   }
   
   // Copy segments from origin_mesh into destination_mesh, assign edge index to segments
   DLL_HEADER void CNg_CopySegments(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh, int edge_index)
   {
       Mesh* origin_mesh = (Mesh*)orig_mesh;
       Mesh* destination_mesh = (Mesh*)dest_mesh;

       //Geometries
       const OCCGeometry& origin_geom = *dynamic_pointer_cast<OCCGeometry>(origin_mesh->GetGeometry());
       const OCCGeometry& destination_geom = *dynamic_pointer_cast<OCCGeometry> (destination_mesh->GetGeometry());
       double eps = 1e-6 * destination_geom.GetBoundingBox().Diam();
			
 
      // map to hold node ids.
      // Example:
	  // node 1 will correspond to node 10, if there are already 9 nodes in destination mesh
      // node 3 will correspond to node 7, if they are equal (say, edge node)
      std::map<PointIndex, PointIndex> node_ids;

      int nNodes = origin_mesh->GetNP();
      for (PointIndex opi: origin_mesh->Points().Range())
      {
          const Point3d& p = (*origin_mesh)[opi];
          bool exists = false;
          for (PointIndex dpi : destination_mesh->Points().Range())
          {
              if (Dist2((*destination_mesh)[dpi], p) < eps * eps)
              {
                  exists = true;
                  node_ids[opi] = dpi;
                  break;
              }
          }

          if (!exists)
          {
              destination_mesh->AddPoint(p);
              node_ids[opi] = PointIndex(destination_mesh->GetNP());
          }
      }

      int nEdgeEl = origin_mesh->GetNSeg();

      //copy segments

	  for (int i = 1; i <= nEdgeEl; i++)
	  {
		  Segment seg = origin_mesh->LineSegment(i);

		  seg[0] = node_ids[seg[0]];
		  seg[1] = node_ids[seg[1]];
		  seg.epgeominfo[0].edgenr = edge_index;
		  seg.epgeominfo[1].edgenr = edge_index;
		  seg.edgenr = destination_mesh->GetNSeg() + 1;
		  seg.cd2i = -1;

		  destination_mesh->AddSegment(seg);
	  }
   }


   // Copy surface elements from origin_mesh into destination_mesh, assign edge index to segments
   DLL_HEADER void CNg_CopySurfaceElements(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh, int face_index)
   {
       Mesh* origin_mesh = (Mesh*)orig_mesh;
       Mesh* destination_mesh = (Mesh*)dest_mesh;

       //Geometries
       const OCCGeometry& origin_geom = *dynamic_pointer_cast<OCCGeometry>(origin_mesh->GetGeometry());
       const OCCGeometry& destination_geom = *dynamic_pointer_cast<OCCGeometry> (destination_mesh->GetGeometry());
       double eps = 1e-6 * destination_geom.GetBoundingBox().Diam();


       // map to hold node ids.
       // Example:
       // node 1 will correspond to node 10, if there are already 9 nodes in destination mesh
       // node 3 will correspond to node 7, if they are equal (say, edge node)
       std::map<PointIndex, PointIndex> node_ids;

       int nNodes = origin_mesh->GetNP();
       for (PointIndex opi : origin_mesh->Points().Range())
       {
           const Point3d& p = (*origin_mesh)[opi];
           bool exists = false;
           for (PointIndex dpi : destination_mesh->Points().Range())
           {
               if (Dist2((*destination_mesh)[dpi], p) < eps * eps)
               {
                   exists = true;
                   node_ids[opi] = dpi;
                   break;
               }
           }
       }

       int nFaceEl = origin_mesh->GetNSE();

       // copy surface elements
        for (int i = 1; i <= nFaceEl; i++)
        {
            //TODO check element type!!!
            Element2d el = origin_mesh->SurfaceElement(i);
            el.PNum(1) = node_ids[el.PNum(1)];
            el.PNum(2) = node_ids[el.PNum(2)];
            el.PNum(3) = node_ids[el.PNum(3)];
            if (el.GetType() == QUAD)
                el.PNum(4) = node_ids[el.PNum(4)];
            el.SetIndex(face_index);
            destination_mesh->AddSurfaceElement(el);
        }
   }

   // Prepare surface meshing
   DLL_HEADER void CNg_PrepareSurfaceMeshing(CNg_Mesh* mesh)
   {
       Mesh* occ_mesh = (Mesh*)mesh;

       const OCCGeometry& occgeom = *dynamic_pointer_cast<OCCGeometry> (occ_mesh->GetGeometry());
	   
       TopoDS_Face face = TopoDS::Face(occgeom.fmap.FindKey(1));
       int uvc = 0;
       for (int e = 1; e <= occgeom.emap.Extent(); e++)
       {
           TopoDS_Edge edge = TopoDS::Edge(occgeom.emap.FindKey(e));
           Handle(Geom2d_Curve) cof;
           double s0, s1;
           cof = BRep_Tool::CurveOnSurface(edge, face, s0, s1);

           for (SegmentIndex si = 0; si < occ_mesh->GetNSeg(); si++)
           {
               if ((*occ_mesh)[si].epgeominfo[1].edgenr == e)
               {
                   gp_Pnt2d p2d;
                   uvc++;
                   p2d = cof->Value((*occ_mesh)[si].epgeominfo[0].dist) ;

                   (*occ_mesh)[si].epgeominfo[0].u = p2d.X();

                   (*occ_mesh)[si].epgeominfo[0].v = p2d.Y();

                   p2d = cof->Value((*occ_mesh)[si].epgeominfo[1].dist);
                   (*occ_mesh)[si].epgeominfo[1].u = p2d.X();
                   (*occ_mesh)[si].epgeominfo[1].v = p2d.Y();
               }

           }

       }
   }

   // Reverse all segments in mesh
   DLL_HEADER void CNg_ReverseSegments(CNg_Mesh* mesh)
   {
       Mesh* occ_mesh = (Mesh*)mesh;

       for (SegmentIndex si = 0; si < occ_mesh->GetNSeg(); si++)
       {
           swap((*occ_mesh)[si][0], (*occ_mesh)[si][1]);
           swap((*occ_mesh)[si].epgeominfo[0].dist, (*occ_mesh)[si].epgeominfo[1].dist);
           swap((*occ_mesh)[si].epgeominfo[0].edgenr, (*occ_mesh)[si].epgeominfo[1].edgenr);
           swap((*occ_mesh)[si].epgeominfo[0].u, (*occ_mesh)[si].epgeominfo[1].u);
           swap((*occ_mesh)[si].epgeominfo[0].v, (*occ_mesh)[si].epgeominfo[1].v);
       }
   }


   // Reverse all faces in mesh
   DLL_HEADER void CNg_ReverseFaces(CNg_Mesh* mesh)
   {
       Mesh* occ_mesh = (Mesh*)mesh;

       for (SurfaceElementIndex si = 0; si < occ_mesh->GetNSE(); si++)
       {
           (*occ_mesh)[si].Invert();
       }
   }

   


   // Save a netgen mesh in the native VOL format 
   //copy of Ng_SaveMesh
   DLL_HEADER void CNg_SaveMesh(CNg_Mesh * mesh, const char* filename)
   {
      ((Mesh*)mesh)->Save(filename);
   }

   DLL_HEADER CNg_LocalH CNg_GetLocalH(CNg_Mesh * mesh)
   {
       ((Mesh*)mesh)->CalcLocalHFromPointDistances(0.4);
       return (CNg_LocalH)((Mesh*)mesh)->GetLocalH().get();
   }

   DLL_HEADER void CNg_CopyLocalH(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh)
   {
       Mesh* origin_mesh = (Mesh*)orig_mesh;
       Mesh* destination_mesh = (Mesh*)dest_mesh;
       destination_mesh->SetLocalH(origin_mesh->GetLocalH());
   }

   DLL_HEADER double CNg_GetMaxH(CNg_Mesh* occ_mesh)
   {
       Mesh* mesh = (Mesh*)occ_mesh;
       mesh->CalcLocalHFromPointDistances(0.4);

       double max_h = 0;
       for (PointIndex opi : mesh->Points().Range())
       {
           const Point3d& p = (*mesh)[opi];
           double point_h = mesh->GetH(p);
           if (point_h > max_h)
               max_h = point_h;
       }
       return max_h;
   }
   

   // Manually add a point to an existing mesh object
   //copy of Ng_AddPoint
   DLL_HEADER void CNg_AddPoint (CNg_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;
      m->AddPoint (Point3d (x[0], x[1], x[2]));
   }


		  
	// Manually add a segment element of a given type to an existing mesh object
   DLL_HEADER void CNg_AddSegmentElement(Ng_Mesh * mesh, int pi1, int pi2, int edgeIndex)
   {
     Mesh * m = (Mesh*)mesh;
	 const Point3d & p1 = m->Point(pi1);
	 const Point3d & p2 = m->Point(pi2); 
	 
     Segment seg;
     seg[0] = pi1;
     seg[1] = pi2;
	    
     seg.epgeominfo[0].edgenr = edgeIndex;
     seg.epgeominfo[1].edgenr = edgeIndex;
	 
 
	 seg.epgeominfo[0].u = p1.X();
	 seg.epgeominfo[0].v = p1.Y();
     seg.epgeominfo[1].u = p2.X();
     seg.epgeominfo[1].v = p2.X();
								
     seg.si = 1;
     seg.edgenr = m->GetNSeg() + 1;
 
      m->AddSegment (seg);
   }
   
   // Manually add a surface element of a given type to an existing mesh object
   DLL_HEADER void CNg_AddSurfaceElementUV(CNg_Mesh* mesh, CNg_Surface_Element_Type et,
       int* pi, int surfIndx, double* uv1, double* uv2, double* uv3)
   {
       Mesh* m = (Mesh*)mesh;
       Element2d el;
       if (et == CNG_TRIG)
           el = Element2d(3);
       else if (et == CNG_QUAD)
           el = Element2d(4);

       el.SetIndex(surfIndx);
       el.PNum(1) = pi[0];
       el.PNum(2) = pi[1];
       el.PNum(3) = pi[2];

       //add fourth node if it is quad
       if (et == CNG_QUAD)
           el.PNum(4) = pi[3];

       el.GeomInfoPi(1).u = uv1[0];
       el.GeomInfoPi(1).v = uv1[1];
       el.GeomInfoPi(2).u = uv2[0];
       el.GeomInfoPi(2).v = uv2[1];
       el.GeomInfoPi(3).u = uv3[0];
       el.GeomInfoPi(3).v = uv3[1];
       m->AddSurfaceElement(el);
   }

   // Manually add a surface element of a given type to an existing mesh object with
   // added surface index
   DLL_HEADER void CNg_AddSurfaceElement (CNg_Mesh * mesh, CNg_Surface_Element_Type et,
                                         int * pi, int surfIndx)
   {
       Mesh* m = (Mesh*)mesh;
       Element2d el(3);
       if (et == CNG_TRIG)
           el = Element2d(3);
       else if (et == CNG_QUAD)
           el = Element2d(4);

       el.SetIndex(surfIndx);
       el.PNum(1) = pi[0];
       el.PNum(2) = pi[1];
       el.PNum(3) = pi[2];

       //add fourth node if it is quad
       if (et == CNG_QUAD)
           el.PNum(4) = pi[3];

       m->AddSurfaceElement(el);
   }
   
   // Manually add a volume element of a given type to an existing mesh object
   // added volume index
   // added prism check
   DLL_HEADER void CNg_AddVolumeElement (CNg_Mesh * mesh, CNg_Volume_Element_Type et,
                                        int * pi, int volIndx)
   {
      Mesh * m = (Mesh*)mesh;
	  int nodeCount = 4;
	  if (et == CNg_Volume_Element_Type::CNG_PRISM)
	  {
		  nodeCount = 6;
	  }

	  Element el (nodeCount);
      el.SetIndex (volIndx);
      el.PNum(1) = pi[0];
      el.PNum(2) = pi[1];
      el.PNum(3) = pi[2];
      el.PNum(4) = pi[3];
	  if (et == CNg_Volume_Element_Type::CNG_PRISM)
	  {
		el.PNum(5) = pi[4];
		el.PNum(6) = pi[5];
	  }
		  
      m->AddVolumeElement (el);
   }
   
  
   // Obtain the number of points in the mesh
   // copy of Ng_GetNP
   DLL_HEADER int CNg_GetNP (CNg_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNP();
   }
   
   
   DLL_HEADER int CNg_GetNSeg(CNg_Mesh * mesh)
   {
      return ((Mesh*)mesh) ->GetNSeg();
   }




   // Obtain the number of surface elements in the mesh
   // copy of Ng_GetNSE
   DLL_HEADER int CNg_GetNSE (CNg_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNSE();
   }




   // Obtain the number of volume elements in the mesh
   // copy of Ng_GetNE
   DLL_HEADER int CNg_GetNE (CNg_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNE();
   }




   //  Return point coordinates of a given point index in the mesh
   // copy of Ng_GetPoint
   DLL_HEADER void CNg_GetPoint (CNg_Mesh * mesh, int num, double * x)
   {
      const Point3d & p = ((Mesh*)mesh)->Point(num);
      x[0] = p.X();
      x[1] = p.Y();
      x[2] = p.Z();
   }


   DLL_HEADER int CNg_GetSegment(Ng_Mesh * mesh, int num, int * pi, int * matnum)
   {
      const Segment & seg = ((Mesh*)mesh)->LineSegment(num);
      pi[0] = seg[0];
      pi[1] = seg[1];

      if (matnum)
         *matnum = seg.edgenr;
      if (seg[2] == 0)
      {
          return 1;
      }
      else
      {
          pi[2] = seg[2];
          return 8;
      }
   }
   
   	

   // Return the surface element at a given index "pi"
   // copy of Ng_GetSurfaceElement
   DLL_HEADER CNg_Surface_Element_Type 
      CNg_GetSurfaceElement (CNg_Mesh * mesh, int num, int * pi)
   {
      const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
      CNg_Surface_Element_Type et;
      switch (el.GetNP())
      {
      case 3: et = CNG_TRIG; break;
      case 4: et = CNG_QUAD; break;
      case 6: et = CNG_TRIG6; break;
      case 8: et = CNG_QUAD8; break;
      default:
         et = CNG_TRIG; break; // for the compiler
      }
      return et;
   }




   // Return the volume element at a given index "pi"
   // copy of Ng_GetVolumeElement
   DLL_HEADER CNg_Volume_Element_Type
      CNg_GetVolumeElement (CNg_Mesh * mesh, int num, int * pi)
   {
      const Element & el = ((Mesh*)mesh)->VolumeElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
      CNg_Volume_Element_Type et;
      switch (el.GetNP())
      {
      case 4: et = CNG_TET; break;
      case 5: et = CNG_PYRAMID; break;
      case 6: et = CNG_PRISM; break;
      case 10: et = CNG_TET10; break;
      default:
         et = CNG_TET; break; // for the compiler
      }
      return et;
   }




   // Set a global limit on the maximum mesh size allowed
   // copy of Ng_RestrictMeshSizeGlobal
   DLL_HEADER void CNg_RestrictMeshSizeGlobal (CNg_Mesh * mesh, double h)
   {
      ((Mesh*)mesh) -> SetGlobalH (h);
   }




   // Set a local limit on the maximum mesh size allowed around the given point
   // copy of Ng_RestrictMeshSizePoint
   DLL_HEADER void CNg_RestrictMeshSizePoint (CNg_Mesh * mesh, double * p, double h)
   {
      ((Mesh*)mesh) -> RestrictLocalH (Point3d (p[0], p[1], p[2]), h);
   }


   // Set a local limit on the maximum mesh size allowed within a given box region
   // copy of Ng_RestrictMeshSizeBox
   DLL_HEADER void CNg_RestrictMeshSizeBox (CNg_Mesh * mesh, double * pmin, double * pmax, double h)
   {
      for (double x = pmin[0]; x < pmax[0]; x += h)
         for (double y = pmin[1]; y < pmax[1]; y += h)
            for (double z = pmin[2]; z < pmax[2]; z += h)
               ((Mesh*)mesh) -> RestrictLocalH (Point3d (x, y, z), h);
   }

   // Set a local limit on the maximum mesh size by copyin max size from origin mesh
   DLL_HEADER void CNg_RestrictMeshSizeMesh(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh)
   {
       Mesh* origin_mesh = (Mesh*)orig_mesh;
       Mesh* destination_mesh = (Mesh*)dest_mesh;

       int nNodes = origin_mesh->GetNP();
       for (PointIndex opi : origin_mesh->Points().Range())
       {
           const Point3d& p = (*origin_mesh)[opi];
           destination_mesh->RestrictLocalH(p, origin_mesh->GetH(p));
       }
   }

   // Set a local limit on the maximum mesh size by copyin max size from origin mesh
   DLL_HEADER void CNg_RestrictMeshSizeLocalH(CNg_Mesh* occ_mesh, CNg_LocalH* localh)
   {
       Mesh* mesh = (Mesh*)occ_mesh;
       for (PointIndex opi : mesh->Points().Range())
       {
           const Point3d& p = (*mesh)[opi];
           mesh->RestrictLocalH(p, ((LocalH*)localh)->GetH(p));
       }
   }



   // Generates volume mesh from an existing surface mesh
   //copy of Ng_GenerateVolumeMesh
   DLL_HEADER CNg_Result CNg_GenerateVolumeMesh (CNg_Mesh * mesh, CNg_Meshing_Parameters * mp)
   {
      Mesh * m = (Mesh*)mesh;

      // Philippose - 30/08/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object 
      //MeshingParameters mparam;
      mp->Transfer_Parameters();

      m->CalcLocalH(mparam.grading);

      try 
      {
          MeshVolume(mparam, *m);
          RemoveIllegalElements(*m);
          OptimizeVolume(mparam, *m);
      }
      catch (NgException e)
      {
          std::cout << "Netgen Exception: " << e.what() << std::endl;
          return CNG_ERROR;
      }

      return CNG_OK;
   }

   DLL_HEADER CNg_Result CNg_GenerateBoundaryLayer(CNg_Mesh* mesh,
       int* surfid_arr, int surfid_count,
       double* heights_arr, int heights_count)
   {
       Mesh* m = (Mesh*)mesh;

       BoundaryLayerParameters blp;
       for (int i = 0; i < surfid_count; i++) { blp.surfid.Append(surfid_arr[i]); }
       for (int i = 0; i < heights_count; i++) { blp.heights.Append(heights_arr[i]); }
       blp.new_mat = "BoundaryLayer";
       blp.domains.SetSize(2);
       blp.domains.Clear();
       blp.domains.SetBit(1);
       blp.outside = false;
       blp.grow_edges = true;

       try
       {
           GenerateBoundaryLayer(*m, blp);
        //   RemoveIllegalElements(*m);
           OptimizeVolume(mparam, *m);
       }
       catch (NgException e)
       {
           std::cout << "Netgen Exception: " << e.what() << std::endl;
           return CNG_ERROR;
       }

       return CNG_OK;
   }

   DLL_HEADER CNg_Result CNg_GenerateBoundaryLayer2(CNg_Mesh* mesh,
       int dom_nr, double* heights_arr, int heights_count)
   {
       Mesh* m = (Mesh*)mesh;

       Array<double> thicknesses;
       for (int i = 0; i < heights_count; i++)
       {
           thicknesses.Append(heights_arr[i]);
       }

       try
       {
           GenerateBoundaryLayer2(*m, dom_nr, thicknesses, false);
       }
       catch (NgException e)
       {
           std::cout << "Netgen Exception: " << e.what() << std::endl;
           return CNG_ERROR;
       }

       return CNG_OK;
   }


#ifdef OCCGEOMETRY
   // --------------------- OCC Geometry / Meshing Utility Functions -------------------
   // Create new OCC Geometry Object
   // copy of Ng_OCC_NewGeometry
   DLL_HEADER CNg_OCC_Geometry * CNg_OCC_NewGeometry ()
   {
      return (CNg_OCC_Geometry*)(void*)new OCCGeometry;
   } 



   // Delete the OCC Geometry Object
   // copy of Ng_OCC_DeleteGeometry
   DLL_HEADER CNg_Result CNg_OCC_DeleteGeometry(CNg_OCC_Geometry * geom)
   {
      if (geom != NULL)
      {
         delete (OCCGeometry*)geom;
         geom = NULL;
         return CNG_OK;
      }
      
      return CNG_ERROR;
   }


   // Locally limit the size of the mesh to be generated at various points 
   // based on the topology of the geometry
   // copy of Ng_OCC_SetLocalMeshSize
   DLL_HEADER CNg_Result CNg_SetLocalMeshSize (CNg_OCC_Geometry * geom,
                                                 CNg_Mesh * mesh,
                                                 CNg_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<OCCGeometry>(occgeom, &NOOP_Deleter) );

      me->geomtype = Mesh::GEOM_OCC;

      mp->Transfer_Parameters();
      
      if(mp->closeedgeenable)
        mparam.closeedgefac = mp->closeedgefact;

      // Delete the mesh structures in order to start with a clean 
      // slate
      me->DeleteMesh();

      OCCSetLocalMeshSize(*occgeom, *me, mparam, occparam);

      return(CNG_OK);
   }


   // Set the geometry / topology of mesh
   DLL_HEADER void CNg_SetGeometry(CNg_OCC_Geometry* geom, CNg_Mesh* mesh)
   {
       OCCGeometry* occgeom = (OCCGeometry*)geom;
       Mesh* me = (Mesh*)mesh;
       me->SetGeometry(shared_ptr<OCCGeometry>(occgeom, &NOOP_Deleter));
   }

   // Mesh single edgeId
   DLL_HEADER CNg_Result CNg_DivideEdges(CNg_OCC_Geometry * geom,
                                                 CNg_Mesh * mesh,
                                                 CNg_Meshing_Parameters * meshparams)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<OCCGeometry>(occgeom, &NOOP_Deleter) );
      meshparams->Transfer_Parameters();
      double eps = 1e-6 * occgeom->GetBoundingBox().Diam();

	// add vertices to mesh as points
	  for (int i = 1; i <= occgeom->vmap.Extent();i++)
      {
          gp_Pnt pnt = BRep_Tool::Pnt (TopoDS::Vertex(occgeom->vmap(i)));
          MeshPoint mpnt( Point<3>(pnt.X(), pnt.Y(), pnt.Z()) );

          bool exists = 0;
          for (PointIndex pi : me->Points().Range())
             if (Dist2 ((*me)[pi], Point<3>(mpnt)) < eps*eps)
             {
                 exists = true;
                 break;
             }

          if (!exists)
              me->AddPoint (mpnt);
      }
	  
	  
      PointIndex first_ep = me->Points().Range().end();
      auto vertexrange = me->Points().Range();
	

      for (int edge_i = 1; edge_i <= occgeom->emap.Extent(); edge_i++)
      {
		  TopoDS_Edge edge = TopoDS::Edge(occgeom->emap(edge_i));
		  if (BRep_Tool::Degenerated(edge))
		  {
			//(*testout) << "ignoring degenerated edge" << endl;
			continue;
		  }

		  int geomedgenr = occgeom->emap.FindIndex(edge);
		  if(geomedgenr < 1) continue;

		  if (occgeom->vmap.FindIndex(TopExp::FirstVertex (edge)) ==
              occgeom->vmap.FindIndex(TopExp::LastVertex (edge)))
		  {
			  GProp_GProps system;
			  BRepGProp::LinearProperties(edge, system);

			  if (system.Mass() < eps)
			  {
				cout << "ignoring edge " << occgeom->emap.FindIndex (edge)
					 << ". closed edge with length < " << eps << endl;
				continue;
			  }
		  }


		  NgArray <MeshPoint> mp;
		  NgArray <double> params;

          DivideEdge(edge, mp, params, *me, mparam);


		  NgArray<PointIndex> pnums(mp.Size()+2);
          gp_Pnt p1 = BRep_Tool::Pnt(TopExp::FirstVertex(edge));
		  Point<3> fp = Point<3>(p1.X(), p1.Y(), p1.Z());

          gp_Pnt p2 = BRep_Tool::Pnt(TopExp::LastVertex(edge));
		  Point<3> lp = Point<3>(p2.X(), p2.Y(), p2.Z());

		  pnums[0] = PointIndex::INVALID;
		  pnums.Last() = PointIndex::INVALID;
		  for (PointIndex pi : vertexrange)
		  {
			  if (Dist2 ((*me)[pi], fp) < eps*eps) pnums[0] = pi;
			  if (Dist2 ((*me)[pi], lp) < eps*eps) pnums.Last() = pi;
		  }

		  for (size_t i = 1; i <= mp.Size(); i++)
		  {
			  bool exists = 0;
              for (PointIndex pi : vertexrange)
			    if (((*me)[pi] -Point<3>(mp[i-1])).Length() < eps)
			    {
				    exists = true;
				    pnums[i] = pi;
				    break;
				}
			  		
			  if (!exists)
			      pnums[i] = me->AddPoint (mp[i-1]);
		  }
		
		  (*testout) << "NP = " << me->GetNP() << endl;

		  for (size_t i = 1; i <= mp.Size()+1; i++)
		  {
			  Segment seg;

			  seg[0] = pnums[i-1];
			  seg[1] = pnums[i];
			  seg.edgenr = me->GetNSeg() + 1;
			  seg.si = 1;
			  seg.cd2i = -1;
			  seg.epgeominfo[0].dist = params[i-1];
			  seg.epgeominfo[1].dist = params[i];
			  seg.epgeominfo[0].edgenr = geomedgenr;
			  seg.epgeominfo[1].edgenr = geomedgenr;

			  if (edge.Orientation() == TopAbs_REVERSED)
			  {
				  swap (seg[0], seg[1]);
				  swap (seg.epgeominfo[0].dist, seg.epgeominfo[1].dist);
				  swap (seg.epgeominfo[0].u, seg.epgeominfo[1].u);
				  swap (seg.epgeominfo[0].v, seg.epgeominfo[1].v);
			  }

			  me->AddSegment (seg);
		  }
		}

	  
      if((me->GetNP()) && (me->GetNSeg()))
      {
         return CNG_OK;
      }
      else
      {
         return CNG_ERROR;
      }
   }
   
   // Mesh the edges and add Face descriptors to prepare for surface meshing
   // copy of Ng_OCC_GenerateEdgeMesh
   DLL_HEADER CNg_Result CNg_GenerateEdgeMesh (CNg_OCC_Geometry * geom,
                                                 CNg_Mesh * mesh,
                                                 CNg_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<OCCGeometry>(occgeom, &NOOP_Deleter) );

      mp->Transfer_Parameters();

      OCCFindEdges(*occgeom, *me, mparam);

      if((me->GetNP()) && (me->GetNFD()))
      {
         return CNG_OK;
      }
      else
      {
         return CNG_ERROR;
      }
   }


   // Mesh the edges and add Face descriptors to prepare for surface meshing
   // copy of Ng_OCC_GenerateSurfaceMesh
   DLL_HEADER CNg_Result CNg_GenerateSurfaceMesh (CNg_OCC_Geometry * geom,
                                                    CNg_Mesh * mesh,
                                                    CNg_Meshing_Parameters * mp)
   {
      int numpoints = 0;
      (*mycout) << "Starting GenerateSurfaceMesh" << endl << flush;

      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      if (!me->GetNFD())
          me->AddFaceDescriptor(FaceDescriptor(1, 1, 0, 0));

      me->SetGeometry( shared_ptr<OCCGeometry>(occgeom, &NOOP_Deleter) );
      (*mycout) << "Geometry is set" << endl << flush;

      // Set the internal meshing parameters structure from the nglib meshing 
      // parameters structure
      mp->Transfer_Parameters();
      (*mycout) << "Parameters are transferred" << endl << flush;

      // Only go into surface meshing if the face descriptors have already been added
      if(!me->GetNFD())
         return CNG_ERROR;

      numpoints = me->GetNP();

      // Initially set up only for surface meshing without any optimisation
      int perfstepsend = MESHCONST_MESHSURFACE;

      // Check and if required, enable surface mesh optimisation step
      if(mp->optsurfmeshenable)
      {
         perfstepsend = MESHCONST_OPTSURFACE;
      }
      (*mycout) << "Start mesh surface" << endl << flush;
      try 
      {
          OCCMeshSurface(*occgeom, *me, mparam);

          OCCOptimizeSurface(*occgeom, *me, mparam);

      }
      catch (NgException e)
      {
          std::cout << "Netgen Exception: " << e.what() << std::endl;
          return CNG_ERROR;
      }
      (*mycout) << "Done meshing surface" << endl << flush;
      (*mycout) << "Optimizing surface" << endl << flush;

      (*mycout) << "Done optimizingsurface" << endl << flush;

      me->CalcSurfacesOfNode();
      (*mycout) << "Done CalcSurfacesOfNode" << endl << flush;

      (*mycout) << me->GetNP() << endl << flush;
      (*mycout) << numpoints << endl << flush;
      (*mycout) << me->GetNSE() << endl << flush;

      if (mp->second_order)
          me->GetGeometry()->GetRefinement().MakeSecondOrder(*me);

      if(me->GetNP() < numpoints)
         return CNG_ERROR;

      if(me->GetNSE() <= 0)
         return CNG_ERROR;

      return CNG_OK;
   }



   // Extract the solid map from the OCC geometry
   // The solid map basically gives an index to each solid in the geometry, 
   // which can be used to access a specific solid
    DLL_HEADER CNg_Result CNg_GetSolidMap(CNg_OCC_Geometry* geom,
        CNg_OCC_TopTools_IndexedMapOfShape* SolidMap)
    {
        OCCGeometry* occgeom = (OCCGeometry*)geom;
        TopTools_IndexedMapOfShape* occsomap = (TopTools_IndexedMapOfShape*)SolidMap;

        // Copy the solid map from the geometry to the given variable
        occsomap->Assign(occgeom->somap);

        if (occsomap->Extent())
        {
            return CNG_OK;
        }
        else
        {
            return CNG_ERROR;
        }
    }
	
   // Extract the face map from the OCC geometry
   // The face map basically gives an index to each face in the geometry, 
   // which can be used to access a specific face
   // copy of Ng_OCC_GetFMap
   DLL_HEADER CNg_Result CNg_GetFaceMap(CNg_OCC_Geometry * geom, 
                                       CNg_OCC_TopTools_IndexedMapOfShape * FaceMap)
   {
      OCCGeometry* occgeom = (OCCGeometry*)geom;
      TopTools_IndexedMapOfShape *occfmap = (TopTools_IndexedMapOfShape *)FaceMap;

      // Copy the face map from the geometry to the given variable
      occfmap->Assign(occgeom->fmap);

      if(occfmap->Extent())
      {
         return CNG_OK;
      }
      else
      {
         return CNG_ERROR;
      }
   }

   // Extract the edge map from the OCC geometry
   // The edge map basically gives an index to each edge in the geometry, 
   // which can be used to access a specific edge
    DLL_HEADER CNg_Result CNg_GetEdgeMap(CNg_OCC_Geometry* geom,
        CNg_OCC_TopTools_IndexedMapOfShape* EdgeMap)
    {
        OCCGeometry* occgeom = (OCCGeometry*)geom;
        TopTools_IndexedMapOfShape* occedgemap = (TopTools_IndexedMapOfShape*)EdgeMap;

        // Copy the edge map from the geometry to the given variable
        occedgemap->Assign(occgeom->emap);

        if (occedgemap->Extent())
        {
            return CNG_OK;
        }
        else
        {
            return CNG_ERROR;
        }
    }
	
   // ------------------ End - OCC Geometry / Meshing Utility Functions ----------------
#endif




   // ------------------ Begin - Meshing Parameters related functions ------------------
   // Constructor for the local nglib meshing parameters class
   // copy of Ng_Meshing_Parameters
   DLL_HEADER CNg_Meshing_Parameters :: CNg_Meshing_Parameters()
   {
      uselocalh = 1;

      maxh = 1000;
      minh = 0.0;

      fineness = 0.5;
      grading = 0.3;

      elementsperedge = 2.0;
      elementspercurve = 2.0;

      closeedgeenable = 0;
      closeedgefact = 2.0;

	  minedgelenenable = 0;
	  minedgelen = 1e-4;

      second_order = 0;
      quad_dominated = 0;

      meshsize_filename = 0;

      optsurfmeshenable = 1;
      optvolmeshenable = 1;

      optsteps_2d = 3;
      optsteps_3d = 3;

      invert_tets = 0;
      invert_trigs = 0;

      check_overlap = 1;
      check_overlapping_boundary = 1;
   }




   // Reset the local meshing parameters to the default values
   // copy of Ng_Meshing_Parameters :: Reset_Parameters
   DLL_HEADER void CNg_Meshing_Parameters :: Reset_Parameters()
   {
      uselocalh = 1;

      maxh = 1000;
      minh = 0;

      fineness = 0.5;
      grading = 0.3;

      elementsperedge = 2.0;
      elementspercurve = 2.0;

      closeedgeenable = 0;
      closeedgefact = 2.0;

  	  minedgelenenable = 0;
	  minedgelen = 1e-4;

      second_order = 0;
      quad_dominated = 0;

      meshsize_filename = 0;

      optsurfmeshenable = 1;
      optvolmeshenable = 1;

      optsteps_2d = 3;
      optsteps_3d = 3;

      invert_tets = 0;
      invert_trigs = 0;

      check_overlap = 1;
      check_overlapping_boundary = 1;
   }




   // copy of Ng_Meshing_Parameters :: Transfer_Parameters()
   DLL_HEADER void CNg_Meshing_Parameters :: Transfer_Parameters()
   {
      mparam.uselocalh = uselocalh;
      
      mparam.maxh = maxh;
      mparam.minh = minh;

      mparam.grading = grading;
      mparam.curvaturesafety = elementspercurve;
      mparam.segmentsperedge = elementsperedge;

      mparam.secondorder = second_order;
      mparam.quad = quad_dominated;

      if (meshsize_filename)
        mparam.meshsizefilename = meshsize_filename;
      else
        mparam.meshsizefilename = "";
      mparam.optsteps2d = optsteps_2d;
      mparam.optsteps3d = optsteps_3d;

      mparam.inverttets = invert_tets;
      mparam.inverttrigs = invert_trigs;

      mparam.checkoverlap = check_overlap;
      mparam.checkoverlappingboundary = check_overlapping_boundary;

   }
   // ------------------ End - Meshing Parameters related functions --------------------






   // ------------------ Begin - Mesh Refinement functions ---------------------
   


#ifdef OCCGEOMETRY
// copy of Ng_OCC_Uniform_Refinement
   DLL_HEADER void CNg_Refine (CNg_Mesh * mesh)
   {

       /*Mesh* me = (Mesh*)mesh;
       me->GetGeometry()->GetRefinement().Refine(*me);
       me->UpdateTopology();*/


       BisectionOptions biopt;
       biopt.usemarkedelements = 1;
       biopt.refine_p = 0;
       biopt.refine_hp = 0;
       Mesh* me = (Mesh*)mesh;

       for (int i = 1; i <= me->GetNSE(); i++)
           me->SurfaceElement(i).SetRefinementFlag(false);

       me->GetGeometry()->GetRefinement().Bisect(*me, biopt);
       //me->GetGeometry()->GetRefinement().Refine(*me);
       me->UpdateTopology();
       me->GetCurvedElements().SetIsHighOrder(false);

	  /* other possible option - check if same result!
      * ( (OCCGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );   
	 */
   }
   
   // Set refinement flag for element
   DLL_HEADER void CNg_SetElementRefinement (CNg_Mesh * mesh, int el_index, bool flag)
   {
       Mesh* me = (Mesh*)mesh;
       if (me->GetDimension() == 3)
       {
           me->VolumeElement(el_index).SetRefinementFlag(flag);
       }
       else
       {
           me->SurfaceElement(el_index).SetRefinementFlag(flag);
       }	 
   }


#endif
   // ------------------ End - Mesh Refinement functions -----------------------
} // End of namespace nglib





// Force linking libinterface to libnglib
#include <../interface/writeuser.hpp>
void MyDummyToForceLinkingLibInterface(Mesh &mesh, NetgenGeometry &geom)
{
  netgen::WriteUserFormat("", mesh, /* geom, */ "");
}

namespace cnglib
{
    DLL_HEADER void CNg_ExportMeshToGmesh2(CNg_Mesh* mesh, const char* filename)
    {
        //Mesh* m = (Mesh*)mesh;
        netgen::Cenos_WriteGmsh2Format(*(Mesh*)mesh, string(filename));
    }

    DLL_HEADER void CNg_ExportMeshToOpenFOAM(CNg_Mesh* mesh, const char* filename)
    {
        //Mesh* m = (Mesh*)mesh;
        netgen::Cenos_WriteOpenFOAM15xFormat(*(Mesh*)mesh, string(filename), false);
    }

    DLL_HEADER CNg_OCC_Geometry* CNg_OCC_ShapeToGeometry(CNg_TopoDS_Shape * s)
    {
        OCCGeometry* geom = netgen::CreateOCCGeometryFromTopoDS((Cenos_TopoDS_Shape*)s);
        return (CNg_OCC_Geometry*)geom;
    }


    DLL_HEADER int CNg_GetEdgeElementIndex(CNg_Mesh* mesh, int num)
    {
//        if (((Mesh*)mesh)->GetDimension() == 3)
//            return ((Mesh*)mesh)->LineSegment(num).edgenr;
//        else
            return ((Mesh*)mesh)->LineSegment(num).epgeominfo[0].edgenr;
    }
	
    DLL_HEADER int CNg_GetSurfaceElementIndex(CNg_Mesh* mesh, int num)
    {
        int fd_id = ((Mesh*)mesh)->SurfaceElement(num).GetIndex();
        return ((Mesh*)mesh)->GetFaceDescriptor(fd_id).SurfNr();
    }
    

    DLL_HEADER int CNg_GetVolumeElementIndex(CNg_Mesh* mesh, int num)
    {
        return ((Mesh*)mesh)->VolumeElement(num).GetIndex();
    }

    DLL_HEADER void CNg_RedirectCout(void* ptr_filestream)
    {
        mycout = (ofstream*)ptr_filestream;
    }

	DLL_HEADER void CNg_CalculateSurfacesOfNode(CNg_Mesh* mesh)
	{
		Mesh* m = (Mesh*)mesh;
		m->CalcSurfacesOfNode();
	}
	
	
    DLL_HEADER void CNg_DumpSegments(CNg_Mesh* mesh, void* ptr_filestream)
    {
        Mesh* m = (Mesh*)mesh;
		ostream* file = (ofstream*)ptr_filestream;
        (*file) << "GetNCD2Names() GetNFD() " << std::endl;

        (*file) << m->GetNCD2Names() << " " << m->GetNFD() << std::endl;
		
        (*file) << "pp1 p2 edgenr singedge_left singedge_right trignum[0] u0 v0 trignum[1] u1 v1 dist1 dist2 si cd2i domin domout tlosurf surfnr1 surfnr2" << std::endl;

        for (int i = 1; i <= m->GetNSeg(); i++)
        {
            Segment seg = m->LineSegment(i);
            (*file) << seg[0] << " " << seg[1] << " " << seg.edgenr << " " << seg.singedge_left <<
                " " << seg.singedge_right << " " << seg.geominfo[0].trignum << " " << seg.epgeominfo[0].u << " " << seg.epgeominfo[0].v << " " <<
                seg.geominfo[1].trignum << " " << seg.epgeominfo[1].u << " " << seg.epgeominfo[1].v << " " << seg.epgeominfo[0].dist << " " << seg.epgeominfo[1].dist 
                << " " << seg.si << " " << seg.cd2i << " " << seg.domin << " " << seg.domout << " " << seg.tlosurf << " " << 
                 seg.surfnr1 << " " << seg.surfnr2   << std::endl;
        }
    }

    
   DLL_HEADER void CNg_AddFaceDescriptor (CNg_Mesh* mesh, int faceId, int dominId, int domOutId, int TLOface)
   {
      Mesh* m = (Mesh*)mesh;
      m->AddFaceDescriptor (FaceDescriptor (faceId, dominId, domOutId, TLOface));
   }
   
   
   DLL_HEADER void CNg_AddEdgeDescriptor (CNg_Mesh* mesh, int edgeId)
   {
      Mesh* m = (Mesh*)mesh;
	  EdgeDescriptor  ed;
	  ed.SetTLOSurface(edgeId);
      m->AddEdgeDescriptor (ed);
   }
}


