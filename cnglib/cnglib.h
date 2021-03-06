#ifndef NGLIB
#define NGLIB

/**************************************************************************/
/* File:   nglib.h                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   7. May. 2000                                                   */
/**************************************************************************/

/*!
   \file cnglib.h
   \brief Library interface to the netgen meshing kernel
   \author Joachim Schoeberl
   \date 7. May 2000

   Modified and adapted for use with Cenos
   This header file provides access to the core functionality of the Netgen 
   Mesher via a library interface, without an interactive User Interface.

   The intention of providing these set of functions is to allow system 
   developers to integrate Netgen into top-level code, to act as the low 
   level mesh generation / optimisation kernel.  
*/

// Philippose - 14.02.2009
// Modifications for creating a DLL in Windows
#ifndef DLL_HEADER
#ifdef WIN32
   #ifdef NGLIB_EXPORTS || nglib_EXPORTS
      #define DLL_HEADER   __declspec(dllexport)
   #else
      #define DLL_HEADER   __declspec(dllimport)
   #endif
#else
  #define DLL_HEADER __attribute__((visibility("default")))
#endif
#endif


// ** Constants used within Netgen *********************
/// Maximum allowed number of nodes per volume element
#define NG_VOLUME_ELEMENT_MAXPOINTS 10

/// Maximum allowed number of nodes per surface element
#define NG_SURFACE_ELEMENT_MAXPOINTS 8



// *** Data-types for accessing Netgen functionality ***
/// Data type for NETGEN mesh
typedef void * CNg_Mesh;

typedef void * CNg_LocalH;

#ifdef OCCGEOMETRY
/// Data type for NETGEN OpenCascade geometry
typedef void * CNg_OCC_Geometry;
typedef void * CNg_OCC_TopTools_IndexedMapOfShape;
typedef void * CNg_TopoDS_Shape;
#endif


// *** Special Enum types used within Netgen ***********
/// Currently implemented surface element types
/// numbering according to gmsh format
enum CNg_Surface_Element_Type 
   { CNG_TRIG = 2, CNG_QUAD = 3, CNG_TRIG6 = 9, CNG_QUAD8 = 16 };

/// Currently implemented volume element types
/// numbering according to gmsh format
enum CNg_Volume_Element_Type 
   { CNG_TET = 4, CNG_HEXA = 5, CNG_PYRAMID = 7, CNG_PRISM = 6, CNG_TET10 = 11 };

/// Values returned by Netgen functions
enum CNg_Result 
   { 
     CNG_ERROR               = -1,   
     CNG_OK                  = 0, 
     CNG_SURFACE_INPUT_ERROR = 1,
     CNG_VOLUME_FAILURE      = 2, 
     CNG_STL_INPUT_ERROR     = 3,
     CNG_SURFACE_FAILURE     = 4,
     CNG_FILE_NOT_FOUND      = 5 
   };



// *** Classes required for use within Netgen **********
/// Netgen Meshing Parameters class
class CNg_Meshing_Parameters 
{
public:
   int uselocalh;                      //!< Switch to enable / disable usage of local mesh size modifiers

   double maxh;                        //!< Maximum global mesh size allowed
   double minh;                        //!< Minimum global mesh size allowed

   double fineness;                    //!< Mesh density: 0...1 (0 => coarse; 1 => fine)
   double grading;                     //!< Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)

   double elementsperedge;             //!< Number of elements to generate per edge of the geometry
   double elementspercurve;            //!< Elements to generate per curvature radius

   int closeedgeenable;                //!< Enable / Disable mesh refinement at close edges
   double closeedgefact;               //!< Factor to use for refinement at close edges (larger => finer)

   int minedgelenenable;			   //!< Enable / Disable user defined minimum edge length for edge subdivision
   double minedgelen;                  //!< Minimum edge length to use while subdividing the edges (default = 1e-4)

   int second_order;                   //!< Generate second-order surface and volume elements
   int quad_dominated;                 //!< Creates a Quad-dominated mesh 

   char * meshsize_filename;           //!< Optional external mesh size file 

   int optsurfmeshenable;              //!< Enable / Disable automatic surface mesh optimization
   int optvolmeshenable;               //!< Enable / Disable automatic volume mesh optimization

   int optsteps_3d;                     //!< Number of optimize steps to use for 3-D mesh optimization
   int optsteps_2d;                     //!< Number of optimize steps to use for 2-D mesh optimization

   // Philippose - 13/09/2010
   // Added a couple more parameters into the meshing parameters list 
   // from Netgen into Nglib
   int invert_tets;                    //!< Invert all the volume elements
   int invert_trigs;                   //!< Invert all the surface triangle elements

   int check_overlap;                  //!< Check for overlapping surfaces during Surface meshing
   int check_overlapping_boundary;     //!< Check for overlapping surface elements before volume meshing


   /*!
      Default constructor for the Mesh Parameters class

      Note: This constructor initialises the variables in the 
      class with the following default values
      - #uselocalh: 1
      - #maxh: 1000.0
      - #fineness: 0.5
      - #grading: 0.3
      - #elementsperedge: 2.0
      - #elementspercurve: 2.0
      - #closeedgeenable: 0
      - #closeedgefact: 2.0
      - #secondorder: 0
      - #meshsize_filename: null
      - #quad_dominated: 0
      - #optsurfmeshenable: 1
      - #optvolmeshenable: 1
      - #optsteps_2d: 3
      - #optsteps_3d: 3
      - #invert_tets: 0
      - #invert_trigs:0 
      - #check_overlap: 1
      - #check_overlapping_boundary: 1
   */
   DLL_HEADER CNg_Meshing_Parameters();



   /*!
       Reset the meshing parameters to their defaults

       This member function resets all the meshing parameters 
       of the object to the default values
   */
   DLL_HEADER void Reset_Parameters();



   /*!
       Transfer local meshing parameters to internal meshing parameters

       This member function transfers all the meshing parameters 
       defined in the local meshing parameters structure of nglib into 
       the internal meshing parameters structure used by the Netgen core
   */
   DLL_HEADER void Transfer_Parameters();
};




// *** Functions Exported by this Library *************

// ------------------------------------------------------------------
// Netgen library initialisation / destruction functions

/*! \brief Initialise the Netgen library and prepare for use

    This function needs to be called by the third-party 
    program before beginning to use the other Netgen 
    specific functions.
*/
DLL_HEADER void CNg_Init ();


/*! \brief Exit the Netgen meshing kernel in a clean manner

    Use this function to exit the meshing sub-system in 
    a clean and orderly manner.
*/
DLL_HEADER void CNg_Exit ();
  

/*! \brief Create a new (and empty) Netgen Mesh Structure

    This function creates a new Netgen Mesh, initialises 
    it, and returns a pointer to the created mesh structure. 

    Use the returned pointer for subsequent operations 
    which involve mesh operations.

    \return Ng_Mesh Pointer to a Netgen Mesh type #Ng_Mesh
*/
DLL_HEADER  CNg_Mesh * CNg_NewMesh ();


/*! \brief Delete an existing Netgen Mesh Structure

    Use this function to delete an existing Netgen mesh 
    structure and release the used memory. 

    \param mesh Pointer to an existing Netgen Mesh structure 
                of type #Ng_Mesh
*/
DLL_HEADER void CNg_DeleteMesh (CNg_Mesh * mesh);


// Copy nodes and elements from origin_mesh into destination_mesh
DLL_HEADER void CNg_MergeMesh(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh, int index = 0);

// Copy segments from origin_mesh into destination_mesh, assign edge index to segments
DLL_HEADER void CNg_CopySegments(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh, int edge_index);

// Copy surface elements from origin_mesh into destination_mesh, assign edge index to segments
DLL_HEADER void CNg_CopySurfaceElements(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh, int face_index);

// Prepare surface meshing
DLL_HEADER void CNg_PrepareSurfaceMeshing(CNg_Mesh* mesh);


// Reverses all segments in mesh
DLL_HEADER void CNg_ReverseSegments(CNg_Mesh* mesh);

// Reverse all faces in mesh
DLL_HEADER void CNg_ReverseFaces(CNg_Mesh* mesh);

/*! \brief Save a Netgen Mesh to disk

    This function allows a generated mesh structure to be saved 
    to disk.

    A Mesh saved using this function, will be written to disk 
    in the Netgen VOL file format.

    \param mesh    Pointer to an existing Netgen Mesh structure 
                   of type #Ng_Mesh
    \param filename Pointer to a character array containing the 
                    name of the file to which the mesh should 
                    be saved
*/
DLL_HEADER void CNg_SaveMesh(CNg_Mesh * mesh, const char* filename);

DLL_HEADER CNg_LocalH CNg_GetLocalH(CNg_Mesh* mesh);

DLL_HEADER void CNg_CopyLocalH(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh);

DLL_HEADER double CNg_GetMaxH(CNg_Mesh* occ_mesh);

DLL_HEADER void CNg_ExportMeshToGmesh2(CNg_Mesh* mesh, const char* filename);

DLL_HEADER void CNg_ExportMeshToOpenFOAM(CNg_Mesh* mesh, const char* filename);


DLL_HEADER CNg_Result CNg_GetSolidMap(CNg_OCC_Geometry* geom,
    CNg_OCC_TopTools_IndexedMapOfShape* SoMap);
	
// Get the face map of an already loaded OCC geometry
DLL_HEADER CNg_Result CNg_GetFaceMap(CNg_OCC_Geometry * geom, 
                                    CNg_OCC_TopTools_IndexedMapOfShape * FMap);
									
DLL_HEADER CNg_Result CNg_GetEdgeMap(CNg_OCC_Geometry* geom,
    CNg_OCC_TopTools_IndexedMapOfShape* EdgeMap);

DLL_HEADER int CNg_GetEdgeElementIndex(CNg_Mesh* mesh, int num);

DLL_HEADER int CNg_GetSurfaceElementIndex(CNg_Mesh* mesh, int num);

DLL_HEADER int CNg_GetVolumeElementIndex(CNg_Mesh* mesh, int num);

DLL_HEADER void CNg_RedirectCout(void* ptr_filestream);

DLL_HEADER void CNg_CalculateSurfacesOfNode(CNg_Mesh* mesh);

DLL_HEADER void CNg_DumpSegments(CNg_Mesh* mesh, void* ptr_filestream);

// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Basic Meshing functions for manually adding points, surface elements 
// and volume elements to a Netgen Mesh structure

/*! \brief Add a point to a given Netgen Mesh Structure

    This function allows points to be directly added to a Netgen 
    mesh structure by providing the co-ordinates.

    Each call to the function allows only one point to be added.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param x    Pointer to an array of type double containing the co-ordinates 
                of the point to be added in the form: \n
                - x[0] = X co-ordinate
                - x[1] = Y co-ordinate
                - x[2] = Z co-ordinate
*/
DLL_HEADER void CNg_AddPoint(CNg_Mesh * mesh, double * x);

DLL_HEADER void CNg_AddSegmentElement(CNg_Mesh * mesh, int pi1, int pi2, int edgeIndex);


/*! \brief Add a surface element to a given Netgen Mesh Structure

    This function allows the top-level code to directly add individual 
    Surface Elements to a Netgen Mesh Structure by providing the type of 
    element to be added and the indices of the points which constitute the 
    element.

    <i>Note:</i> 
    - The points referred to by the surface elements must have been
      added prior to calling this function. 
    - Currently only triangular elements are supported, and the Surface Element 
      Type argument is not used.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param et   Surface Element type provided via the enumerated type 
                #Ng_Surface_Element_Type 
    \param pi   Pointer to an array of integers containing the indices of the 
                points which constitute the surface element being added
*/
DLL_HEADER void CNg_AddSurfaceElementUV(CNg_Mesh* mesh, CNg_Surface_Element_Type et,
       int* pi, int surfIndx, double* uv1, double* uv2, double* uv3);
	   
DLL_HEADER void CNg_AddSurfaceElement(CNg_Mesh * mesh, CNg_Surface_Element_Type et, int * pi, int surfIndx);



/*! \brief Add a volume element to a given Netgen Mesh Structure

    This function allows the top-level code to directly add individual 
    Volume Elements to a Netgen Mesh Structure by providing the type of 
    element to be added and the indices of the points which constitute the 
    element.

    <i>Note:</i> 
    - The points referred to by the volume elements must have been
      added prior to calling this function. 
    - Currently only tetrahedral elements are supported, and the Volume Element 
      Type argument is not used.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param et   Volume Element type provided via the enumerated type 
                #Ng_Volume_Element_Type 
    \param pi   Pointer to an array of integers containing the indices of the 
                points which constitute the volume element being added

*/
DLL_HEADER void CNg_AddVolumeElement(CNg_Mesh * mesh, CNg_Volume_Element_Type et, int * pi, int volIndx);
  
  
  
DLL_HEADER void CNg_AddFaceDescriptor (CNg_Mesh * mesh, int faceId, int dominId, int domOutId, int TLOface);

DLL_HEADER void CNg_AddEdgeDescriptor (CNg_Mesh * mesh, int edgeId);
// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Local Mesh Size restriction / limiting utilities

/*! \brief Apply a global restriction on mesh element size

    This utility allows the user to apply a global mesh element 
    size limitation. 

    During mesh creation, in the absence of an explicit local 
    size restriction around the neighbourhood of a point within 
    the meshing domain, this global size restriction will be 
    utilised.

    <b>Note</b>: This function only limits the <b>Maximum</b> 
    size of an element within the mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param h    Variable of type <i>double</i>, specifying the maximum
                allowable mesh size
*/
DLL_HEADER void CNg_RestrictMeshSizeGlobal (CNg_Mesh * mesh, double h);


/*! \brief Locally restrict the mesh element size at the given point

    Unlike the function #Ng_RestrictMeshSizeGlobal, this function 
    allows the user to locally restrict the maximum allowable mesh 
    size at a given point.

    The point is specified via its three cartesian co-ordinates.

    <b>Note</b>: This function only limits the <b>Maximum</b> size 
    of the elements around the specified point.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param p    Pointer to an Array of type <i>double</i>, containing 
                the three co-ordinates of the point in the form: \n
                - p[0] = X co-ordinate
                - p[1] = Y co-ordinate
                - p[2] = Z co-ordinate
    \param h    Variable of type <i>double</i>, specifying the maximum
                allowable mesh size at that point
*/
DLL_HEADER void CNg_RestrictMeshSizePoint (CNg_Mesh * mesh, double * p, double h);


/*! \brief Locally restrict the mesh element size within a specified box

    Similar to the function #Ng_RestrictMeshSizePoint, this function 
    allows the size of elements within a mesh to be locally limited.

    However, rather than limit the mesh size at a single point, this 
    utility restricts the local mesh size within a 3D Box region, specified 
    via the co-ordinates of the two diagonally opposite points of a cuboid.

    <b>Note</b>: This function only limits the <b>Maximum</b> size 
    of the elements within the specified region.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param pmin Pointer to an Array of type <i>double</i>, containing 
                the three co-ordinates of the first point of the cuboid: \n
                - pmin[0] = X co-ordinate
                - pmin[1] = Y co-ordinate
                - pmin[2] = Z co-ordinate
    \param pmax Pointer to an Array of type <i>double</i>, containing 
                the three co-ordinates of the opposite point of the 
                cuboid: \n
                - pmax[0] = X co-ordinate
                - pmax[1] = Y co-ordinate
                - pmax[2] = Z co-ordinate
    \param h    Variable of type <i>double</i>, specifying the maximum
                allowable mesh size at that point
*/
DLL_HEADER void CNg_RestrictMeshSizeBox (CNg_Mesh * mesh, double * pmin, double * pmax, double h);


DLL_HEADER void CNg_RestrictMeshSizeMesh(CNg_Mesh* orig_mesh, CNg_Mesh* dest_mesh);

DLL_HEADER void CNg_RestrictMeshSizeLocalH(CNg_Mesh* mesh, CNg_LocalH* localh);
// ------------------------------------------------------------------






// ------------------------------------------------------------------
// Basic Mesh information functions

/*! \brief Returns the Number of Points present in the specified Mesh

    Given an already existent Netgen Mesh Structure, this function 
    returns the number of points currently present within the Mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \return 
                Integer Data-type with the number of points in the Mesh
*/
DLL_HEADER int CNg_GetNP (CNg_Mesh * mesh);


/*! \brief Returns the Number of Surface Elements present in the specified Mesh

    Given an already existent Netgen Mesh Structure, this function 
    returns the number of surface elements currently present within 
    the Mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \return 
                Integer Data-type with the number of surface elements in the Mesh
*/
DLL_HEADER int CNg_GetNSE (CNg_Mesh * mesh);


DLL_HEADER int CNg_GetNSeg(CNg_Mesh * mesh);


/*! \brief Returns the Number of Volume Elements present in the specified Mesh

    Given an already existent Netgen Mesh Structure, this function 
    returns the number of volume elements currently present within 
    the Mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \return 
                Integer Data-type with the number of volume elements in the Mesh
*/
DLL_HEADER int CNg_GetNE (CNg_Mesh * mesh);

// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Mesh Topology functions
// Use these functions to extract points, surface / volume elements, 
// perform topological searches, etc..etc...
  
//  Return the Point Coordinates of a specified Point
// The x, y and z co-ordinates are returned in the array pointer as 
// x[0] = x ; x[1] = y ; x[2] = z
DLL_HEADER void CNg_GetPoint (CNg_Mesh * mesh, int num, double * x);


DLL_HEADER int CNg_GetSegment(CNg_Mesh * mesh, int num, int * pi, int * matnum);

   
// return surface and volume element in pi
DLL_HEADER CNg_Surface_Element_Type 
CNg_GetSurfaceElement (CNg_Mesh * mesh, int num, int * pi);

DLL_HEADER CNg_Volume_Element_Type
CNg_GetVolumeElement (CNg_Mesh * mesh, int num, int * pi);

// ------------------------------------------------------------------



#ifdef OCCGEOMETRY

// **********************************************************
// **   OpenCascade Geometry / Meshing Utilities           **
// **********************************************************

// Create new OCC Geometry Object
DLL_HEADER CNg_OCC_Geometry * CNg_OCC_NewGeometry ();

// Delete an OCC Geometry Object
DLL_HEADER CNg_Result CNg_OCC_DeleteGeometry (CNg_OCC_Geometry * geom);

//CENOS
DLL_HEADER CNg_OCC_Geometry * CNg_OCC_ShapeToGeometry(CNg_TopoDS_Shape* shape);

// Set the local mesh size based on geometry / topology
DLL_HEADER CNg_Result CNg_SetLocalMeshSize (CNg_OCC_Geometry * geom,
                                              CNg_Mesh * mesh,
                                              CNg_Meshing_Parameters * mp);

// Set the mesh geometry / topology
DLL_HEADER void CNg_SetGeometry(CNg_OCC_Geometry* geom, CNg_Mesh* mesh);


// Mesh Edges only
DLL_HEADER CNg_Result CNg_DivideEdges(CNg_OCC_Geometry * geom,
                                              CNg_Mesh * mesh,
                                              CNg_Meshing_Parameters * mp);
												 
// Mesh the edges and add Face descriptors to prepare for surface meshing
DLL_HEADER CNg_Result CNg_GenerateEdgeMesh (CNg_OCC_Geometry * geom,
                                              CNg_Mesh * mesh,
                                              CNg_Meshing_Parameters * mp);

// Mesh the surfaces of an OCC geometry
DLL_HEADER CNg_Result CNg_GenerateSurfaceMesh (CNg_OCC_Geometry * geom,
                                                 CNg_Mesh * mesh,
                                                 CNg_Meshing_Parameters * mp); 


/*! \brief Create a 3D Volume Mesh given a Surface Mesh

    After creating a surface mesh, this function can be utilised 
    to automatically generate the corresponding 3D Volume Mesh.

    Mesh generation parameters (such as grading, maximum element size, 
    etc.) are specified via the meshing parameters class which also 
    needs to be passed to this function.

    <b>Note</b>: Currently, Netgen generates pure tetrahedral volume 
    meshes.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param mp   Pointer to a copy of the Meshing Parameters class
                (#Ng_Meshing_Parameters), filled up with the 
                required values

    \return Ng_Result Status of the Mesh Generation routine. More 
                      details regarding the return value can be 
                      found in the description of #Ng_Result
*/
DLL_HEADER CNg_Result CNg_GenerateVolumeMesh (CNg_Mesh * mesh, CNg_Meshing_Parameters * mp);


/*! \brief Create a boundary layers for given volume mesh

    After creating a volume mesh, this function can be utilised
    to generate boundary layers on given surfaces.


    \param mesh         Pointer to an existing Netgen Mesh structure of
                        type #Ng_Mesh
    \param surfid_arr   array of surface ids where layers will be created

    \param surfid_count   number of faces in surfid_arr

    \param heights_arr   array of layer heights

    \param heights_count   number of layers

    \return Ng_Result Status of the Mesh Generation routine. More
                      details regarding the return value can be
                      found in the description of #Ng_Result
*/
DLL_HEADER CNg_Result CNg_GenerateBoundaryLayer(CNg_Mesh* mesh,
    int* surfid_arr, int surfid_count,
    double* heights_arr, int heights_count);


DLL_HEADER CNg_Result CNg_GenerateBoundaryLayer2(CNg_Mesh* mesh,
    int dom_nr, double* heights_arr, int heights_count);
#endif // OCCGEOMETRY



// **********************************************************
// **   Mesh refinement algorithms                         **
// **********************************************************


#ifdef OCCGEOMETRY
DLL_HEADER void CNg_Refine(CNg_Mesh * mesh);

DLL_HEADER void CNg_SetElementRefinement(CNg_Mesh* mesh, int el_index, bool b);

#endif


#endif // NGLIB
