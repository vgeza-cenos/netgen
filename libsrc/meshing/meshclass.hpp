#ifndef MESHCLASS
#define MESHCLASS

/**************************************************************************/
/* File:   meshclass.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

/*
  The mesh class
*/

namespace netgen
{
  enum resthtype { RESTRICTH_FACE, RESTRICTH_EDGE, 
		   RESTRICTH_SURFACEELEMENT, RESTRICTH_POINT, RESTRICTH_SEGMENT };

  class HPRefElement;


  /// 2d/3d mesh
  class Mesh
  {
  public:
    typedef ::netgen::T_POINTS T_POINTS;
    typedef Array<Element, 0, size_t> T_VOLELEMENTS;
    // typedef Array<Element2d, 0, SurfaceElementIndex> T_SURFELEMENTS;
    typedef Array<Element2d, 0, size_t> T_SURFELEMENTS;

  private:
    /// point coordinates
    T_POINTS points;

#ifdef PARALLEL
    // The communicator for this mesh. (more or less dummy for now!)  
    MPI_Comm comm;
#endif
    
    /// line-segments at edges
    Array<Segment, 0, size_t> segments;
    /// surface elements, 2d-inner elements
    T_SURFELEMENTS surfelements;
    /// volume elements
    T_VOLELEMENTS volelements;
    /// points will be fixed forever
    Array<PointIndex> lockedpoints;


    /// surface indices at boundary nodes
    // TABLE<int,PointIndex::BASE> surfacesonnode;
    /// boundary edges  (1..normal bedge, 2..segment)
    INDEX_2_CLOSED_HASHTABLE<int> * boundaryedges;
    ///
    INDEX_2_CLOSED_HASHTABLE<int> * segmentht;
    ///
    INDEX_3_CLOSED_HASHTABLE<int> * surfelementht;

    /// faces of rest-solid
    Array<Element2d> openelements;
    /// open segmenets for surface meshing  
    Array<Segment> opensegments;



    /**
       Representation of local mesh-size h
    */
    LocalH * lochfunc;
    ///
    double hglob;
    ///
    double hmin;
    ///
    Array<double> maxhdomain;
  
    /**
       the face-index of the surface element maps into
       this table.
    */
    Array<FaceDescriptor> facedecoding;

  
    /**
       the edge-index of the line element maps into
       this table.
    */
    Array<EdgeDescriptor> edgedecoding;

    /// sub-domain materials 
    Array<string*> materials;

    /// labels for boundary conditions
    Array<string*> bcnames;

    /// labels for co dim 2 bboundary conditions
    Array<string*> cd2names;

    /// labels for co dim 3 bbboundary conditions
    Array<string*> cd3names;

    /// Periodic surface, close surface, etc. identifications
    Identifications * ident;


    /// number of vertices (if < 0, use np)
    int numvertices;

    /// geometric search tree for interval intersection search
    BoxTree<3> * elementsearchtree;
    /// time stamp for tree
    mutable int elementsearchtreets;

    /// element -> face, element -> edge etc ...
    MeshTopology topology;
    /// methods for high order elements
    class CurvedElements * curvedelems;

    /// nodes identified by close points 
    class AnisotropicClusters * clusters;

    /// space dimension (2 or 3)
    int dimension;
  
    /// changed by every minor modification (addpoint, ...)
    int timestamp;
    /// changed after finishing global algorithm (improve, ...)
    int majortimestamp;

    /// mesh access semaphors.
    NgMutex mutex;
    /// mesh access semaphors.
    NgMutex majormutex;

    SYMBOLTABLE< Array<int>* > userdata_int;
    SYMBOLTABLE< Array<double>* > userdata_double; 


    mutable Array< Point3d > pointcurves;
    mutable Array<int> pointcurves_startpoint;
    mutable Array<double> pointcurves_red,pointcurves_green,pointcurves_blue;


    /// start element for point search (GetElementOfPoint)
    mutable int ps_startelement;


#ifdef PARALLEL
    /// connection to parallel meshes
    class ParallelMeshTopology * paralleltop;

#endif

    
    shared_ptr<NetgenGeometry> geometry;

  private:
    void BuildBoundaryEdges(void);

  public:
    bool PointContainedIn2DElement(const Point3d & p,
				   double lami[3],
				   const int element,
				   bool consider3D = false) const;
    bool PointContainedIn3DElement(const Point3d & p,
				   double lami[3],
				   const int element) const;
    bool PointContainedIn3DElementOld(const Point3d & p,
				      double lami[3],
				      const int element) const;

  public:

    // store coarse mesh before hp-refinement
    Array<HPRefElement> * hpelements;
    Mesh * coarsemesh;
  
  
    /// number of refinement levels
    int mglevels;
    /// refinement hierarchy
    Array<PointIndices<2>,PointIndex::BASE> mlbetweennodes;
    /// parent element of volume element
    Array<int> mlparentelement;
    /// parent element of surface element
    Array<int> mlparentsurfaceelement;



    ///
    DLL_HEADER Mesh();
    ///
    DLL_HEADER ~Mesh();

    Mesh & operator= (const Mesh & mesh2);
  
    ///
    DLL_HEADER void DeleteMesh();
  
    ///
    void ClearSurfaceElements();

    ///
    DLL_HEADER void ClearVolumeElements()
    {
      volelements.SetSize(0); 
      timestamp = NextTimeStamp();
    }

    ///
    DLL_HEADER void ClearSegments()
    { 
      segments.SetSize(0); 
      timestamp = NextTimeStamp();
    }
    
    ///
    bool TestOk () const;

    void SetAllocSize(int nnodes, int nsegs, int nsel, int nel);
    

    DLL_HEADER PointIndex AddPoint (const Point3d & p, int layer = 1);
    DLL_HEADER PointIndex AddPoint (const Point3d & p, int layer, POINTTYPE type);

    int GetNP () const { return points.Size(); }

    // [[deprecated("Use Point(PointIndex) instead of int !")]]        
    MeshPoint & Point(int i) { return points.Elem(i); }
    MeshPoint & Point(PointIndex pi) { return points[pi]; }
    // [[deprecated("Use Point(PointIndex) instead of int !")]]            
    const MeshPoint & Point(int i) const { return points.Get(i); }
    const MeshPoint & Point(PointIndex pi) const { return points[pi]; }

    const MeshPoint & operator[] (PointIndex pi) const { return points[pi]; }
    MeshPoint & operator[] (PointIndex pi) { return points[pi]; }

    const T_POINTS & Points() const { return points; }
    T_POINTS & Points() { return points; }


    DLL_HEADER SegmentIndex AddSegment (const Segment & s);
    void DeleteSegment (int segnr)
    {
      segments.Elem(segnr)[0].Invalidate();
      segments.Elem(segnr)[1].Invalidate();
    }
    /*
    void FullDeleteSegment (int segnr)  // von wem ist das ???
    {
      segments.Delete(segnr-PointIndex::BASE);
    }
    */

    int GetNSeg () const { return segments.Size(); }
    // [[deprecated("Use LineSegment(SegmentIndex) instead of int !")]]                
    Segment & LineSegment(int i) { return segments.Elem(i); }
    // [[deprecated("Use LineSegment(SegmentIndex) instead of int !")]]                    
    const Segment & LineSegment(int i) const { return segments.Get(i); }

    Segment & LineSegment(SegmentIndex si) { return segments[si]; }
    const Segment & LineSegment(SegmentIndex si) const { return segments[si]; }
    const Segment & operator[] (SegmentIndex si) const { return segments[si]; }
    Segment & operator[] (SegmentIndex si) { return segments[si]; }

    /*
    const Array<Segment> & LineSegments() const { return segments; }
    Array<Segment> & LineSegments() { return segments; }
    */
    const auto & LineSegments() const { return segments; }
    auto & LineSegments() { return segments; }
    
    Array<Element0d> pointelements;  // only via python interface

    DLL_HEADER SurfaceElementIndex AddSurfaceElement (const Element2d & el);
    // write to pre-allocated container, thread-safe
    DLL_HEADER void SetSurfaceElement (SurfaceElementIndex sei, const Element2d & el);
    
    // [[deprecated("Use DeleteSurfaceElement(SurfaceElementIndex) instead of int !")]]
    void DeleteSurfaceElement (int eli)
    { 
      surfelements.Elem(eli).Delete();
      surfelements.Elem(eli).PNum(1).Invalidate();
      surfelements.Elem(eli).PNum(2).Invalidate();
      surfelements.Elem(eli).PNum(3).Invalidate();
      timestamp = NextTimeStamp();
    }

    void DeleteSurfaceElement (SurfaceElementIndex eli)
    {
      for (auto & p : surfelements[eli].PNums()) p.Invalidate();
      surfelements[eli].Delete();
      timestamp = NextTimeStamp();
    }

    int GetNSE () const { return surfelements.Size(); }

    // [[deprecated("Use SurfaceElement(SurfaceElementIndex) instead of int !")]]    
    Element2d & SurfaceElement(int i) { return surfelements.Elem(i); }
    // [[deprecated("Use SurfaceElement(SurfaceElementIndex) instead of int !")]]        
    const Element2d & SurfaceElement(int i) const { return surfelements.Get(i); }
    Element2d & SurfaceElement(SurfaceElementIndex i) { return surfelements[i]; }
    const Element2d & SurfaceElement(SurfaceElementIndex i) const { return surfelements[i]; }

    const Element2d & operator[] (SurfaceElementIndex ei) const
    { return surfelements[ei]; }
    Element2d & operator[] (SurfaceElementIndex ei)
    { return surfelements[ei]; }

    const T_SURFELEMENTS & SurfaceElements() const { return surfelements; }
    T_SURFELEMENTS & SurfaceElements() { return surfelements; }

  
    DLL_HEADER void RebuildSurfaceElementLists ();
    DLL_HEADER void GetSurfaceElementsOfFace (int facenr, Array<SurfaceElementIndex> & sei) const;

    DLL_HEADER ElementIndex AddVolumeElement (const Element & el);
    // write to pre-allocated container, thread-safe
    DLL_HEADER void SetVolumeElement (ElementIndex sei, const Element & el);

    int GetNE () const { return volelements.Size(); }

    // [[deprecated("Use VolumeElement(ElementIndex) instead of int !")]]    
    Element & VolumeElement(int i) { return volelements.Elem(i); }
    // [[deprecated("Use VolumeElement(ElementIndex) instead of int !")]]        
    const Element & VolumeElement(int i) const { return volelements.Get(i); }
    Element & VolumeElement(ElementIndex i) { return volelements[i]; }
    const Element & VolumeElement(ElementIndex i) const { return volelements[i]; }

    const Element & operator[] (ElementIndex ei) const 
    { return volelements[ei]; }
    Element & operator[] (ElementIndex ei)
    { return volelements[ei]; }



    ELEMENTTYPE ElementType (ElementIndex i) const 
    { return (volelements[i].flags.fixed) ? FIXEDELEMENT : FREEELEMENT; }

    const auto & VolumeElements() const { return volelements; }
    auto & VolumeElements() { return volelements; }

    ///
    DLL_HEADER double ElementError (int eli, const MeshingParameters & mp) const;

    /// 
    DLL_HEADER void AddLockedPoint (PointIndex pi);
    ///
    void ClearLockedPoints ();

    const auto & LockedPoints() const { return lockedpoints; }

    /// Returns number of domains
    DLL_HEADER int GetNDomains() const;
    ///
    int GetDimension() const { return dimension; }
    void SetDimension (int dim) { dimension = dim; }

    /// sets internal tables
    DLL_HEADER void CalcSurfacesOfNode ();

    /// additional (temporarily) fix points 
    void FixPoints (const BitArray & fixpoints);

    /**
       finds elements without neighbour and
       boundary elements without inner element.
       Results are stored in openelements.
       if dom == 0, all sub-domains, else subdomain dom */
    DLL_HEADER void FindOpenElements (int dom = 0);

  
    /**
       finds segments without surface element,
       and surface elements without neighbours.
       store in opensegmentsy
    */
    DLL_HEADER void FindOpenSegments (int surfnr = 0);
    /**
       remove one layer of surface elements
    */
    DLL_HEADER void RemoveOneLayerSurfaceElements ();


    int GetNOpenSegments () { return opensegments.Size(); }
    const Segment & GetOpenSegment (int nr) { return opensegments.Get(nr); }
  
    /**
       Checks overlap of boundary
       return == 1, iff overlap
    */
    DLL_HEADER int CheckOverlappingBoundary ();
    /**
       Checks consistent boundary
       return == 0, everything ok
    */
    DLL_HEADER int CheckConsistentBoundary () const;

    /*
      checks element orientation
    */
    DLL_HEADER int CheckVolumeMesh () const;


    /**
       finds average h of surface surfnr if surfnr > 0,
       else of all surfaces.
    */
    DLL_HEADER double AverageH (int surfnr = 0) const;
    /// Calculates localh 
    DLL_HEADER void CalcLocalH (double grading);
    ///
    DLL_HEADER void SetLocalH (netgen::Point<3> pmin, netgen::Point<3> pmax, double grading);
    ///
    DLL_HEADER void RestrictLocalH (const Point3d & p, double hloc);
    ///
    DLL_HEADER void RestrictLocalHLine (const Point3d & p1, const Point3d & p2, 
			     double hloc);
    /// number of elements per radius
    DLL_HEADER void CalcLocalHFromSurfaceCurvature(double grading, double elperr);
    ///
    DLL_HEADER void CalcLocalHFromPointDistances(double grading);
    ///
    DLL_HEADER void RestrictLocalH (resthtype rht, int nr, double loch);
    ///
    DLL_HEADER void LoadLocalMeshSize (const string & meshsizefilename);
    ///
    DLL_HEADER void SetGlobalH (double h);
    ///
	DLL_HEADER void SetMinimalH (double h);
    ///
	DLL_HEADER double MaxHDomain (int dom) const;
    ///
	DLL_HEADER void SetMaxHDomain (const Array<double> & mhd);
    ///
    DLL_HEADER double GetH (const Point3d & p) const;
    ///
    double GetMinH (const Point3d & pmin, const Point3d & pmax);
    ///
    bool HasLocalHFunction () { return lochfunc != nullptr; }
    ///
    LocalH & LocalHFunction () { return * lochfunc; }
    ///
    bool LocalHFunctionGenerated(void) const { return (lochfunc != NULL); }

    /// Find bounding box
    DLL_HEADER void GetBox (Point3d & pmin, Point3d & pmax, int dom = -1) const;

    /// Find bounding box of points of typ ptyp or less
    DLL_HEADER void GetBox (Point3d & pmin, Point3d & pmax, POINTTYPE ptyp ) const;

    ///
    int GetNOpenElements() const
    { return openelements.Size(); }
    ///
    const Element2d & OpenElement(int i) const
    { return openelements.Get(i); }


    /// are also quads open elements
    bool HasOpenQuads () const;

    /// split into connected pieces
	DLL_HEADER void SplitIntoParts ();

    /// 
	DLL_HEADER void SplitSeparatedFaces ();

    /// Refines mesh and projects points to true surface
    // void Refine (int levels, const CSGeometry * geom);
  
    
    bool BoundaryEdge (PointIndex pi1, PointIndex pi2) const
    {
      if(!boundaryedges)
	const_cast<Mesh *>(this)->BuildBoundaryEdges();

      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return boundaryedges->Used (i2);
    }

    bool IsSegment (PointIndex pi1, PointIndex pi2) const
    {
      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return segmentht->Used (i2);
    }

    SegmentIndex SegmentNr (PointIndex pi1, PointIndex pi2) const
    {
      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return segmentht->Get (i2);
    }


    /**
       Remove unused points. etc.
    */
    DLL_HEADER void Compress ();

    /// first vertex has lowest index
    void OrderElements(); 

    ///
	DLL_HEADER void Save (ostream & outfile) const;
    ///
	DLL_HEADER void Load (istream & infile);
    ///
	DLL_HEADER void Merge (istream & infile, const int surfindex_offset = 0);
    ///
	DLL_HEADER void Save (const string & filename) const;
    ///
	DLL_HEADER void Load (const string & filename);
    ///
	DLL_HEADER void Merge (const string & filename, const int surfindex_offset = 0);


    DLL_HEADER void DoArchive (Archive & archive);
    ///
	DLL_HEADER void ImproveMesh (const MeshingParameters & mp, OPTIMIZEGOAL goal = OPT_QUALITY);

    ///
    void ImproveMeshJacobian (const MeshingParameters & mp, OPTIMIZEGOAL goal = OPT_QUALITY, const BitArray * usepoint = NULL);
    ///
    void ImproveMeshJacobianOnSurface (const MeshingParameters & mp,
				       const BitArray & usepoint, 
				       const Array< Vec<3>* > & nv,
				       OPTIMIZEGOAL goal = OPT_QUALITY,
				       const Array< Array<int,PointIndex::BASE>* > * idmaps = NULL);
    /**
       free nodes in environment of openelements 
       for optimiztion
    */
    void FreeOpenElementsEnvironment (int layers);

    ///
    bool LegalTet (Element & el) const
    {
      if (el.IllegalValid())
	return !el.Illegal();
      return LegalTet2 (el);
    }
    ///
    bool LegalTet2 (Element & el) const;


    ///
    bool LegalTrig (const Element2d & el) const;
    /**
       if values non-null, return values in 4-double array:
       triangle angles min/max, tetangles min/max
       if null, output results on cout
    */
	DLL_HEADER void CalcMinMaxAngle (double badellimit, double * retvalues = NULL);

    /*
      Marks elements which are dangerous to refine
      return: number of illegal elements
    */
	DLL_HEADER int MarkIllegalElements ();

    /// orient surface mesh, for one sub-domain only
	DLL_HEADER void SurfaceMeshOrientation ();

    /// convert mixed element mesh to tet-mesh
	DLL_HEADER void Split2Tets();


    /// build box-search tree
    DLL_HEADER void BuildElementSearchTree ();

    void SetPointSearchStartElement(const int el) const {ps_startelement = el;}

    /// gives element of point, barycentric coordinates
    int GetElementOfPoint (const netgen::Point<3> & p,
			   double * lami,
			   bool build_searchtree = 0,
			   const int index = -1,
			   const bool allowindex = true) const;
    int GetElementOfPoint (const netgen::Point<3> & p,
			   double * lami,
			   const Array<int> * const indices,
			   bool build_searchtree = 0,
			   const bool allowindex = true) const;
    int GetSurfaceElementOfPoint (const netgen::Point<3> & p,
				  double * lami,
				  bool build_searchtree = 0,
				  const int index = -1,
				  const bool allowindex = true) const;
    int GetSurfaceElementOfPoint (const netgen::Point<3> & p,
				  double * lami,
				  const Array<int> * const indices,
				  bool build_searchtree = 0,
				  const bool allowindex = true) const;

    /// give list of vol elements which are int the box(p1,p2)
    void GetIntersectingVolEls(const Point3d& p1, const Point3d& p2, 
			       Array<int> & locels) const;

    ///
    int AddFaceDescriptor(const FaceDescriptor& fd)
    { facedecoding.Append(fd); return facedecoding.Size(); }

    int AddEdgeDescriptor(const EdgeDescriptor & fd)
    { edgedecoding.Append(fd); return edgedecoding.Size() - 1; }

    ///
    DLL_HEADER void SetMaterial (int domnr, const string & mat);
    ///
    DLL_HEADER const string & GetMaterial (int domnr) const;
    DLL_HEADER static string defaultmat;
    const string * GetMaterialPtr (int domnr) const // 1-based
    {
      return domnr <= materials.Size() ? materials.Get(domnr) : &defaultmat;
    }
    
    DLL_HEADER void SetNBCNames ( int nbcn );

    DLL_HEADER void SetBCName ( int bcnr, const string & abcname );

    DLL_HEADER const string & GetBCName ( int bcnr ) const;

    DLL_HEADER void SetNCD2Names (int ncd2n);
    DLL_HEADER void SetCD2Name (int cd2nr, const string & abcname);

    DLL_HEADER const string & GetCD2Name (int cd2nr ) const;
    DLL_HEADER static string cd2_default_name;
    string * GetCD2NamePtr (int cd2nr ) const
    {
      if (cd2nr < cd2names.Size() && cd2names[cd2nr]) return cd2names[cd2nr];
      return &cd2_default_name;
    }
    size_t GetNCD2Names() const { return cd2names.Size(); }

    DLL_HEADER void SetNCD3Names (int ncd3n);
    DLL_HEADER void SetCD3Name (int cd3nr, const string & abcname);

    DLL_HEADER const string & GetCD3Name (int cd3nr ) const;
    DLL_HEADER static string cd3_default_name;
    string * GetCD3NamePtr (int cd3nr ) const
    {
      if (cd3nr < cd3names.Size() && cd3names[cd3nr]) return cd3names[cd3nr];
      return &cd3_default_name;
    }
    size_t GetNCD3Names() const { return cd3names.Size(); }

    DLL_HEADER static string default_bc;
    string * GetBCNamePtr (int bcnr) const
    { return (bcnr < bcnames.Size() && bcnames[bcnr]) ? bcnames[bcnr] : &default_bc; }

    ///
    void ClearFaceDescriptors()
    { facedecoding.SetSize(0); }

    ///
    int GetNFD () const
    { return facedecoding.Size(); }

    const FaceDescriptor & GetFaceDescriptor (int i) const
    { return facedecoding.Get(i); }

    const EdgeDescriptor & GetEdgeDescriptor (int i) const
    { return edgedecoding[i]; }


    ///
    FaceDescriptor & GetFaceDescriptor (int i) 
    { return facedecoding.Elem(i); }

    // #ifdef NONE
    //   /*
    //     Identify points pi1 and pi2, due to
    //     identification nr identnr
    //   */
    //   void AddIdentification (int pi1, int pi2, int identnr);

    //   int GetIdentification (int pi1, int pi2) const;
    //   int GetIdentificationSym (int pi1, int pi2) const;
    //   ///
    //   INDEX_2_HASHTABLE<int> & GetIdentifiedPoints () 
    //   { 
    //     return *identifiedpoints; 
    //   }

    //   ///
    //   void GetIdentificationMap (int identnr, Array<int> & identmap) const;
    //   ///
    //   void GetIdentificationPairs (int identnr, Array<INDEX_2> & identpairs) const;
    //   ///
    //   int GetMaxIdentificationNr () const
    //   { 
    //     return maxidentnr; 
    //   }
    // #endif

    /// return periodic, close surface etc. identifications
    Identifications & GetIdentifications () { return *ident; }
    /// return periodic, close surface etc. identifications
    const Identifications & GetIdentifications () const { return *ident; }
    ///
    bool HasIdentifications() const { return ident != nullptr; }

    void InitPointCurve(double red = 1, double green = 0, double blue = 0) const;
    void AddPointCurvePoint(const Point3d & pt) const;
    int GetNumPointCurves(void) const;
    int GetNumPointsOfPointCurve(int curve) const;
    Point3d & GetPointCurvePoint(int curve, int n) const;
    void GetPointCurveColor(int curve, double & red, double & green, double & blue) const;




    /// find number of vertices
    void ComputeNVertices ();
    /// number of vertices (no edge-midpoints)
    int GetNV () const;
    /// remove edge points
    void SetNP (int np);

  


    DLL_HEADER bool PureTrigMesh (int faceindex = 0) const;
    DLL_HEADER bool PureTetMesh () const;


    const MeshTopology & GetTopology () const
    { return topology; }

    DLL_HEADER void UpdateTopology (TaskManager tm = &DummyTaskManager,
                                    Tracer tracer = &DummyTracer);
  
    class CurvedElements & GetCurvedElements () const
    { return *curvedelems; }
    
    DLL_HEADER void BuildCurvedElements  (const class Refinement * ref, int aorder, bool arational = false);
    DLL_HEADER void BuildCurvedElements  (int aorder);

    const class AnisotropicClusters & GetClusters () const
    { return *clusters; }


    class CSurfaceArea
    {
      const Mesh & mesh;
      bool valid;
      double area;
    public:
      CSurfaceArea (const Mesh & amesh) 
	: mesh(amesh), valid(false) { ; }

      void Add (const Element2d & sel)
      {
	if (sel.GetNP() == 3)
	  area += Cross ( mesh[sel[1]]-mesh[sel[0]],
			  mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;
	else
	  area += Cross (Vec3d (mesh[sel.PNum(1)], mesh[sel.PNum(3)]),
			 Vec3d (mesh[sel.PNum(1)], mesh[sel.PNum(4)])).Length() / 2;;
      }
      void ReCalc ()
      {
	area = 0;
	for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	  Add (mesh[sei]);
	valid = true;
      }

      operator double () const { return area; }
      bool Valid() const { return valid; }
    };

    CSurfaceArea surfarea;
    CSurfaceArea & SurfaceArea() { return surfarea; }
    const CSurfaceArea & SurfaceArea() const { return surfarea; }



    int GetTimeStamp() const { return timestamp; }
    void SetNextTimeStamp() 
    { timestamp = NextTimeStamp(); }

    int GetMajorTimeStamp() const { return majortimestamp; }
    void SetNextMajorTimeStamp() 
    { majortimestamp = timestamp = NextTimeStamp(); }


    /// return mutex
    NgMutex & Mutex ()   { return mutex; }
    NgMutex & MajorMutex ()   { return majormutex; }


    shared_ptr<NetgenGeometry> GetGeometry() const
    { 
      return geometry; 
    }
    void SetGeometry (shared_ptr<NetgenGeometry> geom) 
    {
      geometry = geom;
    }

    ///
    void SetUserData(const char * id, Array<int> & data);
    ///
    bool GetUserData(const char * id, Array<int> & data, int shift = 0) const;
    ///
    void SetUserData(const char * id, Array<double> & data);
    ///
    bool GetUserData(const char * id, Array<double> & data, int shift = 0) const;

    ///
    friend void OptimizeRestart (Mesh & mesh3d);
    ///
    void PrintMemInfo (ostream & ost) const;
    /// 
    friend class Meshing3;


    enum GEOM_TYPE { NO_GEOM = 0, GEOM_2D = 1, GEOM_CSG = 10, GEOM_STL = 11, GEOM_OCC = 12, GEOM_ACIS = 13 };
    GEOM_TYPE geomtype;
  


#ifdef PARALLEL
    /// returns parallel topology
    class ParallelMeshTopology & GetParallelTopology () const
    { return *paralleltop; }


    /// distributes the master-mesh to local meshes
    void Distribute ();
    void Distribute (Array<int> & volume_weights, Array<int> & surface_weights, 
		     Array<int> & segment_weights);


    /// find connection to parallel meshes
    //   void FindExchangePoints () ;

    //   void FindExchangeEdges ();
    //   void FindExchangeFaces ();

    /// use metis to decompose master mesh 
    void ParallelMetis (); //  Array<int> & neloc );
    void ParallelMetis (Array<int> & volume_weights, Array<int> & surface_weights, 
			Array<int> & segment_weights); 

    void PartHybridMesh (); //  Array<int> & neloc );
    void PartDualHybridMesh (); //  Array<int> & neloc );
    void PartDualHybridMesh2D ();  // ( Array<int> & neloc );


    /// send mesh from master to local procs
    void SendRecvMesh ();

    /// send mesh to parallel machine, keep global mesh at master 
    void SendMesh ( ) const;   // Mesh * mastermesh, Array<int> & neloc) const;
    /// loads a mesh sent from master processor
    void ReceiveParallelMesh ();

#endif


  };

  inline ostream& operator<<(ostream& ost, const Mesh& mesh)
  {
    ost << "mesh: " << endl;
    mesh.Save(ost);
    return ost;
  }

}

#endif


