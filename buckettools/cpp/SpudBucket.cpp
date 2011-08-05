
#include "SpudBucket.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "BoostTypes.h"
#include "BucketDolfinBase.h"
//#include "GenericDetectors.h"
//#include "PointDetectors.h"
//#include "PythonDetectors.h"
//#include "PythonExpression.h"
//#include "InitialConditionExpression.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

// Default constructor for spudbucket derived class
SpudBucket::SpudBucket() : Bucket()
{
  // Do nothing
}

// Default constructor for spudbucket derived class
SpudBucket::SpudBucket(std::string name) : Bucket(name)
{
  // Do nothing
}

// Default constructor for spudbucket derived class
SpudBucket::SpudBucket(std::string name, std::string optionpath) : optionpath_(optionpath), Bucket(name)
{
  // Do nothing
}

// Default destructor for spudbucket derived class
SpudBucket::~SpudBucket()
{
  empty_();
}

// Fill the bucket with stuff based on options tree provided by Spud
void SpudBucket::fill()
{
  // A buffer to put option paths in
  std::stringstream buffer;
  Spud::OptionError serr;
  
  // Get the dimension to pass it down to all systems
  // (we currently assume this is the length of some
  //  of things, which may need to be generalized in the
  //  future)
  buffer.str(""); buffer << optionpath() << "/geometry/dimension";  
  serr = Spud::get_option(buffer.str(), dimension_); spud_err(buffer.str(), serr);

  // Put meshes into the bucket
  buffer.str(""); buffer << optionpath() << "/geometry/mesh";  
  int nmeshes = Spud::option_count(buffer.str());
  for (uint i = 0; i<nmeshes; i++)
  {
    buffer.str(""); buffer << optionpath() << "/geometry/mesh[" << i << "]";
    meshes_fill_(buffer.str());
  }
  
  // Having just registered the system uflsymbols, quickly loop through
  // the rest of the fields recording their symbol and name
  buffer.str(""); buffer << optionpath() << "/system";
  int nsystems = Spud::option_count(buffer.str());
  for (uint i = 0; i<nsystems; i++)
  {
    buffer.str(""); buffer << optionpath() << "/system[" << i << "]";
    uflnames_fill_(buffer.str());
  }

  // Put systems into the bucket
  for (uint i = 0; i<nsystems; i++)
  {
    buffer.str(""); buffer << optionpath() << "/system[" << i << "]";
    systems_fill_(buffer.str());
  }

  for (SystemBucket_it sys_it = systems_begin(); sys_it != systems_end(); sys_it++)
  {
    (*boost::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).funccoeffs_fill();
  }
  
  uflsymbols_fill_();

  for (SystemBucket_it sys_it = systems_begin(); sys_it != systems_end(); sys_it++)
  {
    (*(*sys_it).second).attach_and_initialize();
  }
  
//  // Put detectors in the bucket
//  detectors_fill_();

}

// Insert a mesh (with given optionpath) into the bucket
void SpudBucket::meshes_fill_(const std::string &optionpath)
{
  // A buffer to put option paths in
  std::stringstream buffer;
  Spud::OptionError serr;
  
  // Get the name of the mesh
  std::string meshname;
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), meshname); spud_err(buffer.str(), serr);
  
  std::string source;
  buffer.str(""); buffer << optionpath << "/source/name";
  serr = Spud::get_option(buffer.str(), source); spud_err(buffer.str(), serr);

  Mesh_ptr mesh;

  if (source=="File")
  {
    // Get the name of the mesh file
    std::string basename;
    buffer.str(""); buffer << optionpath << "/source/file";
    serr = Spud::get_option(buffer.str(), basename); spud_err(buffer.str(), serr);

    // Use DOLFIN to read in the mesh
    std::stringstream filename;
    filename.str(""); filename << basename << ".xml";
    mesh.reset(new dolfin::Mesh(filename.str()));
    (*mesh).init();

    std::ifstream file;
    
    // Register the mesh functions (in dolfin::MeshData associated with the mesh - saves us having a separate set of mesh functions in
    // the Bucket, which would need to be associated with a particular mesh in the case of multiple meshes!):
    // - for the edge ids
    filename.str(""); filename << basename << "_edge_subdomain.xml";
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)
    {
      file.close();
      MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains");
      dolfin::MeshFunction<dolfin::uint> edgeids(*mesh, filename.str());
      *meshfuncedgeids = edgeids;
    }

    // - for the cell ids
    filename.str(""); filename << basename << "_cell_subdomain.xml";
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)
    {
      file.close();
      MeshFunction_uint_ptr meshfunccellids = (*mesh).data().create_mesh_function("cell_domains");
      dolfin::MeshFunction<dolfin::uint> cellids(*mesh, filename.str());
      *meshfunccellids = cellids;
    }

  }
  else if (source=="UnitInterval")
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitInterval(cells) );

    MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains", 0);
    Side left(0, 0.0);
    Side right(0, 1.0);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
  }
  else if (source=="Interval")
  {
    double leftx;
    buffer.str(""); buffer << optionpath << "/source/left";
    serr = Spud::get_option(buffer.str(), leftx); spud_err(buffer.str(), serr);

    double rightx;
    buffer.str(""); buffer << optionpath << "/source/right";
    serr = Spud::get_option(buffer.str(), rightx); spud_err(buffer.str(), serr);

    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::Interval(cells, leftx, rightx) );

    MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains", 0);
    Side left(0, leftx);
    Side right(0, rightx);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
  }
  else if (source=="UnitSquare")
  {
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitSquare(cells[0], cells[1], diagonal) );

    MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains", 1);
    Side left(0, 0.0);
    Side right(0, 1.0);
    Side bottom(1, 0.0);
    Side top(1, 1.0);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
    bottom.mark(*meshfuncedgeids, 3);
    top.mark(*meshfuncedgeids, 4);
  }
  else if (source=="Rectangle")
  {
    std::vector<double> lowerleft;
    buffer.str(""); buffer << optionpath << "/source/lower_left";
    serr = Spud::get_option(buffer.str(), lowerleft); spud_err(buffer.str(), serr);
    
    std::vector<double> upperright;
    buffer.str(""); buffer << optionpath << "/source/upper_right";
    serr = Spud::get_option(buffer.str(), upperright); spud_err(buffer.str(), serr);
    
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::Rectangle(lowerleft[0], lowerleft[1], 
                                      upperright[0], upperright[1], 
                                      cells[0], cells[1], diagonal) );

    MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains", 1);
    Side left(0, lowerleft[0]);
    Side right(0, upperright[0]);
    Side bottom(1, lowerleft[1]);
    Side top(1, upperright[1]);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
    bottom.mark(*meshfuncedgeids, 3);
    top.mark(*meshfuncedgeids, 4);
  }
  else if (source=="UnitCircle")
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); spud_err(buffer.str(), serr);
    
    std::string transformation;
    buffer.str(""); buffer << optionpath << "/source/transformation";
    serr = Spud::get_option(buffer.str(), transformation); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitCircle(cells, diagonal, transformation) );

  }
  else if (source=="UnitCube")
  {
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitCube(cells[0], cells[1], cells[2]) );

    MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains", 2);
    Side left(0, 0.0);
    Side right(0, 1.0);
    Side bottom(2, 0.0);
    Side top(2, 1.0);
    Side back(1, 0.0);
    Side front(1, 1.0);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
    bottom.mark(*meshfuncedgeids, 3);
    top.mark(*meshfuncedgeids, 4);
    back.mark(*meshfuncedgeids, 5);
    front.mark(*meshfuncedgeids, 6);
  }
  else if (source=="Box")
  {
    std::vector<double> lowerbackleft;
    buffer.str(""); buffer << optionpath << "/source/lower_back_left";
    serr = Spud::get_option(buffer.str(), lowerbackleft); spud_err(buffer.str(), serr);
    
    std::vector<double> upperfrontright;
    buffer.str(""); buffer << optionpath << "/source/upper_front_right";
    serr = Spud::get_option(buffer.str(), upperfrontright); spud_err(buffer.str(), serr);
    
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::Box(lowerbackleft[0], lowerbackleft[1], lowerbackleft[2],
                                upperfrontright[0], upperfrontright[1], upperfrontright[2], 
                                cells[0], cells[1], cells[2]) );

    MeshFunction_uint_ptr meshfuncedgeids = (*mesh).data().create_mesh_function("exterior_facet_domains", 2);
    Side left(0, lowerbackleft[0]);
    Side right(0, upperfrontright[0]);
    Side bottom(2, lowerbackleft[2]);
    Side top(2, upperfrontright[2]);
    Side back(1, lowerbackleft[1]);
    Side front(1, upperfrontright[1]);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
    bottom.mark(*meshfuncedgeids, 3);
    top.mark(*meshfuncedgeids, 4);
    back.mark(*meshfuncedgeids, 5);
    front.mark(*meshfuncedgeids, 6);
  }
  else if (source=="UnitSphere")
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitSphere(cells) );

  }
  else
  {
    dolfin::error("Unknown mesh source.");
  }

  // Put the mesh into the bucket
  register_mesh(mesh, meshname, optionpath);
}

void SpudBucket::uflnames_fill_(const std::string &optionpath)
{
  std::stringstream buffer;
  Spud::OptionError serr;
  std::string uflsymbol, type;

  buffer.str("");  buffer << optionpath << "/coefficient";
  int ncoeffs = Spud::option_count(buffer.str());
  // Loop over the coefficients and register their ufsymbols
  for (uint i = 0; i < ncoeffs; i++)
  {
    buffer.str(""); buffer << optionpath << "/coefficient[" << i << "]/type/name";
    serr = Spud::get_option(buffer.str(), type); spud_err(buffer.str(), serr);
    if (type=="Function")
    {
      buffer.str(""); buffer << optionpath << "/coefficient[" << i << "]/ufl_symbol";
      serr = Spud::get_option(buffer.str(), uflsymbol); spud_err(buffer.str(), serr);
      register_uflname(uflsymbol, uflsymbol);
      register_uflname(uflsymbol, uflsymbol+"_i");
      register_uflname(uflsymbol, uflsymbol+"_n");
    }
  }
}

// Insert a system (with optionpath) into the bucket
void SpudBucket::systems_fill_(const std::string &optionpath)
{
  SpudSystemBucket_ptr system( new SpudSystemBucket(optionpath, this) );

  (*system).fill();

  register_system(system, (*system).name());
}

// Register a pointer to a DOLFIN mesh (and an optionpath)
void SpudBucket::register_mesh(Mesh_ptr mesh, std::string name, std::string optionpath)
{
  // First check if a mesh with this name already exists
  Mesh_it m_it = meshes_.find(name);
  if (m_it != meshes_end())
  {
    // if it does, issue an error
    dolfin::error("Mesh named \"%s\" already exists in spudbucket.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    meshes_[name]           = mesh;
    mesh_optionpaths_[name] = optionpath;
  }
}

// Fetch the option path for a mesh
std::string SpudBucket::fetch_mesh_optionpath(const std::string name)
{
  // Check if the mesh exists in the bucket
  string_it s_it = mesh_optionpaths_.find(name);
  if (s_it == mesh_optionpaths_end())
  {
    // if it doesn't then throw an error
    dolfin::error("Mesh named \"%s\" does not exist in spudbucket.", name.c_str());
  }
  else
  {
    // if it does then return the pointer to the mesh
    return (*s_it).second;
  }
}

// Public iterator access
string_it SpudBucket::mesh_optionpaths_begin()
{
  return mesh_optionpaths_.begin();
}

// Public iterator access
string_const_it SpudBucket::mesh_optionpaths_begin() const
{
  return mesh_optionpaths_.begin();
}

// Public iterator access
string_it SpudBucket::mesh_optionpaths_end()
{
  return mesh_optionpaths_.end();
}

// Public iterator access
string_const_it SpudBucket::mesh_optionpaths_end() const
{
  return mesh_optionpaths_.end();
}

// Describe the contents of the mesh_optionpaths_ map
const std::string SpudBucket::meshes_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = mesh_optionpaths_begin(); s_it != mesh_optionpaths_end(); s_it++ )
  {
    s << indentation << "Mesh " << (*s_it).first << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

void SpudBucket::run()
{
  Spud::OptionError serr;
  std::stringstream buffer;

  buffer.str(""); buffer << "/timestepping";
  if (Spud::have_option(buffer.str()))
  {
    timestep_run_(); 
  }
  else
  {
    steady_run_();
  }

}

void SpudBucket::timestep_run_()
{
  Spud::OptionError serr;
  std::stringstream buffer;

  double time;
  buffer.str(""); buffer << "/timestepping/current_time";
  serr = Spud::get_option(buffer.str(), time); spud_err(buffer.str(), serr);

  double finish_time;
  buffer.str(""); buffer << "/timestepping/finish_time";
  serr = Spud::get_option(buffer.str(), finish_time); spud_err(buffer.str(), serr);

  Constant_ptr timestep;
  std::string systemname;
  buffer.str(""); buffer << "/timestepping/timestep/system/name";
  serr = Spud::get_option(buffer.str(), systemname); spud_err(buffer.str(), serr);
  std::string coeffname;
  buffer.str(""); buffer << "/timestepping/timestep/coefficient/name";
  serr = Spud::get_option(buffer.str(), coeffname); spud_err(buffer.str(), serr);
  timestep = boost::dynamic_pointer_cast< dolfin::Constant >((*fetch_system(systemname)).fetch_coeff(coeffname));

  int nints;
  buffer.str(""); buffer << "/nonlinear_systems/nonlinear_iterations";
  serr = Spud::get_option(buffer.str(), nints, 1); spud_err(buffer.str(), serr);

  int timestep_count = 0;

  while(time < finish_time)
  {
    dolfin::error("Timestepping not implemented.");
    for (uint i = 0; i < nints; i++)
    {
      solve();
    }

    time += double(*timestep);
    timestep_count++;
  }

}

void SpudBucket::steady_run_()
{
  Spud::OptionError serr;
  std::stringstream buffer;

  int nints;
  buffer.str(""); buffer << "/nonlinear_systems/nonlinear_iterations";
  serr = Spud::get_option(buffer.str(), nints, 1); spud_err(buffer.str(), serr);

  for (uint i = 0; i < nints; i++)
  {

    solve();

  }

}

// Empty the spudbucket
void SpudBucket::empty_()
{
  mesh_optionpaths_.clear();
}

//void SpudBucket::detectors_fill_()
//{
//  // Initialise a pointer to a generic detector - so that it can be reset in the loops
//  GenericDetectors_ptr det(new GenericDetectors());
//  // A string buffer for option paths
//  std::stringstream buffer;
//  
//  // Find out how many point detectors there are
//  int npdets = Spud::option_count("/io/detectors/point");
//  // and loop over them
//  for (uint i=0; i<npdets; i++)
//  {
//    // Set up the base path for a point detector
//    std::stringstream detectorpath;
//    detectorpath.str(""); detectorpath << "/io/detectors/point[" << i << "]";
//    
//    // Get the name of the detector
//    std::string detname;
//    buffer.str(detectorpath.str()); buffer << "/name";
//    Spud::get_option(buffer.str(), detname);
//
//    // Get the location of the detector
//    std::vector<double> point;
//    buffer.str(detectorpath.str());
//    Spud::get_option(buffer.str(), point);
//    
//    // Initialise and register the point detector
//    det.reset(new PointDetectors(&point, detname));
//    register_detector(det, detname);
//  }
//  
//  // Find out how many detector arrays there are
//  int nadets = Spud::option_count("/io/detectors/array");
//  // and loop over them
//  for (uint i=0; i<nadets; i++)
//  {
//    // Set up the base path for a detector array
//    std::stringstream detectorpath;
//    detectorpath.str(""); detectorpath << "/io/detectors/array[" << i << "]";
//    
//    // Get the number of points in this array
//    int no_det;
//    buffer.str(detectorpath.str()); buffer << "/number_of_detectors";
//    Spud::get_option(buffer.str(), no_det);
//
//    // Get the detector array name
//    std::string detname;
//    buffer.str(detectorpath.str()); buffer << "/name";
//    Spud::get_option(buffer.str(), detname);
//
//    // Get the python function that describes the positions of the detectors
//    std::string function;
//    buffer.str(detectorpath.str()); buffer << "/python";
//    Spud::get_option(buffer.str(), function);
//    
//    // Get the dimension of this problem
//    // - here we assume that there is only one dimension in this problem
//    int dimension;
//    Spud::get_option("/geometry/dimension", dimension);
//
//    // Create the detector and register it in the bucket
//    det.reset(new PythonDetectors(no_det, dimension, function, detname));
//    register_detector(det, detname);
//  }  
//  
//}
