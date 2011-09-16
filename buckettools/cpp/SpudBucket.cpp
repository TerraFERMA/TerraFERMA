
#include "SpudBucket.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "BoostTypes.h"
#include "BucketDolfinBase.h"
//#include "GenericDetectors.h"
#include "PointDetectors.h"
#include "PythonDetectors.h"
//#include "PythonExpression.h"
//#include "InitialConditionExpression.h"
#include "DiagnosticsFile.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
SpudBucket::SpudBucket() : Bucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudBucket::SpudBucket(const std::string &name) : Bucket(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudBucket::SpudBucket(const std::string &name, const std::string &optionpath) : 
                                optionpath_(optionpath), Bucket(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudBucket::~SpudBucket()
{
  empty_();                                                          // empty the data structures
}

//*******************************************************************|************************************************************//
// run the model described by this bucket
//*******************************************************************|************************************************************//
void SpudBucket::run()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  std::string basename;                                              // get the output base name
  buffer.str(""); buffer << "/io/output_base_name";
  serr = Spud::get_option(buffer.str(), basename); 
  spud_err(buffer.str(), serr);

  DiagnosticsFile diagfile(basename+".stat");
  diagfile.write_header(*this);
  diagfile.write_data(*this);

  output();

  do {                                                               // loop over time

    std::cout << "Time: " << current_time() << std::endl << std::flush;
    std::cout << "Timestep number: " << timestep_count() << std::endl << std::flush;

    for (*iteration_count_ = 0; \
         *iteration_count_ < nonlinear_iterations(); 
         (*iteration_count_)++)                                      // loop over the nonlinear iterations
    {
      solve();                                                       // solve all systems in the bucket
    }

    update();                                                        // update all functions in the bucket
    *current_time_ += timestep();                                    // increment time with the timestep
    (*timestep_count_)++;                                            // increment the number of timesteps taken

    diagfile.write_data(*this);
    output();

  } while (*current_time_ < finish_time());                          // syntax ensures at least one solve

  diagfile.close();                                                  // close the diagnostics

}

//*******************************************************************|************************************************************//
// fill the bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudBucket::fill()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  buffer.str(""); buffer << optionpath() << "/geometry/dimension";   // geometry dimension set in the bucket to pass it down to all
  serr = Spud::get_option(buffer.str(), dimension_);                 // systems (we assume this is the length of things that do
  spud_err(buffer.str(), serr);                                      // not have them independently specified)

  timestepping_fill_();                                              // fill in the timestepping options (if there are any)

  buffer.str(""); buffer << optionpath() << "/geometry/mesh";        // put the meshes into the bucket
  int nmeshes = Spud::option_count(buffer.str());
  for (uint i = 0; i<nmeshes; i++)                                   // loop over the meshes defined in the options file
  {
    buffer.str(""); buffer << optionpath() << "/geometry/mesh[" 
                                                        << i << "]";
    meshes_fill_(buffer.str());
  }
  
  buffer.str(""); buffer << optionpath() << "/system";
  int nsystems = Spud::option_count(buffer.str());
  for (uint i = 0; i<nsystems; i++)                                  // loop over the systems registering the base uflsymbols of any
  {                                                                  // coefficient functions contained in them
    buffer.str(""); buffer << optionpath() << "/system[" << i << "]";
    baseuflsymbols_fill_(buffer.str());
  }

  for (uint i = 0; i<nsystems; i++)                                  // loop over the systems *again*, this time filling all the
  {                                                                  // systems data into the bucket data structures
    buffer.str(""); buffer << optionpath() << "/system[" << i << "]";
    systems_fill_(buffer.str());
  }

  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *third* time, this time filling
                                  sys_it != systems_end(); sys_it++) // in the data for the coefficient functions contained within
  {                                                                  // them
    (*boost::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).funccoeffs_fill();
  }                                                                  // we couldn't do this before because we might not have had
                                                                     // the right functionspace available
  
  uflsymbols_fill_();                                                // now all the functions in the systems are complete we can 
                                                                     // register them in the bucket so it's easy to attach them
                                                                     // to the forms and functionals

  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *fourth* time, attaching the
                                  sys_it != systems_end(); sys_it++) // coefficients to the forms and functionals and initializing
  {                                                                  // the matrices by performing a preassembly step on them
    (*(*sys_it).second).attach_and_initialize();
  }
  
  detectors_fill_();                                                 // put the detectors in the bucket

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a dolfin mesh in the bucket (and spudbucket) data maps with a spud optionpath
//*******************************************************************|************************************************************//
void SpudBucket::register_mesh(Mesh_ptr mesh, 
                               const std::string &name, 
                               const std::string &optionpath)
{
  Mesh_it m_it = meshes_.find(name);                                 // check if a mesh with this name already exists
  if (m_it != meshes_end())
  {
    dolfin::error("Mesh named \"%s\" already exists in spudbucket.", // if it does, issue an error
                                                      name.c_str());
  }
  else
  {
    meshes_[name]           = mesh;                                  // if not register it in the map
    mesh_optionpaths_[name] = optionpath;                            // also register its optionpath
  }
}

//*******************************************************************|************************************************************//
// return a string containing the mesh optionpath in the spudbucket data maps
//*******************************************************************|************************************************************//
std::string SpudBucket::fetch_mesh_optionpath(const std::string &name)
{
  string_it s_it = mesh_optionpaths_.find(name);                     // check if a mesh with this name exists
  if (s_it == mesh_optionpaths_end())
  {
    dolfin::error("Mesh named \"%s\" does not exist in spudbucket.", // if it doesn't, issue an error
                                                      name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return the optionpath
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::mesh_optionpaths_begin()
{
  return mesh_optionpaths_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::mesh_optionpaths_begin() const
{
  return mesh_optionpaths_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::mesh_optionpaths_end()
{
  return mesh_optionpaths_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::mesh_optionpaths_end() const
{
  return mesh_optionpaths_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a detector in the bucket (and spudbucket) data maps with a spud optionpath
//*******************************************************************|************************************************************//
void SpudBucket::register_detector(GenericDetectors_ptr detector, 
                               const std::string &name, 
                               const std::string &optionpath)
{
  GenericDetectors_it d_it = detectors_.find(name);                  // check if a detector set with this name already exists
  if (d_it != detectors_end())
  {
    dolfin::error(
              "Detector named \"%s\" already exists in spudbucket.", // if it does, issue an error
                                                      name.c_str());
  }
  else
  {
    detectors_[name]            = detector;                          // if not register it in the map
    detector_optionpaths_[name] = optionpath;                        // also register its optionpath
  }
}

//*******************************************************************|************************************************************//
// return a string containing the detector optionpath in the spudbucket data maps
//*******************************************************************|************************************************************//
std::string SpudBucket::fetch_detector_optionpath(const std::string &name)
{
  string_it s_it = detector_optionpaths_.find(name);                 // check if a mesh with this name exists
  if (s_it == detector_optionpaths_end())
  {
    dolfin::error(
              "Detector named \"%s\" does not exist in spudbucket.",
                                                      name.c_str()); // if it doesn't, issue an error
  }
  else
  {
    return (*s_it).second;                                           // if it does, return the optionpath
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the detector_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::detector_optionpaths_begin()
{
  return detector_optionpaths_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the detector_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::detector_optionpaths_begin() const
{
  return detector_optionpaths_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the detector_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::detector_optionpaths_end()
{
  return detector_optionpaths_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::detector_optionpaths_end() const
{
  return detector_optionpaths_.end();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the mesh_optionpaths_ data structure
//*******************************************************************|************************************************************//
const std::string SpudBucket::meshes_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = mesh_optionpaths_begin(); 
                            s_it != mesh_optionpaths_end(); s_it++ )
  {
    s << indentation << "Mesh " << (*s_it).first 
                    << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the detector_optionpaths_ data structure
//*******************************************************************|************************************************************//
const std::string SpudBucket::detectors_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = detector_optionpaths_begin(); 
                            s_it != detector_optionpaths_end(); s_it++ )
  {
    s << indentation << "Detector set " << (*s_it).first 
                    << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// fill in any timestepping data or set up dummy values instead (zero essentially)
//*******************************************************************|************************************************************//
void SpudBucket::timestepping_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error
  
  timestep_count_.reset( new int );
  *timestep_count_ = 0;                                              // the number of timesteps taken
  iteration_count_.reset( new int );
  *iteration_count_ = 0;                                             // the number of iterations taken

  current_time_.reset( new double );
  buffer.str(""); buffer << "/timestepping/current_time";            // get the current time
  serr = Spud::get_option(buffer.str(), *current_time_, 0.0);        // may be non-zero at the start of simulations (e.g.
  spud_err(buffer.str(), serr);                                      // checkpoints) but assumed zero for steady simulations

  finish_time_.reset( new double );
  buffer.str(""); buffer << "/timestepping/finish_time";             // get the finish time (assumed zero for steady simulations)
  serr = Spud::get_option(buffer.str(), *finish_time_, 0.0); 
  spud_err(buffer.str(), serr);

  nonlinear_iterations_.reset( new int );
  buffer.str(""); buffer                                             // find out if we're doing nonlinear iterations
                    << "/nonlinear_systems/nonlinear_iterations";    // between the systems
  serr = Spud::get_option(buffer.str(), *nonlinear_iterations_, 1); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << "/timestepping";
  if (Spud::have_option(buffer.str()))
  {
    buffer.str(""); 
    buffer << "/timestepping/timestep/coefficient::Timestep/ufl_symbol";
    serr = Spud::get_option(buffer.str(), timestep_.first);
    spud_err(buffer.str(), serr);

    buffer.str(""); 
    buffer << "/timestepping/timestep/coefficient::Timestep/type::Constant/rank::Scalar/value::WholeMesh";
    timestep_.second = boost::dynamic_pointer_cast< dolfin::Constant >(initialize_expression(buffer.str()));
  }
  else
  {
    timestep_.first = "";
    timestep_.second.reset( new dolfin::Constant(0.0) );
  }
  
}

//*******************************************************************|************************************************************//
// create or get from a file a dolfin mesh object and insert it into the bucket data structures
//*******************************************************************|************************************************************//
void SpudBucket::meshes_fill_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  std::string meshname;                                              // get the name of the mesh
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), meshname); 
  spud_err(buffer.str(), serr);
  
  std::string source;                                                // get the source of the mesh (i.e. from file or internal)
  buffer.str(""); buffer << optionpath << "/source/name";
  serr = Spud::get_option(buffer.str(), source); 
  spud_err(buffer.str(), serr);

  Mesh_ptr mesh;                                                     // initialize the pointer

  if (source=="File")                                                // source is a file
  {
    std::string basename;                                            // get the base file name (without the .xml)
    buffer.str(""); buffer << optionpath << "/source/file";
    serr = Spud::get_option(buffer.str(), basename); 
    spud_err(buffer.str(), serr);

    std::stringstream filename;                                      // read in the mesh (using the dolfin constructor)
    filename.str(""); filename << basename << ".xml";
    mesh.reset(new dolfin::Mesh(filename.str()));
    (*mesh).init();                                                  // initialize the mesh (maps between dimensions etc.)

    std::ifstream file;                                              // a dummy file stream to test if files exist
                                                                     // (better way of doing this?)
    
    filename.str(""); filename << basename << "_edge_subdomain.xml"; // check if the edge subdomain mesh function file exists
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)                                                        // if it does then attach it to the dolfin MeshData structure 
    {                                                                // using the dolfin reserved name for exterior facets
      file.close();
      MeshFunction_uint_ptr meshfuncedgeids = 
        (*mesh).data().create_mesh_function("exterior_facet_domains");
      dolfin::MeshFunction<dolfin::uint> edgeids(*mesh, filename.str());
      *meshfuncedgeids = edgeids;
    }

    filename.str(""); filename << basename << "_cell_subdomain.xml"; // check if the edge subdomain mesh function file exists
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)                                                        // if it does then attach it to the dolfin MeshData structure 
    {                                                                // using the dolfin reserved name for cell domains
      file.close();
      MeshFunction_uint_ptr meshfunccellids = 
        (*mesh).data().create_mesh_function("cell_domains");
      dolfin::MeshFunction<dolfin::uint> cellids(*mesh, filename.str());
      *meshfunccellids = cellids;
    }

  }
  else if (source=="UnitInterval")                                   // source is an internally generated dolfin mesh
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitInterval(cells) );

    MeshFunction_uint_ptr meshfuncedgeids = 
      (*mesh).data().create_mesh_function("exterior_facet_domains", 0);
    Side left(0, 0.0);
    Side right(0, 1.0);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
  }
  else if (source=="Interval")                                       // source is an internally generated dolfin mesh
  {
    double leftx;
    buffer.str(""); buffer << optionpath << "/source/left";
    serr = Spud::get_option(buffer.str(), leftx); 
    spud_err(buffer.str(), serr);

    double rightx;
    buffer.str(""); buffer << optionpath << "/source/right";
    serr = Spud::get_option(buffer.str(), rightx); 
    spud_err(buffer.str(), serr);

    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::Interval(cells, leftx, rightx) );

    MeshFunction_uint_ptr meshfuncedgeids = 
      (*mesh).data().create_mesh_function("exterior_facet_domains", 0);
    Side left(0, leftx);
    Side right(0, rightx);
    (*meshfuncedgeids).set_all(0);
    left.mark(*meshfuncedgeids, 1);
    right.mark(*meshfuncedgeids, 2);
  }
  else if (source=="UnitSquare")                                     // source is an internally generated dolfin mesh
  {
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitSquare(cells[0], cells[1], diagonal) );

    MeshFunction_uint_ptr meshfuncedgeids = 
      (*mesh).data().create_mesh_function("exterior_facet_domains", 1);
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
  else if (source=="Rectangle")                                      // source is an internally generated dolfin mesh
  {
    std::vector<double> lowerleft;
    buffer.str(""); buffer << optionpath << "/source/lower_left";
    serr = Spud::get_option(buffer.str(), lowerleft); 
    spud_err(buffer.str(), serr);
    
    std::vector<double> upperright;
    buffer.str(""); buffer << optionpath << "/source/upper_right";
    serr = Spud::get_option(buffer.str(), upperright); 
    spud_err(buffer.str(), serr);
    
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); 
    spud_err(buffer.str(), serr);
    
    mesh.reset(new dolfin::Rectangle(lowerleft[0], lowerleft[1], 
                                    upperright[0], upperright[1], 
                                    cells[0], cells[1], diagonal));

    MeshFunction_uint_ptr meshfuncedgeids = 
      (*mesh).data().create_mesh_function("exterior_facet_domains", 1);
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
  else if (source=="UnitCircle")                                     // source is an internally generated dolfin mesh
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); 
    spud_err(buffer.str(), serr);
    
    std::string transformation;
    buffer.str(""); buffer << optionpath << "/source/transformation";
    serr = Spud::get_option(buffer.str(), transformation); 
    spud_err(buffer.str(), serr);
    
    mesh.reset(new dolfin::UnitCircle(cells, 
                                          diagonal, transformation));

  }
  else if (source=="UnitCube")                                       // source is an internally generated dolfin mesh
  {
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitCube(cells[0], 
                                     cells[1], 
                                     cells[2]) );

    MeshFunction_uint_ptr meshfuncedgeids = 
      (*mesh).data().create_mesh_function("exterior_facet_domains", 2);
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
  else if (source=="Box")                                            // source is an internally generated dolfin mesh
  {
    std::vector<double> lowerbackleft;
    buffer.str(""); buffer << optionpath 
                              << "/source/lower_back_left";
    serr = Spud::get_option(buffer.str(), lowerbackleft); 
    spud_err(buffer.str(), serr);
    
    std::vector<double> upperfrontright;
    buffer.str(""); buffer << optionpath 
                            << "/source/upper_front_right";
    serr = Spud::get_option(buffer.str(), upperfrontright); 
    spud_err(buffer.str(), serr);
    
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::Box(lowerbackleft[0], 
                                lowerbackleft[1], 
                                lowerbackleft[2],
                                upperfrontright[0], 
                                upperfrontright[1], 
                                upperfrontright[2], 
                                cells[0], cells[1], cells[2]) );

    MeshFunction_uint_ptr meshfuncedgeids = 
      (*mesh).data().create_mesh_function("exterior_facet_domains", 2);
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
  else if (source=="UnitSphere")                                     // source is an internally generated dolfin mesh
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitSphere(cells) );

  }
  else                                                               // source is unrecognised
  {
    dolfin::error("Unknown mesh source.");
  }

  register_mesh(mesh, meshname, optionpath);                         // put the new mesh in the bucket
}

//*******************************************************************|************************************************************//
// create a new system, fill it and put it into the bucket
//*******************************************************************|************************************************************//
void SpudBucket::systems_fill_(const std::string &optionpath)
{
  SpudSystemBucket_ptr system(new SpudSystemBucket(optionpath,       // create a new system (assumed to be a spudsystem with this 
                                                            this));  // bucket as a parent)

  (*system).fill();                                                  // fill the system

  register_system(system, (*system).name());                         // put the system in the bucket
}

//*******************************************************************|************************************************************//
// loop over the systems defined in the options dictionary and register the base uflsymbols for all coefficient functions
//*******************************************************************|************************************************************//
void SpudBucket::baseuflsymbols_fill_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::string uflsymbol, type;                                       // base ufl symbol and coefficient type

  buffer.str("");  buffer << optionpath << "/coefficient";           // the optionpath is for a system so append the coefficient tag
  int ncoeffs = Spud::option_count(buffer.str());                    // and find out how many there are
  for (uint i = 0; i < ncoeffs; i++)                                 // loop over them
  {
    buffer.str(""); buffer << optionpath << "/coefficient[" << i 
                                                    << "]/type/name";
    serr = Spud::get_option(buffer.str(), type);                     // find out the coefficient type
    spud_err(buffer.str(), serr);
    if (type=="Function")                                            // if it's a function then register its derived ufl symbols
    {
      buffer.str(""); buffer << optionpath << "/coefficient[" << i 
                                                   << "]/ufl_symbol";
      serr = Spud::get_option(buffer.str(), uflsymbol); 
      spud_err(buffer.str(), serr);
      register_baseuflsymbol(uflsymbol, uflsymbol);
      register_baseuflsymbol(uflsymbol, uflsymbol+"_i");
      register_baseuflsymbol(uflsymbol, uflsymbol+"_n");
    }
  }
}

//*******************************************************************|************************************************************//
// loop over the detectors defined in the options dictionary and set up the requested detectors
//*******************************************************************|************************************************************//
void SpudBucket::detectors_fill_()
{
  GenericDetectors_ptr det(new GenericDetectors());                  // initialize a pointer to a generic detector
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  int npdets = Spud::option_count("/io/detectors/point");            // number of point detectors
  for (uint i=0; i<npdets; i++)                                      // loop over point detectors
  {
    std::stringstream detectorpath;
    detectorpath.str(""); detectorpath << "/io/detectors/point[" 
                                                        << i << "]";
    
    std::string detname;                                             // detector name
    buffer.str(""); buffer << detectorpath.str() << "/name";
    serr = Spud::get_option(buffer.str(), detname);
    spud_err(buffer.str(), serr);

    std::vector<double> point;                                       // point location
    buffer.str(""); buffer << detectorpath.str();
    serr = Spud::get_option(buffer.str(), point);
    spud_err(buffer.str(), serr);
    
    det.reset(new PointDetectors(point, detname));                   // create point detector
    register_detector(det, detname, detectorpath.str());             // register detector
  }
  
  int nadets = Spud::option_count("/io/detectors/array");            // number of array detectors
  for (uint i=0; i<nadets; i++)                                      // loop over array detectors
  {
    std::stringstream detectorpath;
    detectorpath.str(""); detectorpath << "/io/detectors/array[" 
                                                        << i << "]";
    
    int no_det;                                                      // number of detectors in array
    buffer.str(""); buffer << detectorpath.str() << "/number_of_detectors";
    serr = Spud::get_option(buffer.str(), no_det);
    spud_err(buffer.str(), serr);

    std::string detname;                                             // detector array name
    buffer.str(""); buffer << detectorpath.str() << "/name";
    serr = Spud::get_option(buffer.str(), detname);
    spud_err(buffer.str(), serr);

    std::string function;                                            // python function describing the detector positions
    buffer.str(""); buffer << detectorpath.str() << "/python";
    serr = Spud::get_option(buffer.str(), function);
    spud_err(buffer.str(), serr);
    
                                                                     // create python detectors array
    det.reset(new PythonDetectors(no_det, dimension(), function, detname));
    register_detector(det, detname, detectorpath.str());             // register detector
  }  
  
}

//*******************************************************************|************************************************************//
// empty the data structures in the spudbucket
//*******************************************************************|************************************************************//
void SpudBucket::empty_()
{
  mesh_optionpaths_.clear();
  Bucket::empty_();
}

