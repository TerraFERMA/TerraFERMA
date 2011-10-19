
#include "Bucket.h"
// #include "GenericDetectors.h"
// #include "PythonDetectors.h"
#include <dolfin.h>
#include <string>
#include "SignalHandler.h"
#include "EventHandler.h"
#include "StatisticsFile.h"
#include <signal.h>
#include <time.h>

using namespace buckettools;

time_t Bucket::start_walltime_ = time(NULL);                         // initialize global static variable
boost::timer Bucket::timer_;                                         // start timing the simulation

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
Bucket::Bucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
Bucket::Bucket(const std::string &name) : name_(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
Bucket::~Bucket()
{
  empty_();                                                          // empty the data structures (unnecessary?)
}

//*******************************************************************|************************************************************//
// run the model described by this bucket
//*******************************************************************|************************************************************//
void Bucket::run()
{
  std::stringstream buffer;                                          // optionpath buffer

  update();

  output(OUTPUT_START);

  solve_at_start_();

  dolfin::log(dolfin::INFO, "Entering timeloop.");
  bool continue_timestepping = true;
  while (continue_timestepping) 
  {                                                                  // loop over time

    *old_time_ = *current_time_;                                     // old time is now the previous time??
    *current_time_ += timestep();                                    // increment time with the timestep 
                                                                     // (we do this now so that time dependent expressions are
                                                                     // evaluating at the right time, i.e. symbol_i, which is
                                                                     // attached to current_time, is at the next
                                                                     // time level, and symbol_n, which is attached to old_time
                                                                     // is at the previous time level throughout the whole timestep)
    (*timestep_count_)++;                                            // increment the number of timesteps taken

    dolfin::log(dolfin::INFO, "Timestep numbers: %d -> %d", timestep_count()-1, timestep_count());
    dolfin::log(dolfin::INFO, "Times: %f -> %f", old_time(), current_time());
    dolfin::log(dolfin::INFO, "Timestep: %f", timestep());

    solve_in_timeloop_();                                            // this is where the magic happens

    if(complete())                                                   // this must be called before the update as it checks if a
    {                                                                // steady state has been attained
      output(OUTPUT_END);                                            // force an output at the end
      continue_timestepping = false;                                 // signal to stop timestepping
    }
    else
    {
      output(OUTPUT_TIMELOOP);                                       // standard timeloop output
    }

    update_timestep();                                               // update the timestep

    update();                                                        // update all functions in the bucket

  }                                                                  // syntax ensures at least one solve
  dolfin::log(dolfin::INFO, "Finished timeloop.");

}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling solve on each of them
//*******************************************************************|************************************************************//
void Bucket::solve(const int &location)
{
  for (int_SystemBucket_const_it s_it = orderedsystems_begin(); 
                              s_it != orderedsystems_end(); s_it++)
  {
    if((*(*s_it).second).solve_location()==location)
    {
      (*(*s_it).second).solve();
    }
  }
}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling update on each of them
//*******************************************************************|************************************************************//
void Bucket::update()
{
  for (int_SystemBucket_const_it s_it = orderedsystems_begin(); 
                                 s_it != orderedsystems_end(); s_it++)
  {
    (*(*s_it).second).update();
  }
}

//*******************************************************************|************************************************************//
// loop over the constraints, solving the listed systems and updating the timestep based on them
//*******************************************************************|************************************************************//
void Bucket::update_timestep()
{
  if (timestep_constraints_.size()==0)                               // this indicates we don't have adaptive timestepping turned on
  {
    return;                                                          // so return
  }

  bool zero_init_dt = ( (std::abs(timestep()) < DOLFIN_EPS)          // if the timestep is zero and
                                   && (timestep_count()==1) );       // this is the end of the first timestep
                                                                     // then we consider this a zero initial timestep simulation
                                                                     // and force an adapt (and solve) regardless of the period

  bool adapt_dt = perform_action_(timestepadapt_period_, 
                                  timestepadapt_time_, 
                                  timestepadapt_period_timesteps_);  // do we want to update the timestep now

  if ((!adapt_dt) && (!zero_init_dt))
  {
    return;
  }

  double new_dt = HUGE_VAL;
  dolfin::log(dolfin::INFO, "In update_timestep()");

  if (zero_init_dt)                                                  // the timestep is zero initially so we'll get a 0.0 Courant
  {                                                                  // like number... so set a dummy timestep of 1.0 for this
    *(timestep_.second) = 1.0;                                       // calculation
    dolfin::log(dolfin::INFO, "Forcing a timestep adapt as initial timestep is 0.0.");
  }

                                 //  system name, field name   , requested maxval
  std::vector< std::pair< std::pair< std::string, std::string >, double > >::const_iterator c_it;
  for (c_it =  timestep_constraints_.begin(); 
       c_it != timestep_constraints_.end(); 
       c_it++)
  {
    double suggested_dt;

    SystemBucket_ptr sysbucket = fetch_system((*c_it).first.first);
    if((!(*sysbucket).solved()) || zero_init_dt )                    // force a solve of the requested system if it hasn't already
    {                                                                // been solved for this timestep or we've just set a dummy
      (*sysbucket).solve();                                          // timestep of 1.0, which will have to be incorporated into the
    }                                                                // new calculation

    FunctionBucket_ptr funcbucket = (*sysbucket).fetch_field((*c_it).first.second);
    dolfin::Function func =                                          // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >((*funcbucket).function());

    double maxval = func.vector().max();                             // work out the current maximum value

    if (maxval==0.0)
    {
      suggested_dt = HUGE_VAL;
    }
    else
    { 
      suggested_dt = ((*c_it).second*timestep())/maxval;
    }
    new_dt = std::min(suggested_dt, new_dt);
  }

  if (timestep_increasetol_)
  {
    new_dt = std::min(timestep()*(*timestep_increasetol_), new_dt);
  }
 
  *(timestep_.second) = new_dt;

}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling update_nonlinear on each of them
//*******************************************************************|************************************************************//
void Bucket::update_nonlinear()
{
  for (int_SystemBucket_const_it s_it = orderedsystems_begin(); 
                             s_it != orderedsystems_end(); s_it++)
  {
    (*(*s_it).second).update_nonlinear();
  }
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the simulation has finished or not (normally for reasons other than the time being complete)
//*******************************************************************|************************************************************//
bool Bucket::complete()
{
  bool completed = false;
  
  if (current_time() >= finish_time())
  {
    dolfin::log(dolfin::WARNING, "Finish time reached, terminating timeloop.");
    completed = true;
  }

  if (steadystate_())
  {
    dolfin::log(dolfin::WARNING, "Steady state attained, terminating timeloop.");
    completed = true;
  }

  if ((*(*SignalHandler::instance()).return_handler(SIGINT)).received())
  {
    dolfin::log(dolfin::ERROR, "SigInt received, terminating timeloop.");
    completed = true;
  }

  return completed;
}

//*******************************************************************|************************************************************//
// loop over the selected forms and attach the coefficients they request using the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::attach_coeffs(Form_it f_begin, Form_it f_end)
{
  for (Form_it f_it = f_begin; f_it != f_end; f_it++)                // loop over the forms
  {
    attach_coeffs((*f_it).second);
  }
}

//*******************************************************************|************************************************************//
// loop over a form and attach the coefficients they request using the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::attach_coeffs(Form_ptr form)
{
  uint ncoeff = (*form).num_coefficients();                          // find out how many coefficients this form requires
  for (uint i = 0; i < ncoeff; i++)                                  // loop over the required coefficients
  {
    std::string uflsymbol = (*form).coefficient_name(i);             // get the (possibly derived) ufl symbol of this coefficient
    GenericFunction_ptr function = fetch_uflsymbol(uflsymbol);       // grab the corresponding function (possibly old, iterated etc.
                                                                     // if ufl symbol is derived)

    (*form).set_coefficient(uflsymbol, function);                    // attach that function as a coefficient
  }
}

//*******************************************************************|************************************************************//
// make a partial copy of the provided bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void Bucket::copy_diagnostics(Bucket_ptr &bucket) const
{

  if(!bucket)
  {
    bucket.reset( new Bucket );
  }

  (*bucket).name_ = name_;

  (*bucket).start_time_ = start_time_;
  (*bucket).old_time_ = old_time_;
  (*bucket).current_time_ = current_time_;
  (*bucket).finish_time_ = finish_time_;
  (*bucket).timestep_count_ = timestep_count_;
  (*bucket).timestep_ = timestep_;
  (*bucket).nonlinear_iterations_ = nonlinear_iterations_;
  (*bucket).iteration_count_ = iteration_count_;

  (*bucket).meshes_ = meshes_;

  for (SystemBucket_const_it sys_it = systems_begin();               // loop over the systems
                           sys_it != systems_end(); sys_it++)
  {                                                                  
    SystemBucket_ptr system;                                         // create a new system
    
    (*(*sys_it).second).copy_diagnostics(system, bucket);

    (*bucket).register_system(system, (*system).name());             // put the system in the bucket
  }                                                                  

  (*bucket).detectors_ = detectors_;

  (*bucket).steadystate_tol_ = steadystate_tol_;

}

//*******************************************************************|************************************************************//
// return the timestep count
//*******************************************************************|************************************************************//
const int Bucket::timestep_count() const
{
  assert(timestep_count_);
  return *timestep_count_;
}

//*******************************************************************|************************************************************//
// return the start time
//*******************************************************************|************************************************************//
const double Bucket::start_time() const
{
  assert(start_time_);
  return *start_time_;
}

//*******************************************************************|************************************************************//
// return the old time (from the previous timestep)
//*******************************************************************|************************************************************//
const double Bucket::old_time() const
{
  assert(old_time_);
  return *old_time_;
}

//*******************************************************************|************************************************************//
// return the current time
//*******************************************************************|************************************************************//
const double Bucket::current_time() const
{
  assert(current_time_);
  return *current_time_;
}

//*******************************************************************|************************************************************//
// return the finish time
//*******************************************************************|************************************************************//
const double Bucket::finish_time() const
{
  assert(finish_time_);
  return *finish_time_;
}

//*******************************************************************|************************************************************//
// return the timestep (as a double)
//*******************************************************************|************************************************************//
const double Bucket::timestep() const
{
  assert(timestep_.second);
  return double(*(timestep_.second));
}

//*******************************************************************|************************************************************//
// return the number of nonlinear iterations requested
//*******************************************************************|************************************************************//
const int Bucket::nonlinear_iterations() const
{
  assert(nonlinear_iterations_);
  return *nonlinear_iterations_;
}

//*******************************************************************|************************************************************//
// return the number of nonlinear iterations taken
//*******************************************************************|************************************************************//
const int Bucket::iteration_count() const
{
  return *iteration_count_;
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a dolfin mesh in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_mesh(Mesh_ptr mesh, const std::string &name)
{
  Mesh_it m_it = meshes_.find(name);                                 // check if a mesh with this name already exists
  if (m_it != meshes_end())
  {
    dolfin::error("Mesh named \"%s\" already exists in bucket.",     // if it does, issue an error
                                                    name.c_str());
  }
  else
  {
    meshes_[name] = mesh;                                            // if not, insert it into the meshes_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a dolfin mesh in the bucket data maps
//*******************************************************************|************************************************************//
Mesh_ptr Bucket::fetch_mesh(const std::string &name)
{
  Mesh_it m_it = meshes_.find(name);                                 // check if this mesh exists in the meshes_ map
  if (m_it == meshes_end())
  {
    dolfin::error("Mesh named \"%s\" does not exist in bucket.",     // if it doesn't, issue an error
                                                    name.c_str());
  }
  else
  {
    return (*m_it).second;                                           // if it does, return a (boost shared) pointer to it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_it Bucket::meshes_begin()
{
  return meshes_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_const_it Bucket::meshes_begin() const
{
  return meshes_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_it Bucket::meshes_end()
{
  return meshes_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_const_it Bucket::meshes_end() const
{
  return meshes_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a system bucket in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_system(SystemBucket_ptr system, 
                                         const std::string &name)
{
  SystemBucket_it s_it = systems_.find(name);                        // check if a system with this name already exists
  if (s_it != systems_end())
  {
    dolfin::error(
            "SystemBucket named \"%s\" already exists in bucket",    // if it does, issue an error
                                  name.c_str());
  }
  else
  {
    systems_[name] = system;                                         // if not insert it into the systems_ map
    orderedsystems_[(int) systems_.size()] = system;                 // also insert it in the order it was registered in the 
                                                                     // orderedsystems_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a system bucket in the bucket data maps
//*******************************************************************|************************************************************//
SystemBucket_ptr Bucket::fetch_system(const std::string &name)
{
  SystemBucket_it s_it = systems_.find(name);                        // check if a system with this name exists in the bucket
  if (s_it == systems_end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "SystemBucket named \"%s\" does not exist in bucket.", 
                                name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return a pointer to it
  }
}

//*******************************************************************|************************************************************//
// return a constant (boost shared) pointer to a system bucket in the bucket data maps
//*******************************************************************|************************************************************//
const SystemBucket_ptr Bucket::fetch_system(const std::string &name) 
                                                              const
{
  SystemBucket_const_it s_it = systems_.find(name);                  // check if a system with this name exists in the bucket
  if (s_it == systems_end())
  {
    dolfin::error(
              "SystemBucket named \"%s\" does not exist in bucket.", // if it doesn't, throw an error
                                                      name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return a pointer to it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_it Bucket::systems_begin()
{
  return systems_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_const_it Bucket::systems_begin() const
{
  return systems_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_it Bucket::systems_end()
{
  return systems_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_const_it Bucket::systems_end() const
{
  return systems_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_it Bucket::orderedsystems_begin()
{
  return orderedsystems_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_const_it Bucket::orderedsystems_begin() const
{
  return orderedsystems_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_it Bucket::orderedsystems_end()
{
  return orderedsystems_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_const_it Bucket::orderedsystems_end() const
{
  return orderedsystems_.end();
}

//*******************************************************************|************************************************************//
// register a base ufl symbol (associated with the derived ufl symbol) in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_baseuflsymbol(const std::string &baseuflsymbol, 
                                    const std::string &uflsymbol)
{
  string_it s_it = baseuflsymbols_.find(uflsymbol);                  // check if this ufl symbol already exists
  if (s_it != baseuflsymbols_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
            "Name with ufl symbol \"%s\" already exists in system.", 
                                              uflsymbol.c_str());
  }
  else
  {
    baseuflsymbols_[uflsymbol] = baseuflsymbol;                      // if it doesn't, assign the baseuflsymbol to the maps
  }
}

//*******************************************************************|************************************************************//
// return a string containing the base ufl symbol for a given ufl symbol from the bucket data maps
//*******************************************************************|************************************************************//
const std::string Bucket::fetch_baseuflsymbol(
                                const std::string &uflsymbol) const
{
  string_const_it s_it = baseuflsymbols_.find(uflsymbol);            // check if this ufl symbol exists
  if (s_it == baseuflsymbols_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "Name with uflsymbol \"%s\" does not exist in system.", 
                                                uflsymbol.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return the string containing the base ufl symbol
  }
}

//*******************************************************************|************************************************************//
// check if a ufl symbol exists in the baseuflsymbols_ bucket data map
//*******************************************************************|************************************************************//
const bool Bucket::contains_baseuflsymbol(
                              const std::string &uflsymbol) const
{
  string_const_it s_it = baseuflsymbols_.find(uflsymbol);
  return s_it != baseuflsymbols_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a function with the given ufl symbol in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_uflsymbol(const std::pair< std::string, GenericFunction_ptr > &uflfunctionpair)
{
  register_uflsymbol(uflfunctionpair.second, uflfunctionpair.first);
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a function with the given ufl symbol in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_uflsymbol(GenericFunction_ptr function, 
                                const std::string &uflsymbol)
{
  GenericFunction_it g_it = uflsymbols_.find(uflsymbol);             // check if the ufl symbol already exists
  if (g_it != uflsymbols_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
    "GenericFunction with ufl symbol \"%s\" already exists in system.", 
                                  uflsymbol.c_str());
  }
  else
  {
    uflsymbols_[uflsymbol] = function;                               // if not, register the pointer in the maps
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to the function associated with the given ufl symbol
//*******************************************************************|************************************************************//
GenericFunction_ptr Bucket::fetch_uflsymbol(
                                const std::string &uflsymbol) const
{
  GenericFunction_const_it g_it = uflsymbols_.find(uflsymbol);       // check if the ufl symbol exists
  if (g_it == uflsymbols_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
    "GenericFunction with uflsymbol \"%s\" does not exist in system.", 
                                          uflsymbol.c_str());
  }
  else
  {
    return (*g_it).second;                                           // if it does, return a pointer to the associated function
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a functionspace for a coefficient with the given ufl symbol in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_coefficientspace(
                                FunctionSpace_ptr coefficientspace, 
                                const std::string &uflsymbol)
{
  FunctionSpace_it f_it = coefficientspaces_.find(uflsymbol);        // check if this ufl symbol already exists
  if (f_it != coefficientspaces_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
    "FunctionSpace with uflsymbol \"%s\" already exists in system coefficientspaces.", 
                                              uflsymbol.c_str());
  }
  else
  {
    coefficientspaces_[uflsymbol] = coefficientspace;                // if not, register the pointer in the coefficientspaces_ map
  }
}

//*******************************************************************|************************************************************//
// check if a functionspace for a coefficient with the given ufl symbol exists in the bucket data maps
//*******************************************************************|************************************************************//
const bool Bucket::contains_coefficientspace(
                                  const std::string &uflsymbol) const
{
  FunctionSpace_const_it f_it = coefficientspaces_.find(uflsymbol);
  return f_it != coefficientspaces_.end();
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a functionspace for a coefficient with the given ufl symbol from the bucket data maps
//*******************************************************************|************************************************************//
FunctionSpace_ptr Bucket::fetch_coefficientspace(
                                  const std::string &uflsymbol) const
{
  FunctionSpace_const_it f_it = coefficientspaces_.find(uflsymbol);  // check if a functionspace with this uflsymbol already exists
  if (f_it == coefficientspaces_.end())
  {
    dolfin::log(dolfin::ERROR, coefficientspaces_str());             // if it doesn't, output which coefficientspaces are available
    dolfin::error(                                                   // and issue an error
    "FunctionSpace with uflsymbol \"%s\" doesn't exist in system coefficientspaces.", 
                                    uflsymbol.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return a pointer to the functionspace
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a detector set in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_detector(GenericDetectors_ptr detector, 
                                            const std::string &name)
{
  GenericDetectors_it d_it = detectors_.find(name);                  // check if a mesh with this name already exists
  if (d_it != detectors_end())
  {
    dolfin::error(
          "Detector set named \"%s\" already exists in bucket.",     // if it does, issue an error
                                                    name.c_str());
  }
  else
  {
    detectors_[name] = detector;                                     // if not, insert it into the detectors_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a detector set from the bucket data maps
//*******************************************************************|************************************************************//
GenericDetectors_ptr Bucket::fetch_detector(const std::string &name)
{
  GenericDetectors_it d_it = detectors_.find(name);                  // check if this detector exists in the detectors_ map
  if (d_it == detectors_end())
  {
    dolfin::error(
          "Detector set named \"%s\" does not exist in bucket.",     // if it doesn't, issue an error
                                                    name.c_str());
  }
  else
  {
    return (*d_it).second;                                           // if it does, return a (boost shared) pointer to it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the detectors_ map
//*******************************************************************|************************************************************//
GenericDetectors_it Bucket::detectors_begin()
{
  return detectors_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the detectors_ map
//*******************************************************************|************************************************************//
GenericDetectors_const_it Bucket::detectors_begin() const
{
  return detectors_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the detectors_ map
//*******************************************************************|************************************************************//
GenericDetectors_it Bucket::detectors_end()
{
  return detectors_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the detectors_ map
//*******************************************************************|************************************************************//
GenericDetectors_const_it Bucket::detectors_end() const
{
  return detectors_.end();
}

//*******************************************************************|************************************************************//
// loop over the systems in the bucket, telling each to output diagnostic data
//*******************************************************************|************************************************************//
void Bucket::output(const int &location)
{

  bool write_vis = (perform_action_(visualization_period_, 
                                    visualization_dumptime_, 
                                  visualization_period_timesteps_) ||// are we writing pvd files?
                     location==OUTPUT_START || location==OUTPUT_END);// force output at start and end

  bool write_stat = (perform_action_(statistics_period_, 
                                     statistics_dumptime_,           // are we writing the stat file? 
                                     statistics_period_timesteps_) ||          
                     location==OUTPUT_START || location==OUTPUT_END);// force output at start and end

  bool write_steady = ((perform_action_(steadystate_period_, 
                            steadystate_dumptime_,                   // are we writing the steady state file?
                            steadystate_period_timesteps_) ||
                      location==OUTPUT_END) &&                       // force output at end
                      location !=OUTPUT_START);                      // prevent output at start

  bool write_det = (perform_action_(detectors_period_, detectors_dumptime_, 
                         detectors_period_timesteps_) ||             // are we writing the detectors file?
                     location==OUTPUT_START || location==OUTPUT_END);// force output at start and end

  if (!write_vis && !write_stat && !write_steady && !write_det)      // if there's nothing to do, just return
  {
    return;
  }  

  if (location != OUTPUT_START)                                      // if we're not at the start of the simulation
  {
    for (int_SystemBucket_const_it s_it = orderedsystems_begin();    // loop over the systems (in order)
                                s_it != orderedsystems_end(); s_it++)
    {
      if((*(*s_it).second).solve_location()==SOLVE_DIAGNOSTICS)      // find if any are meant to be solved before output
      {                                                              // check it has fields included in the current output
        if( (write_vis    && (*(*s_it).second).include_in_visualization()) ||
            (write_stat   && (*(*s_it).second).include_in_statistics())    ||
            (write_steady && (*(*s_it).second).include_in_steadystate())   ||
            (write_det    && (*(*s_it).second).include_in_detectors())        )
        {
          (*(*s_it).second).solve();                                 // solve for those fields
        }
      }
    }
  }


  if (write_stat)
  {
    (*statfile_).write_data();                                       // write data to the statistics file
  }

  if (detfile_ && write_det)
  {
    (*detfile_).write_data();                                        // write data to the detectors file
  }

  if (steadyfile_ && write_steady)
  {
    (*steadyfile_).write_data();                                     // write data to the steady state file
  }

  for (SystemBucket_it s_it = systems_begin(); s_it != systems_end();// loop over the systems
                                                              s_it++)
  {
    (*(*s_it).second).output(write_vis);                             // and output pvd files
  }

}

//*******************************************************************|************************************************************//
// return a string describing the contents of the bucket
//*******************************************************************|************************************************************//
const std::string Bucket::str() const 
{
  std::stringstream s;
  int indent = 1;
  s << "Bucket " << name() << std::endl;
  s << uflsymbols_str(indent);
  s << coefficientspaces_str(indent);
  s << meshes_str(indent);
  s << systems_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the meshes_ data structure
//*******************************************************************|************************************************************//
const std::string Bucket::meshes_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Mesh_const_it m_it = meshes_begin(); m_it != meshes_end(); 
                                                            m_it++ )
  {
    s << indentation << "Mesh " << (*m_it).first  << std::endl;
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the systems_ data structure
// (loop over the systems producing appending strings for each)
//*******************************************************************|************************************************************//
const std::string Bucket::systems_str(const int &indent) const
{
  std::stringstream s;
  for ( SystemBucket_const_it s_it = systems_begin(); 
                                      s_it != systems_end(); s_it++ )
  {
    s << (*(*s_it).second).str(indent);
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing what functionspaces are registered for coefficients in the bucket
//*******************************************************************|************************************************************//
const std::string Bucket::coefficientspaces_str(const int &indent) 
                                                              const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionSpace_const_it f_it = coefficientspaces_.begin(); 
                          f_it != coefficientspaces_.end(); f_it++ )
  {
    s << indentation << "CoefficientSpace for " << (*f_it).first  
                                                      << std::endl;
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing which uflsymbols have functions associated
//*******************************************************************|************************************************************//
const std::string Bucket::uflsymbols_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( GenericFunction_const_it g_it = uflsymbols_.begin(); 
                                  g_it != uflsymbols_.end(); g_it++ )
  {
    if ((*g_it).second)
    {
      s << indentation << "UFLSymbol " << (*g_it).first 
                                      << " associated" << std::endl;
    }
    else
    {
      s << indentation << "UFLSymbol " << (*g_it).first 
                                  << " not associated" << std::endl;
    }
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// after having filled the system and function buckets loop over them and register their functions with their uflsymbols 
//*******************************************************************|************************************************************//
void Bucket::uflsymbols_fill_()
{

  if (!timestep_.first.empty())                                      // if we've registered a timestep (i.e. this is a dynamic
  {                                                                  // simulation)
    register_uflsymbol(timestep_);                                   // register the timestep ufl symbol and function
  }

  for (SystemBucket_const_it s_it = systems_begin();                 // loop over the systems
                                      s_it != systems_end(); s_it++)
  {
    SystemBucket_ptr system = (*s_it).second;
    register_uflsymbol((*system).function(),                         // current system function
                                      (*system).uflsymbol());
    register_uflsymbol((*system).oldfunction(),                      // old system function
                                      (*system).uflsymbol()+"_n");
    register_uflsymbol((*system).iteratedfunction(),                 // iterated system function
                                      (*system).uflsymbol()+"_i");
    for (FunctionBucket_const_it f_it = (*system).fields_begin();    // loop over the fields in this system
                            f_it != (*system).fields_end(); f_it++)
    {
      FunctionBucket_ptr field = (*f_it).second;
      register_uflsymbol((*system).function(),                       // current field
                                      (*field).uflsymbol());
      register_uflsymbol((*system).oldfunction(),                    // old field
                                      (*field).uflsymbol()+"_n");
      register_uflsymbol((*system).iteratedfunction(),               // iterated field
                                      (*field).uflsymbol()+"_i");
    }
    for (FunctionBucket_const_it f_it = (*system).coeffs_begin();    // loop over the coefficients in this system
                            f_it != (*system).coeffs_end(); f_it++)
    {
      FunctionBucket_ptr coeff = (*f_it).second;
      register_uflsymbol((*coeff).function(),                        // current coefficient
                                      (*coeff).uflsymbol());
      register_uflsymbol((*coeff).oldfunction(),                     // old coefficient
                                      (*coeff).uflsymbol()+"_n");
      register_uflsymbol((*coeff).iteratedfunction(),                // iterated coefficient
                                      (*coeff).uflsymbol()+"_i");
    }
  }
}

//*******************************************************************|************************************************************//
// empty the data structures in the bucket
//*******************************************************************|************************************************************//
void Bucket::empty_()
{
  meshes_.clear();
  systems_.clear();
  orderedsystems_.clear();
  baseuflsymbols_.clear();
  uflsymbols_.clear();
  coefficientspaces_.clear();
  detectors_.clear();

  if(statfile_)
  {  
    (*statfile_).close();
  }
}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling solve on each that has requested a solve at the start
//*******************************************************************|************************************************************//
void Bucket::solve_at_start_()
{
  bool systems_solved = false;

  for (int_SystemBucket_const_it s_it = orderedsystems_begin(); 
                              s_it != orderedsystems_end(); s_it++)
  {
    if((*(*s_it).second).solve_location()==SOLVE_START)
    {
      (*(*s_it).second).solve();
      systems_solved = true;
    }
  }

  if(systems_solved)
  {
    output(OUTPUT_START);
  }
}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling solve on each that has requested a solve in the timeloop (within a nonlinear
// systems iteration loop)
//*******************************************************************|************************************************************//
void Bucket::solve_in_timeloop_()
{
  for (*iteration_count_ = 0; \
       *iteration_count_ < nonlinear_iterations(); 
       (*iteration_count_)++)                                        // loop over the nonlinear iterations
  {
    solve(SOLVE_TIMELOOP);                                           // solve all systems in the bucket

    //update_nonlinear();
  }
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the simulation has reached a steady state or not
//*******************************************************************|************************************************************//
bool Bucket::steadystate_()
{
  bool steady = false;

  if (steadystate_tol_)
  {
    double maxchange = 0.0;
    for (SystemBucket_it s_it = systems_begin(); 
                                    s_it != systems_end(); s_it++)
    {
      double systemchange = (*(*s_it).second).maxchange();
      dolfin::log(dolfin::DBG, "  steady state systemchange = %f", systemchange);
      maxchange = std::max(systemchange, maxchange);
    }
    dolfin::log(dolfin::INFO, "steady state maxchange = %f", maxchange);
    steady = (maxchange < *steadystate_tol_);
  }

  return steady;
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the simulation should output or not
//*******************************************************************|************************************************************//
bool Bucket::perform_action_(double_ptr action_period, 
                             double_ptr previous_action_time, 
                             int_ptr action_period_timesteps)
{
  bool performing = true;

  if(action_period)
  {
    assert(previous_action_time);
    performing = ((current_time()-*previous_action_time) >= *action_period);
    if (performing)
    {
      *previous_action_time = current_time();
    }
  }
  else if(action_period_timesteps)
  {
    performing = (timestep_count()%(*action_period_timesteps)==0);
  }  

  return performing;
}

