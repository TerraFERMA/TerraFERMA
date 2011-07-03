import libspud
import sys

def get_function_details(optionpath, system_cell):
  function = {}
  function["name"]   = libspud.get_option(optionpath+"/name")
  function["symbol"] = libspud.get_option(optionpath+"/ufl_symbol")
  function["type"]   = libspud.get_option(optionpath+"/type/name")
  if function["type"] == "Aliased":
    aliasedsystem_name   = libspud.get_option(optionpath+"/type/system")
    aliasedfunction_name = libspud.get_option(optionpath+"/type/generic_function")
    
    aliasedsystem_optionpath = "/system::"+aliasedsystem_name
    if libspud.have_option(aliasedsystem_optionpath+"/field::"+aliasedfunction_name):
      aliasedfunction_optionpath = aliasedsystem_optionpath+"/field::"+aliasedfunction_name
    elif libspud.have_option(aliasedsystem_optionpath+"/coefficient::"+aliasedfunction_name):
      aliasedfunction_optionpath = aliasedsystem_optionpath+"/coefficient::"+aliasedfunction_name
    else:
      print "Unable to find aliased generic_function as either a field or a coefficient."
      sys.exit(1)
    function["type"] = libspud.get_option(aliasedfunction_optionpath+"/type/name")
    if function["type"]=="Aliased":
      print "Can't alias to an aliased function."
      sys.exit(1)
     
    # Perform a check that the meshes at least have the same cell (a necessary but not sufficient check)
    aliasedsystemmesh_name       = libspud.get_option(aliasedsystem_optionpath+"/mesh/name")
    aliasedsystemmesh_optionpath = "/geometry/mesh::"+aliasedsystemmesh_name
    aliasedsystem_cell           = libspud.get_option(aliasedsystemmesh_optionpath+"/cell")
    assert(aliasedsystem_cell == system_cell)
    
    newoptionpath = aliasedfunction_optionpath
  else:
    newoptionpath = optionpath

  function["optionpath"] = optionpath
  function["rank"]   = libspud.get_option(newoptionpath+"/type/rank/name")
  function["family"]   = None
  function["degree"] = None
  if function["type"] != "Constant":
    function["family"] = libspud.get_option(newoptionpath+"/type/rank/element/family")
    function["degree"] = libspud.get_option(newoptionpath+"/type/rank/element/degree")
  
  function["size"]     = None
  function["shape"]    = None
  function["symmetry"] = None
  if function["rank"] == "Vector":
    if libspud.have_option(newoptionpath+"/type/rank/element/size"):
      function["size"] = libspud.get_option(newoptionpath+"/type/rank/element/size")
  elif function["rank"] == "Tensor":
    if libspud.have_option(newoptionpath+"/type/rank/element/shape"):
      function["shape"] = libspud.get_option(newoptionpath+"/type/rank/element/shape")
    if libspud.have_option(newoptionpath+"/type/rank/element/symmetry"):
      function["symmetry"] = True

  return function

def get_diagnosticfunctional_details(optionpath):
  functional = {}
  functional["name"]   = libspud.get_option(optionpath+"/name")
  functional["symbol"] = libspud.get_option(optionpath+"/ufl_symbol")
  functional["form"]   = libspud.get_option(optionpath)+"\n"
  return functional
  
