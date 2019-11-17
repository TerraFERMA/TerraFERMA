# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

from buckettools.base import *
import os

class SystemBucket:
  """A class that stores all the information necessary to write the ufl for a system (i.e. mixed function space)."""

  def __init__(self):
    """Define the expected members of the system class."""
    self.cell = None
    self.mesh_name = None
    self.name = None
    self.symbol = None
    self.fields = None
    self.coeffs = None
    self.special_coeffs = None  # these are constants declared outside systems but still need system information
    self.solvers = None
    self.functionals = None
    self.bucket = None

  def functions_ufl(self):
    """Write an array of ufl strings describing all the functions (fields and coefficients) within a system."""
    ufl = []
    ufl.append(declaration_comment("Function elements", "System", self.name))
    for field in self.fields:
      ufl += field.element_ufl()
    ufl += self.element_ufl()
    ufl.append(os.linesep)
    ufl += self.test_ufl()
    ufl.append(os.linesep)
    ufl += self.trial_ufl()
    ufl.append(os.linesep)
    ufl += self.function_ufl()
    ufl.append(os.linesep)
    ufl += self.iterate_ufl()
    ufl.append(os.linesep)
    ufl += self.old_ufl()
    ufl.append(os.linesep)
    for coeff in self.coeffs:
      if coeff.type == "Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in function_uflsymbol_suffixes():
          ufl.append(coeff.constant_ufl(suffix=suffix))
      else:
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in function_uflsymbol_suffixes():
          ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append(os.linesep)
    # special_coeffs are only added for this system *not the other systems*
    ufl.append(comment("Declaring special coefficients, such as the timestep."))
    ufl.append(os.linesep)
    for coeff in self.special_coeffs:
      if coeff.type == "Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in function_uflsymbol_suffixes():
          ufl.append(coeff.constant_ufl(suffix=suffix))
      else:
        if coeff.type == "Function":
          print("coefficient functionspaces not output for special coefficient functions")
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in function_uflsymbol_suffixes():
          ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append(os.linesep)
    ufl.append(comment("Finished declaring functions for this system, start on other systems."))
    ufl.append(os.linesep)
    for system in self.bucket.systems:
      if system.name == self.name: continue
      ufl.append(declaration_comment("Function elements", "System", system.name))
      for field in system.fields:
        ufl += field.element_ufl()
      ufl += system.element_ufl()
      ufl.append(os.linesep)
      ufl += system.function_ufl()
      ufl.append(os.linesep)
      ufl += system.iterate_ufl()
      ufl.append(os.linesep)
      ufl += system.old_ufl()
      ufl.append(os.linesep)
      for coeff in system.coeffs:
        if coeff.type == "Constant":
          ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
          for suffix in function_uflsymbol_suffixes():
            ufl.append(coeff.constant_ufl(suffix=suffix))
        else:
          ufl += coeff.element_ufl()
          ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
          for suffix in function_uflsymbol_suffixes():
            ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append(os.linesep)
    ufl.append(comment("Finished declaring functions for all other systems, start on forms."))
    ufl.append(os.linesep)
    return ufl

  def element_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) element of a system."""
    ufl = []
    if len(self.fields)==0: return ufl
    if len(self.fields)==1:
      ufl.append(comment("System element is not mixed"))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_e")+os.linesep)
    else:
      ufl.append(declaration_comment("Mixed element", "System", self.name))
      ufl_line = self.symbol+"_e = "
      ufl_line += "MixedElement(["
      for s in range(len(self.fields)-1):
        ufl_line += self.fields[s].symbol+"_e, "
      ufl_line += self.fields[-1].symbol+"_e"
      ufl_line += "])"+os.linesep
      ufl.append(ufl_line)
    return ufl

  def test_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) test space of a system."""
    ufl = []
    if len(self.fields)==0: return ufl
    if len(self.fields)==1:
      ufl.append(declaration_comment("Test space", self.fields[0].type, self.fields[0].name))
      ufl.append(testfunction_ufl(self.fields[0].symbol))
      ufl.append(declaration_comment("Test space", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_t")+os.linesep)
    else:
      ufl.append(declaration_comment("Test space", "System", self.name))
      ufl.append(testfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Test spaces", "functions", ""))
      ufl += self.split_ufl(suffix="_t")
    return ufl

  def trial_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) trial space of a system."""
    ufl = []
    if len(self.fields)==0: return ufl
    if len(self.fields)==1:
      ufl.append(declaration_comment("Trial space", self.fields[0].type, self.fields[0].name))
      ufl.append(trialfunction_ufl(self.fields[0].symbol))
      ufl.append(declaration_comment("Trial space", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_a")+os.linesep)
    else:
      ufl.append(declaration_comment("Trial space", "System", self.name))
      ufl.append(trialfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Trial spaces", "functions", ""))
      ufl += self.split_ufl(suffix="_a")
    return ufl

  def function_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) field values of a system."""
    ufl = []
    if len(self.fields)==0: return ufl
    if len(self.fields)==1:
      ufl.append(declaration_comment("Value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol))
      ufl.append(declaration_comment("Value", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol)+os.linesep)
    else:
      ufl.append(declaration_comment("Value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol))
      ufl.append(declaration_comment("Values", "functions", ""))
      ufl += self.split_ufl()
    return ufl

  def iterate_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) iterated field values of a system."""
    ufl = []
    if len(self.fields)==0: return ufl
    if len(self.fields)==1:
      ufl.append(declaration_comment("Last iteration value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_i"))
      ufl.append(declaration_comment("Last iteration value", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_i")+os.linesep)
    else:
      ufl.append(declaration_comment("Last iteration value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol, suffix="_i"))
      ufl.append(declaration_comment("Last iteration values", "functions", ""))
      ufl += self.split_ufl(suffix="_i")
    return ufl

  def old_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) old field values of a system."""
    ufl = []
    if len(self.fields)==0: return ufl
    if len(self.fields)==1:
      ufl.append(declaration_comment("Previous time-level value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_n"))
      ufl.append(declaration_comment("Previous time-level value", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_n")+os.linesep)
    else:
      ufl.append(declaration_comment("Previous time-level value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol, suffix="_n"))
      ufl.append(declaration_comment("Previous time-level values", "functions", ""))
      ufl += self.split_ufl(suffix="_n")
    return ufl

  def split_ufl(self, suffix=""):
    """Write an array of ufl strings splitting a mixed function or element into components for its constitutive fields."""
    ufl = []
    for s in range(len(self.fields)):
      ufl.append(comment(" - "+self.fields[s].name))
    ufl_line = "("
    for s in range(len(self.fields)-1):
      ufl_line += self.fields[s].symbol+suffix+", "
    ufl_line += self.fields[-1].symbol+suffix
    ufl_line += ") = split("+self.symbol+suffix+")"+os.linesep
    ufl.append(ufl_line)
    return ufl

  def include_systemfunctionals_cpp(self):
    """Write an array of cpp strings including the namespaces of the system ufls."""
    cpp = []
    for coeff in self.coeffs:
      if coeff.functional:
        cpp.append("#include \""+coeff.functional.namespace()+".h\""+os.linesep)
    for functional in self.functionals:
      cpp.append("#include \""+functional.namespace()+".h\""+os.linesep)
    return cpp

  def include_systemsolvers_cpp(self):
    """Write an array of cpp strings including the namespaces of the system ufls."""
    cpp = []
    for solver in self.solvers:
      cpp.append("#include \""+solver.namespace()+".h\""+os.linesep)
    return cpp

  def include_systemexpressions_cpp(self):
    """Write an array of cpp strings including the namespaces of the system ufls."""
    cpp = []
    for field in self.fields:
      for cppexpression in field.cpp:
        cpp.append("#include \""+cppexpression.namespace()+".h\""+os.linesep)
    for coeff in self.coeffs:
      for cppexpression in coeff.cpp:
        cpp.append("#include \""+cppexpression.namespace()+".h\""+os.linesep)
    return cpp

  def functionspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the functionspaces in the ufls (assuming a single solver)."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    cpp.append("      // All solvers within a system should return the same functionspace so just take the first one"+os.linesep)
    if len(self.solvers)==0:
      cpp.append("      tf_err(\"Unknown system functionspace in ufc_fetch_functionspace\", \"System name: %s\", systemname.c_str());"+os.linesep)
    else:
      cpp += self.solvers[0].functionspace_cpp_no_if()
    cpp.append("    }"+os.linesep)
    return cpp

  def solverfunctionspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the functionspaces in the ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    cpp.append("      // All solvers within a system should return the same functionspace"+os.linesep)
    for s in range(len(self.solvers)):
      cpp += self.solvers[s].functionspace_cpp(index=s)
    if len(self.solvers)==0:
      cpp.append("      tf_err(\"No solvers in system in ufc_fetch_functionspace\", \"Solver name: %s\", solvername.c_str());"+os.linesep)
    else:
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown solvername in ufc_fetch_functionspace\", \"Solver name: %s\", solvername.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
    cpp.append("    }"+os.linesep)
    return cpp

  def functionalcoefficientspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    if len(self.functionals)==0:
      cpp.append("        tf_err(\"No functionals in ufc_fetch_coefficientspace_from_functional\", \"Functional name: %s\", functionalname.c_str());"+os.linesep)
    else:
      for f in range(len(self.functionals)):
        if f == 0:
          cpp.append("        if (functionalname ==  \""+self.functionals[f].name+"\")"+os.linesep)
        else:
          cpp.append("        else if (functionalname ==  \""+self.functionals[f].name+"\")"+os.linesep)
        cpp.append("        {"+os.linesep)
        symbols_in_form = form_symbols(self.functionals[f].namespace())
        symbols_found = 0
        for c2 in range(len(self.coeffs)):
          for suffix in function_uflsymbol_suffixes():
            if self.coeffs[c2].symbol+suffix in symbols_in_form: 
              cpp += self.functionals[f].coefficientspace_cpp(self.coeffs[c2], index=symbols_found, suffix=suffix)
              symbols_found =+ 1
              break
        if symbols_found == 0:
          cpp.append("          tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
          cpp.append("        }"+os.linesep)
        else:
          cpp.append("          else"+os.linesep)
          cpp.append("          {"+os.linesep)
          cpp.append("            tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
          cpp.append("          }"+os.linesep)
          cpp.append("        }"+os.linesep)
      cpp.append("        else"+os.linesep)
      cpp.append("        {"+os.linesep)
      cpp.append("          tf_err(\"Unknown functionalname in ufc_fetch_coefficientspace_from_functional\", \"Functional name: %s\", functionalname.c_str());"+os.linesep)
      cpp.append("        }"+os.linesep)
    cpp.append("    }"+os.linesep)
    return cpp

  def constantfunctionalcoefficientspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the constant functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    if (len(self.coeffs)==0):
      cpp.append("      tf_err(\"No coefficients in ufc_fetch_coefficientspace_from_constant_functional\", \"Coefficient name: %s\", coefficientname.c_str());"+os.linesep)
      cpp.append("    }"+os.linesep)
    else:
      for c in range(len(self.coeffs)):
        if c == 0:
          cpp.append("      if (coefficientname ==  \""+self.coeffs[c].name+"\")"+os.linesep)
        else:
          cpp.append("      else if (coefficientname ==  \""+self.coeffs[c].name+"\")"+os.linesep)
        cpp.append("      {"+os.linesep)
        if self.coeffs[c].functional:
          symbols_in_form = form_symbols(self.coeffs[c].functional.namespace())
          symbols_found = 0
          for c2 in range(len(self.coeffs)):
            for suffix in function_uflsymbol_suffixes():
              if self.coeffs[c2].symbol+suffix in symbols_in_form: 
                cpp += self.coeffs[c].functional.coefficientspace_cpp(self.coeffs[c2], index=symbols_found, suffix=suffix)
                symbols_found =+ 1
                break
          if symbols_found == 0:
            cpp.append("        tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_constant_functional\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
            cpp.append("      }"+os.linesep)
          else:
            cpp.append("          else"+os.linesep)
            cpp.append("          {"+os.linesep)
            cpp.append("            tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_constant_functional\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
            cpp.append("          }"+os.linesep)
            cpp.append("        }"+os.linesep)
        else:
          cpp.append("        tf_err(\"No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional\", \"Coefficient name: %s\", coefficientname.c_str());"+os.linesep)
          cpp.append("      }"+os.linesep)
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown coefficientname in ufc_fetch_coefficientspace_from_constant_functional\", \"Coefficient name: %s\", coefficientname.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
      cpp.append("    }"+os.linesep)
    return cpp

  def solvercoefficientspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    for s in range(len(self.solvers)):
      if s == 0:
        cpp.append("      if (solvername ==  \""+self.solvers[s].name+"\")"+os.linesep)
      else:
        cpp.append("      else if (solvername ==  \""+self.solvers[s].name+"\")"+os.linesep)
      cpp.append("      {"+os.linesep)
      if len(self.coeffs)==0:
        cpp.append("        tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
        cpp.append("      }"+os.linesep)
      else:
        symbols_in_forms = form_symbols(self.solvers[s].namespace())
        symbols_found = 0
        for c in range(len(self.coeffs)):
          for suffix in function_uflsymbol_suffixes():
            if self.coeffs[c].symbol+suffix in symbols_in_forms: 
              cpp += self.coeffs[c].solvercoefficientspace_cpp(self.solvers[s].name, index=symbols_found, suffix=suffix)
              symbols_found =+ 1
              break
        if symbols_found == 0:
          cpp.append("        tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
          cpp.append("      }"+os.linesep)
        else:
          cpp.append("        else"+os.linesep)
          cpp.append("        {"+os.linesep)
          cpp.append("          tf_err(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver\", \"UFL symbol: %s\", uflsymbol.c_str());"+os.linesep)
          cpp.append("        }"+os.linesep)
          cpp.append("      }"+os.linesep)
    if len(self.solvers)==0:
      cpp.append("      tf_err(\"No solver in system in ufc_fetch_coefficientspace_from_solver\", \"Solver name: %s\", solvername.c_str());"+os.linesep)
    else:
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown solvername in ufc_fetch_coefficientspace_from_solver\", \"Solver name: %s\", solvername.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
    cpp.append("    }"+os.linesep)
    return cpp

  def form_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the forms in the ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    for s in range(len(self.solvers)):
      if s == 0:
        cpp.append("      if (solvername ==  \""+self.solvers[s].name+"\")"+os.linesep)
      else:
        cpp.append("      else if (solvername ==  \""+self.solvers[s].name+"\")"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        if (solvertype == \""+self.solvers[s].type+"\")"+os.linesep)
      cpp.append("        {"+os.linesep)
      cpp += self.solvers[s].form_cpp()
      cpp.append("          else"+os.linesep)
      cpp.append("          {"+os.linesep)
      cpp.append("            tf_err(\"Unknown formname in ufc_fetch_form\", \"Form name: %s\", formname.c_str());"+os.linesep)
      cpp.append("          }"+os.linesep)
      cpp.append("        }"+os.linesep)
      cpp.append("        else"+os.linesep)
      cpp.append("        {"+os.linesep)
      cpp.append("          tf_err(\"Unknown solvertype in ufc_fetch_form\", \"Form name: %s\", formname.c_str());"+os.linesep)
      cpp.append("        }"+os.linesep)
      cpp.append("      }"+os.linesep)
    if len(self.solvers)==0:
      cpp.append("      tf_err(\"No solver in system in ufc_fetch_form\", \"Solver name: %s\", solvername.c_str());"+os.linesep)
    else:
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown systemname in ufc_fetch_form\", \"System name: %s\", systemname.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
    cpp.append("    }"+os.linesep)
    return cpp

  def functional_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    for a in range(len(self.functionals)):
      cpp += self.functionals[a].cpp(index=a)
    if len(self.functionals)==0:
      cpp.append("        tf_err(\"No functionals in ufc_fetch_functional\", \"Functional name: %s\", functionalname.c_str());"+os.linesep)
    else:
      cpp.append("        else"+os.linesep)
      cpp.append("        {"+os.linesep)
      cpp.append("          tf_err(\"Unknown functionalname in ufc_fetch_functional\", \"Functional name: %s\", functionalname.c_str());"+os.linesep)
      cpp.append("        }"+os.linesep)
    cpp.append("    }"+os.linesep)
    return cpp

  def constantfunctional_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the constant coefficient functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    functionals_found = 0
    for c in range(len(self.coeffs)):
      if self.coeffs[c].functional:
        if functionals_found == 0:
          cpp.append("      if (coefficientname ==  \""+self.coeffs[c].name+"\")"+os.linesep)
        else:
          cpp.append("      else if (coefficientname ==  \""+self.coeffs[c].name+"\")"+os.linesep)
        functionals_found =+ 1
        cpp.append("      {"+os.linesep)
        cpp.append("        functional.reset(new "+self.coeffs[c].functional.namespace()+"::Form_"+self.coeffs[c].functional.symbol+"(mesh));"+os.linesep)
        cpp.append("      }"+os.linesep)
    if functionals_found==0:
      cpp.append("      tf_err(\"Unknown coefficientname in ufc_fetch_constant_functional\", \"Coefficient name: %s\", coefficientname.c_str());"+os.linesep)
      cpp.append("    }"+os.linesep)
    else:
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown coefficientname in ufc_fetch_constant_functional\", \"Coefficient name: %s\", coefficientname.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
      cpp.append("    }"+os.linesep)
    return cpp

  def cppexpression_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the cpp expressions."""
    cpp = []  
    
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    functions = 0
    for c in range(len(self.fields)):
      cpp += self.fields[c].cppexpression_cpp(index=c)
    for c in range(len(self.coeffs)):
      cpp += self.coeffs[c].cppexpression_cpp(index=c+len(self.fields))
    if (len(self.fields)+len(self.coeffs))==0:
      cpp.append("      tf_err(\"Unknown functionname in cpp_fetch_expression.\", \"Function name: %s\", functionname.c_str());"+os.linesep)
      cpp.append("    }"+os.linesep)
    else:
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown functionname in cpp_fetch_expression.\", \"Function name: %s\", functionname.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
      cpp.append("    }"+os.linesep)
    
    return cpp

  def cppexpression_init(self, index=0):
    """Write an array of cpp strings describing the namespace of the cpp expressions."""
    cpp = []  
    
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("    {"+os.linesep)
    expressions_found = 0
    for c in range(len(self.fields)):
      cpp += self.fields[c].cppexpression_init(index=c)
    for c in range(len(self.coeffs)):
      cpp += self.coeffs[c].cppexpression_init(index=c+len(self.fields))
    if (len(self.fields)+len(self.coeffs))==0:
      cpp.append("      tf_err(\"Unknown functionname in cpp_init_expression.\", \"Function name: %s\", functionname.c_str());"+os.linesep)
      cpp.append("    }"+os.linesep)
    else:
      cpp.append("      else"+os.linesep)
      cpp.append("      {"+os.linesep)
      cpp.append("        tf_err(\"Unknown functionname in cpp_init_expression.\", \"Function name: %s\", functionname.c_str());"+os.linesep)
      cpp.append("      }"+os.linesep)
      cpp.append("    }"+os.linesep)
    
    return cpp

