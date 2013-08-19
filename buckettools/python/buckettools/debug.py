# (originally part of fluidity, modified for buckettools by Cian Wilson)

# Copyright (C) 2006 Imperial College London and others.
# 
# Please see the AUTHORS file in the main source directory for a full list
# of copyright holders.
#
# Prof. C Pain
# Applied Modelling and Computation Group
# Department of Earth Science and Engineering
# Imperial College London
#
# amcgsoftware@imperial.ac.uk
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation,
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

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


def shell():
  '''
  shell()

  Return ipython shell. To actually start the shell, invoke the function
  returned by this function.

  This is particularly useful for debugging embedded
  python or for crawling over the data when something has gone wrong.
  '''
  import sys
  
  if not hasattr(sys,"argv"):
    sys.argv=[]

  try:
    from IPython.Shell import IPShellEmbed
  except ImportError:
    sys.stderr.write(
      """
      *****************************************************
      *** Failed to import IPython. This probably means ***
      *** you don't have it installed. Please install   ***
      *** IPython and try again.                        ***
      *****************************************************
      """)
    raise

  banner = """
  This is an IPython shell embedded in buckettools. You can use it to examine
  or even set variables. Press CTRL+d to exit and return to your program.
  """

  ipshell = IPShellEmbed(banner=banner)

  return ipshell
  
