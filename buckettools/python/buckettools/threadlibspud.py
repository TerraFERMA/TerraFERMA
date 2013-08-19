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

import threading
import libspud

class ThreadSpudLock:
  '''A thread-safe way of loading libspud.'''
  def __init__(self):
    self.lock=threading.Lock()

  def load_options(self, filename):
    self.lock.acquire()
    libspud.load_options(filename)

  def clear_options(self):
    libspud.clear_options()
    self.lock.release()

threadlibspud = ThreadSpudLock()

