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

