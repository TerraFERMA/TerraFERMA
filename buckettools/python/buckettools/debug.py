
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
  
