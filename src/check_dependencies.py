import sys, os, argparse, copy

#--------------------------------------------------------------------------------------------------
def findHeaderFiles(directory, header_ext):
  """Find all files with extensions matching one of a list of designated for header files."""
  header_files = []
  for root, _, files in os.walk(directory):
    for file in files:
      for ext in header_ext:
        if file.endswith(ext):
          header_files.append(os.path.join(root, file))
  return header_files

#--------------------------------------------------------------------------------------------------
def seekIncludes(filename, header_exts):
  """Seek the inclusion pragmas within a library header file."""
  incl_list = []
  with open(filename, 'r') as f:
    for line in f:
      if (line.startswith('#') and ('include' in line) and ('\"' in line)):
        found_ext = False
        for ext in header_exts:
          if (ext in line):
            found_ext = True
        if (found_ext):
          dep_start = 0
          dep_end = 0
          depcon = 0
          on_dep = False
          for letter in line:
            if (not on_dep and letter == '\"'):
              on_dep = True
              dep_start = depcon + 1
            if (on_dep and letter == '\"'):
              dep_end = depcon
            depcon += 1
          if (dep_end > dep_start):
            incl_list.append(line[dep_start:dep_end])
  return incl_list

#--------------------------------------------------------------------------------------------------
def buildDependencyList(header_files, header_exts):
  """Build the list of dependencies for each header file."""
  result = []
  for item in header_files:
    result.append((item, seekIncludes(item, header_exts)))
  return result

#--------------------------------------------------------------------------------------------------
def getMajorPathElements(filename):
  """Assuming that a string is the path to a file, break it into its major elements as a list  """
  """of directory names, terminating in the file name itself.  Remove leading 'this directory' """
  """specifications.                                                                           """
  tmp_path = filename.split(os.sep)
  fi_path = []
  for item in tmp_path:
    if (item != '.'):
      fi_path.append(item)
  return fi_path
  
#--------------------------------------------------------------------------------------------------
def codifyDependencies(header_files, deps_list, permit_cases):
  """Codify each header files with a unique integer based on their order of appearance in the """
  """master file list, then assign each dependcy within any file an integer corresponding to  """
  """the best match found among all known header files.                                       """
  result = []
  delimited_hdrs = []
  delimited_incls = []
  osc = os.sep
  for icon, item in enumerate(header_files):
    delimited_hdrs.append(getMajorPathElements(item))
  for icon, item in enumerate(deps_list):
    t_incls = []
    for incl in item[1]:
      t_incls.append(getMajorPathElements(incl))
    delimited_incls.append(t_incls)
    result.append([ -1 ] * len(item[1]))
  n_head = len(header_files)
  head_con = 0
  for head_con, dep_set in enumerate(delimited_incls):
    for depn_con, depn in enumerate(dep_set):
      n_depn_elem = len(depn)
      
      # Match this dependency to one of the known header files.  A match will be declared when
      # each major part of the dependency's path is identified in the delimited parts of a known
      # header file's path, working back to front.
      matched = False
      if (depn[-1] in permit_cases):
        matched = True
        continue
      i = 0
      while (not matched and i < n_head):
        n_head_elem = len(delimited_hdrs[i])
        if (n_head_elem > n_depn_elem):
          n_comp = n_depn_elem
          n_jump = n_head_elem - n_depn_elem
          matched = True
          for j in range(n_comp):
            if (depn[j] != delimited_hdrs[i][j + n_jump]):
              matched = False
          if (matched):
            result[head_con][depn_con] = i
        else:
          n_comp = n_head_elem
          n_jump = n_depn_elem - n_head_elem
          matched = True
          for j in range(n_comp):
            if (depn[j + n_jump] != delimited_hdrs[i][j]):
              matched = False
          if (matched):
            result[head_con][depn_con] = i
        i += 1
      if (not matched):
        print("Dependency matched no known header files:")
        print(depn)
  return result

#--------------------------------------------------------------------------------------------------
def traceBackOneLevel(deps_index, deps_trace, header_files):

  # Check for null dependencies.  These relate to "permitted cases" in which the dependency in
  # question is not expected to lie within the known header files.
  if (deps_trace[-1] == -1):
    return
  
  # Go to the last dependency on the stack and loop over each of its own dependencies
  next_dep = deps_trace[-1]
  n_options = len(deps_index[next_dep])

  # CHECK
#  print("Trace is:")
#  print(deps_trace)
#  print()
#  print("Next comes...")
#  print(deps_index[next_dep])
#  print()
  # END CHECK
  
  for i in range(n_options):

    # Detect a circular dependency
    if (deps_index[next_dep][i] in deps_trace):
      print("Circular dependency detected:")
      for idx in deps_trace:
        cont_str = ""
        if (idx < len(deps_trace) - 1):
          cont_str = " ->"
        print("  " + header_files[idx] + cont_str)
      print("  This then includes " + header_files[deps_index[next_dep][i]])
      quit()
    deps_trace.append(deps_index[next_dep][i])
    traceBackOneLevel(deps_index, deps_trace, header_files)

    # Remove the last element and try a new one
    deps_trace.pop()

#--------------------------------------------------------------------------------------------------
if (__name__ == "__main__"):
  """Main"""
  parser = argparse.ArgumentParser()
  parser.add_argument('-dir', help = 'Restrict or redirect the path to scan.  The default ' +
                      'behavior is to scan the entire current directory and any sub-directories.',
                      default = '.')
  parser.add_argument('-header', '--header', help = 'Add a new file extension to be considered ' +
                      'as a header file.  The default behavior is to consider files ending in ' +
                      '.h and .cuh.', default = ['.h', '.cuh'], action = 'append');
  parser.add_argument('-permit', '--permit', help = 'Specify the name of a dependency that may ' +
                      'not exist in the working tree but will be found through compiler-level ' +
                      'include arguments.', default = [ 'pocketfft_hdronly.h', 'cublas_v2.h' ],
                      action = 'append')
  args = parser.parse_args()
  header_exts = args.header
  permit_cases = args.permit
  header_files = findHeaderFiles(args.dir, header_exts)
  deps_list = buildDependencyList(header_files, header_exts)
  deps_index = codifyDependencies(header_files, deps_list, permit_cases)

  # Search the codified dependencies for circular references
  n_hdr = len(deps_index)
  n_dep = []
  for i in range(n_hdr):
    n_dep.append(len(deps_index[i]))
  for i in range(n_hdr):
    for j in range(n_dep[i]):

      # Initiate a chain
      deps_trace = [ deps_index[i][j] ]
      
      # Call a recursive function to search one dependency up the tree with each call on the stack.
      traceBackOneLevel(deps_index, deps_trace, header_files)
