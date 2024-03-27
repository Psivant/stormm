// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_INTAKE_H
#define STORMM_ATOMGRAPH_INTAKE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "atomgraph_enumerators.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

/// \brief Read and return one or more AtomGraph objects based on a list of file names.  If any of
///        the files are not found, indicate the problem.
///
/// Overloaded:
///   - Try to read in a single file and return the resulting AtomGraph object
///   - Try to read in multiple files and return a list of the resulting AtomGraph objects
///
/// \param file_name    Name of the file to read as the topology
/// \param file_names   List of names to read as topologies
/// \param priority     Indicate the action to take if files are not found
/// \param files_found  Indicator of whether the files were found (modified and returned if not
///                     the null pointer)
/// \{
AtomGraph loadTopology(const std::string &file_name, bool *files_found = nullptr,
                       ExceptionResponse priority = ExceptionResponse::WARN,
                       TopologyKind engine_format = TopologyKind::AMBER);

std::vector<AtomGraph> loadTopology(const std::vector<std::string> &file_names,
                                    bool *files_found = nullptr,
                                    ExceptionResponse priority = ExceptionResponse::WARN,
                                    TopologyKind engine_format = TopologyKind::AMBER);
/// \}
  
} // namespace topology
} // namespace stormm

#endif
