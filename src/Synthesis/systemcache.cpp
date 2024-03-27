#include "copyright.h"
#include "systemcache.h"
#include "FileManagement/file_listing.h"
#include "Math/series_ops.h"
#include "Math/summation.h"
#include "Parsing/parse.h"
#include "Potential/scorecard.h"
#include "Potential/valence_potential.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"

namespace stormm {
namespace synthesis {

using diskutil::getBaseName;
using diskutil::osSeparator;
using diskutil::splitPath;
using energy::evaluateBondTerms;
using energy::evaluateAngleTerms;
using energy::ScoreCard;
using stmath::incrementingSeries;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::sum;
using namelist::MoleculeSystem;
using parse::findStringInVector;
using structure::readStructureDataFile;
using topology::UnitCellType;
using topology::ValenceKit;
using trajectory::detectCoordinateFileKind;
using trajectory::getEnumerationName;
  
//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache(const ExceptionResponse policy_in,
                         const MapRotatableGroups map_chemfe_rotators,
                         const PrintSituation expectation_in) :
    policy{policy_in}, expectation{expectation_in}, system_count{0}, topology_count{0},
    label_count{0}, restraint_count{0}, file_merger_protocol{TrajectoryFusion::AUTO},
    topology_cache{}, coordinates_cache{}, sdf_cache{}, label_cache{}, features_cache{},
    restraints_cache{}, static_masks_cache{}, forward_masks_cache{}, topology_indices{},
    label_indices{}, label_degeneracy{}, trajectory_subindices{}, trajectory_degeneracy{},
    checkpoint_subindices{}, checkpoint_degeneracy{}, restraint_indices{}, example_indices{},
    topology_cases{}, topology_case_bounds{}, system_input_coordinate_names{},
    system_trajectory_names{}, system_checkpoint_names{}, system_labels{}, label_cases{},
    label_case_bounds{}, system_input_coordinate_kinds{}, system_trajectory_kinds{},
    system_checkpoint_kinds{}
{}

//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache(const FilesControls &fcon, const std::vector<RestraintControls> &rstcon,
                         const ExceptionResponse policy_in,
                         const MapRotatableGroups map_chemfe_rotators,
                         const PrintSituation expectation_in,
                         StopWatch *timer_in) :
    SystemCache(policy_in, map_chemfe_rotators, expectation_in)
{
  // Read critical details from the namelist that will guide the I/O
  file_merger_protocol = fcon.getFileFusionProtocol();
  
  // Read all free topologies, using a try-catch block to filter out things that may not work.
  int n_free_top = fcon.getFreeTopologyCount();
  topology_cache.reserve(n_free_top);
  for (int i = 0; i < n_free_top; i++) {
    try {
      topology_cache.push_back(AtomGraph(fcon.getFreeTopologyName(i)));
    }
    catch (std::runtime_error) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The format of topology " + fcon.getFreeTopologyName(i) +
              " could not be understood.", "SystemCache");
      case ExceptionResponse::WARN:
        rtWarn("The format of topology " + fcon.getFreeTopologyName(i) +
               " could not be understood.  The file will be skipped.", "SystemCache");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
  n_free_top = topology_cache.size();
  
  // Read all free coordinate sets, using file type detection to filter the list
  int n_free_crd = fcon.getFreeCoordinatesCount();
  std::vector<CoordinateFrame> tmp_coordinates_cache;
  tmp_coordinates_cache.reserve(n_free_crd);
  std::vector<int> tmp_coordinates_frame_count;
  tmp_coordinates_frame_count.reserve(n_free_crd);
  std::vector<CoordinateFileKind> tmp_coordinates_kind;
  tmp_coordinates_kind.reserve(n_free_crd);
  for (int i = 0; i < n_free_crd; i++) {
    const std::string crd_name = fcon.getFreeCoordinateName(i);
    const CoordinateFileKind kind = detectCoordinateFileKind(crd_name);
    switch (kind) {
    case CoordinateFileKind::AMBER_CRD:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Amber .crd format trajectories such as " + crd_name + " cannot be read without a "
              "means of providing the atom count.  Use the -sys keyword to couple such "
              "trajectories to a specific topology file so that multiple frames can be read from "
              "such a format.", "SystemCache");
      case ExceptionResponse::WARN:
        rtWarn("Amber .crd format trajectories such as " + crd_name + " cannot be read without a "
               "means of providing the atom count.  Use the -sys keyword to couple such "
               "trajectories to a specific topology file so that multiple frames can be read from "
               "such a format.  These frames will be skipped.", "SystemCache");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    case CoordinateFileKind::AMBER_INPCRD:
    case CoordinateFileKind::AMBER_ASCII_RST:
      try {
        tmp_coordinates_cache.push_back(CoordinateFrame(crd_name, kind));
        tmp_coordinates_frame_count.push_back(1);
        tmp_coordinates_kind.push_back(kind);
      }
      catch (std::runtime_error) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The file " + crd_name + " could not be added to the temporary coordinates cache "
                "in order to pair these orphan coordinate systems to a topology.", "SystemCache");
        case ExceptionResponse::WARN:
          rtWarn("Frames from the SD file " + crd_name + " could not be added to the temporary "
                 "coordinates cache in order to pair these orphan coordinate systems to a "
                 "topology.  These frames will be skipped.", "SystemCache");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      break;
    case CoordinateFileKind::AMBER_NETCDF:
    case CoordinateFileKind::AMBER_NETCDF_RST:
      break;
    case CoordinateFileKind::SDF:
      try {
        const std::vector<MdlMol> all_frames = readStructureDataFile(crd_name);
        const size_t n_entries = all_frames.size();
        for (size_t i = 0; i < n_entries; i++) {
          tmp_coordinates_cache.push_back(all_frames[i].exportCoordinateFrame());
          tmp_coordinates_frame_count.push_back(1);
          tmp_coordinates_kind.push_back(kind);
        }
      }
      catch (std::runtime_error) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("Frames from the SD file " + crd_name + " could not be added to the temporary "
                "coordinates cache in order to pair these orphan coordinate systems to a "
                "topology.", "SystemCache");
        case ExceptionResponse::WARN:
          rtWarn("Frames from the SD file " + crd_name + " could not be added to the temporary "
                 "coordinates cache in order to pair these orphan coordinate systems to a "
                 "topology.  These frames will be skipped.", "SystemCache");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      break;
    case CoordinateFileKind::UNKNOWN:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The format of " + crd_name + " could not be understood.", "SystemCache");
      case ExceptionResponse::WARN:
        rtWarn("The format of " + crd_name + " could not be understood.  The file will be skipped",
               "SystemCache");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
  n_free_crd = tmp_coordinates_cache.size();

  // Filter the unique topologies and match them to systems.  List the atom counts of each.
  int max_match = 0;
  std::vector<int> topology_atom_counts(n_free_top);
  std::vector<int> coordinate_atom_counts(n_free_crd);
  for (int i = 0; i < n_free_top; i++) {
    topology_atom_counts[i] = topology_cache[i].getAtomCount();
  }
  for (int i = 0; i < n_free_crd; i++) {
    coordinate_atom_counts[i] = tmp_coordinates_cache[i].getAtomCount();
  }

  // Make a table of the unique topology atom counts.  Track which coordinate sets might
  // be tied to each group of topologies.
  std::vector<bool> topology_covered(n_free_top, false);
  std::vector<bool> coordinates_covered(n_free_crd, false);
  std::vector<int> topology_series(n_free_top);
  std::vector<int> coordinate_series(n_free_top);
  std::vector<int> unique_topology_sizes;
  std::vector<int> unique_topology_size_bounds(1, 0);
  std::vector<int> unique_coordinate_size_bounds(1, 0);
  int top_series_counter = 0;
  int crd_series_counter = 0;
  for (int i = 0; i < n_free_top; i++) {
    if (topology_covered[i]) {
      continue;
    }

    // First, scan to see if any other topologies share the same number of atoms.
    topology_covered[i] = true;
    int n_samesize_topology = 1;
    topology_series[top_series_counter] = i;
    top_series_counter++;
    const int iatom_count = topology_atom_counts[i];
    for (int j = i + 1; j < n_free_top; j++) {
      if (topology_atom_counts[j] == iatom_count) {
        topology_covered[j] = true;
        topology_series[top_series_counter] = j;
        top_series_counter++;
      }
    }
    for (int j = 0; j < n_free_crd; j++) {
      if (coordinates_covered[j] == false && coordinate_atom_counts[j] == iatom_count) {
        coordinates_covered[j] = true;
        coordinate_series[crd_series_counter] = j;
        crd_series_counter++;
      }
    }
    unique_topology_sizes.push_back(iatom_count);
    unique_topology_size_bounds.push_back(top_series_counter);
    unique_coordinate_size_bounds.push_back(crd_series_counter);
  }

  // Check that all coordinates were covered
  std::vector<std::string> orphan_coordinates;
  for (int i = 0; i < n_free_crd; i++) {
    if (coordinates_covered[i] == false) {
      orphan_coordinates.push_back(tmp_coordinates_cache[i].getFileName());
    }
  }
  const int n_orphan = orphan_coordinates.size();
  if (n_orphan > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      {
        std::string orphan_errmsg;
        if (n_orphan > 6) {
          for (int i = 0; i < 3; i++) {
            orphan_errmsg += "  " + orphan_coordinates[i] + '\n';
          }
          orphan_errmsg += "  (... more files, list truncated ...)\n";
          for (int i = n_orphan - 3; i < n_orphan; i++) {
            orphan_errmsg += "  " + orphan_coordinates[i] + '\n';
          }
        }
        else {
          for (int i = 0; i < 6; i++) {
            orphan_errmsg += "  " + orphan_coordinates[i] + '\n';
          }
        }
        rtWarn("A total of " + std::to_string(orphan_coordinates.size()) + " free coordinate sets "
               "do not correspond to any of the provided topologies, due to atom count "
               "mismatches.\n\nOrphan topologies:\n" + orphan_errmsg, "SystemCache");
      }
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Pull the original list of systems out of the &files namelist
  std::vector<MoleculeSystem> sysvec;
  int nsys = fcon.getSystemDefinitionCount();
  sysvec.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    sysvec.push_back(fcon.getSystem(i));
    const std::string traj_str = sysvec[i].getTrajectoryFileName();
    if (traj_str.size() == 0LLU) {
      sysvec[i].setTrajectoryFileName(fcon.getTrajectoryFileName());
    }
    const std::string chkp_str = sysvec[i].getCheckpointFileName();
    if (chkp_str.size() == 0LLU) {
      sysvec[i].setCheckpointFileName(fcon.getCheckpointFileName());
    }
  }
  
  // Test each coordinate set with respect to topologies of the appropriate size.  Evaluate
  // valence term energy, and take the best scoring result as the indicator of which topology
  // describes which coordinate set.  Make new systems based on each coordinate set.
  const int n_unique_sizes = unique_topology_sizes.size();
  ScoreCard sc(1);
  std::string trajectory_base, trajectory_ext, restart_base, restart_ext;
  splitPath(fcon.getTrajectoryFileName(), &trajectory_base, &trajectory_ext);
  splitPath(fcon.getCheckpointFileName(), &restart_base, &restart_ext);
  std::vector<bool> topology_in_use(n_free_top, false);
  int n_paired_systems = 0;
  for (int i = 0; i < n_unique_sizes; i++) {

    // Loop over all coordinates in this size group.  Try interpreting them with each topology.
    // Store each result as a unique MoleculeSystem and expand the list.
    const int j_llim = unique_coordinate_size_bounds[i];
    const int j_hlim = unique_coordinate_size_bounds[i + 1];
    const int k_llim = unique_topology_size_bounds[i];
    const int k_hlim = unique_topology_size_bounds[i + 1];
    for (int j = j_llim; j < j_hlim; j++) {
      const int icrdj = coordinate_series[j];
      const CoordinateFrameReader cfr(tmp_coordinates_cache[icrdj].data());
      int best_topology;
      double min_bondang_e, min_bond_e, min_angl_e;
      for (int k = k_llim; k < k_hlim; k++) {
        const int tpk = topology_series[k];
        const ValenceKit<double> vk = topology_cache[tpk].getDoublePrecisionValenceKit();
        const double bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
        const double angl_e = evaluateAngleTerms(vk, cfr, &sc, 0);
        const double bondang_e = bond_e + angl_e;
        if (k == k_llim || bondang_e < min_bondang_e) {
          best_topology = tpk;
          min_bondang_e = bondang_e;
          min_bond_e = bond_e;
          min_angl_e = angl_e;
        }
      }

      // Construct the trajectory and restart file names for this system based on generic paths.
      // Hack a solution in the odd event that the user has stuffed their file differentiation
      // behind the final dot.
      const std::string orig_crd_file = tmp_coordinates_cache[icrdj].getFileName();
      std::string orig_base, orig_ext;
      splitPath(getBaseName(orig_crd_file), &orig_base, &orig_ext);
      std::string trajectory_middle("_");
      trajectory_middle += (orig_base.size() > 0) ? orig_base + "." : orig_ext + ".";
      
      // Add this pair to the list of systems
      const std::string syslabel = std::string("PairedSystem") + std::to_string(n_paired_systems);
      sysvec.push_back(MoleculeSystem(topology_cache[best_topology].getFileName(),
                                      tmp_coordinates_cache[icrdj].getFileName(),
                                      trajectory_base + trajectory_middle + trajectory_ext,
                                      restart_base + trajectory_middle + restart_ext, syslabel, 0,
                                      tmp_coordinates_frame_count[icrdj], 1,
                                      tmp_coordinates_kind[icrdj],
                                      fcon.getOutputCoordinateFormat(),
                                      fcon.getCheckpointFormat()));
      n_paired_systems++;

      // Note that the topology is used
      topology_in_use[best_topology] = true;
    }
  }

  // Loop back over systems and make sure that each has a unique restart file (even if these files
  // are not to be written).  Systems with similar trajectory file names must all share the same
  // original topology.
  nsys = sysvec.size();
  bool name_collision;
  do {
    name_collision = false;
    std::vector<std::string> unique_restart_names;
    std::vector<int> name_copies;
    unique_restart_names.reserve(nsys);
    name_copies.reserve(nsys);
    int n_unique = 0;
    for (int i = 0; i < nsys; i++) {
      const int pos = findStringInVector(unique_restart_names, sysvec[i].getCheckpointFileName());
      if (pos != n_unique) {
        name_copies[pos] += 1;
        name_collision = true;
      }
      else {
        unique_restart_names.push_back(sysvec[i].getCheckpointFileName());
        name_copies.push_back(1);
        n_unique++;
      }
    }
    if (name_collision) {
      for (int i = 0; i < n_unique; i++) {
        if (name_copies[i] == 1) {
          continue;
        }
        const std::string overused_name = unique_restart_names[i];
        std::string overused_base, overused_ext;
        splitPath(overused_name, &overused_base, &overused_ext);
        const bool base_substantial = (overused_base.size() > 0LLU);
        const bool ext_substantial = (overused_ext.size() > 0LLU);
        const bool both_substantial = (base_substantial && ext_substantial);
        int overused_count = 0;
        for (int j = 0; j < nsys; j++) {
          if (sysvec[j].getCheckpointFileName() == overused_name) {
            if (both_substantial) {
              sysvec[j].setCheckpointFileName(overused_base + "_" +
                                              std::to_string(overused_count) + "." + overused_ext);
            }
            else if (base_substantial) {
              sysvec[j].setCheckpointFileName(overused_base + "_" +
                                              std::to_string(overused_count));
            }
            else {
              sysvec[j].setCheckpointFileName(overused_ext + "_" +
                                              std::to_string(overused_count));
            }
            overused_count++;
          }
        }
      }
    }
  } while (name_collision);
  do {
    name_collision = false;
    std::vector<std::string> unique_trajectory_names;
    std::vector<std::string> unique_trajectory_basis;
    std::vector<int> divergent_topologies;
    unique_trajectory_names.reserve(nsys);
    unique_trajectory_basis.reserve(nsys);
    divergent_topologies.reserve(nsys);
    int n_unique = 0;
    for (int i = 0; i < nsys; i++) {
      const int pos = findStringInVector(unique_trajectory_names,
                                         sysvec[i].getTrajectoryFileName());
      if (pos != n_unique) {

        // Two systems name the same trajectory, but is that a problem?  Not if the topologies
        // describing each system are the same.
        if (unique_trajectory_basis[pos] != sysvec[i].getTopologyFileName()) {
          divergent_topologies[pos] += 1;
          name_collision = true;
        }
      }
      else {
        unique_trajectory_names.push_back(sysvec[i].getTrajectoryFileName());
        unique_trajectory_basis.push_back(sysvec[i].getTopologyFileName());
        divergent_topologies.push_back(1);
        n_unique++;
      }
    }
    if (name_collision) {
      for (int i = 0; i < n_unique; i++) {
        if (divergent_topologies[i] == 1) {
          continue;
        }
        const std::string overused_name = unique_trajectory_names[i];
        std::string overused_base, overused_ext;
        splitPath(overused_name, &overused_base, &overused_ext);
        const bool base_substantial = (overused_base.size() > 0LLU);
        const bool ext_substantial = (overused_ext.size() > 0LLU);
        const bool both_substantial = (base_substantial && ext_substantial);
        int overused_count = 0;
        std::vector<std::string> new_trajectory_basis;
        std::vector<std::string> new_trajectory_names;
        for (int j = 0; j < nsys; j++) {
          if (sysvec[j].getTrajectoryFileName() != overused_name) {
            continue;
          }
          const std::string jtop_name = sysvec[j].getTopologyFileName();
          if (jtop_name == unique_trajectory_basis[i]) {
            continue;
          }
          bool new_basis_found = false;
          for (int k = 0; k < overused_count; k++) {
            if (jtop_name == new_trajectory_basis[k]) {
              sysvec[j].setTrajectoryFileName(new_trajectory_names[k]);
              new_basis_found = true;
            }
          }
          if (new_basis_found == false) {
            std::string ntraj_name;
            if (both_substantial) {
              ntraj_name = overused_base + "_" + std::to_string(overused_count) + "." +
                           overused_ext;
            }
            else if (base_substantial) {
              ntraj_name = overused_base + "_" + std::to_string(overused_count);
            }
            else {
              ntraj_name = overused_ext + "_" + std::to_string(overused_count);
            }
            sysvec[j].setTrajectoryFileName(ntraj_name);
            new_trajectory_names.push_back(ntraj_name);
            new_trajectory_basis.push_back(jtop_name);
            overused_count++;
          }
        }
      }      
    }
  } while (name_collision);

  // Loop back over the systems (now representing all entries, including the paired free topologies
  // and coordinate sets).  If the topology has already been read, don't read it again.  Read
  // coordinates (and perhaps velocities, if available) into phase space objects.
  std::vector<std::string> current_topology_holdings;
  current_topology_holdings.reserve(n_free_top);
  int cache_track = 0;
  for (int i = 0; i < n_free_top; i++) {
    if (topology_in_use[i]) {
      current_topology_holdings.push_back(topology_cache[cache_track].getFileName());
      cache_track++;
    }
    else {
      topology_cache.erase(topology_cache.begin() + cache_track);
    }
  }
  sdf_cache.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    int top_idx = findStringInVector(current_topology_holdings, sysvec[i].getTopologyFileName());
    bool topology_ok = false;
    if (top_idx >= current_topology_holdings.size()) {
      try {
        topology_cache.push_back(AtomGraph(sysvec[i].getTopologyFileName()));
        current_topology_holdings.push_back(sysvec[i].getTopologyFileName());
        top_idx = topology_cache.size() - 1LLU;
        topology_ok = true;
      }
      catch (std::runtime_error) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The format of topology " + sysvec[i].getTopologyFileName() +
                " for system " + std::to_string(i) + " could not be understood.", "SystemCache");
        case ExceptionResponse::WARN:
          rtWarn("The format of topology " + sysvec[i].getTopologyFileName() +
                 " could not be understood.  System " + std::to_string(i + 1) +
                 " will be skipped.", "SystemCache");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
    }
    else {
      topology_ok = true;
    }
    if (topology_ok) {
      bool coordinates_ok = false;
      const CoordinateFileKind icrd_kind =
        detectCoordinateFileKind(sysvec[i].getInputCoordinateFileName(), "SystemCache");
      switch (icrd_kind) {
      case CoordinateFileKind::AMBER_CRD:
        break;
      case CoordinateFileKind::AMBER_INPCRD:
      case CoordinateFileKind::AMBER_ASCII_RST:
        {
          PhaseSpace next_ps;
          try {
            next_ps.buildFromFile(sysvec[i].getInputCoordinateFileName(), icrd_kind);
            coordinates_ok = true;
          }
          catch (std::runtime_error) {
            switch (policy) {
            case ExceptionResponse::DIE:
              rtErr("The format of coordinate file " + sysvec[i].getInputCoordinateFileName() +
                    " was not readable under the detected " + getEnumerationName(icrd_kind) +
                    " file format.", "SystemCache");
            case ExceptionResponse::WARN:
              rtWarn("The format of coordinate file " + sysvec[i].getInputCoordinateFileName() +
                     " was not readable under the detected " + getEnumerationName(icrd_kind) +
                     " file format.  This coordinate set will be skipped.", "SystemCache");
              break;
            case ExceptionResponse::SILENT:
              break;
            }
          }
          if (coordinates_ok) {
            for (int j = 0; j < sysvec[i].getReplicaCount(); j++) {
              coordinates_cache.push_back(next_ps);
              topology_indices.push_back(top_idx);
              system_input_coordinate_names.push_back(sysvec[i].getInputCoordinateFileName());
              system_trajectory_names.push_back(sysvec[i].getTrajectoryFileName());
              system_checkpoint_names.push_back(sysvec[i].getCheckpointFileName());
              system_labels.push_back(sysvec[i].getLabel());
              system_input_coordinate_kinds.push_back(icrd_kind);
              system_trajectory_kinds.push_back(sysvec[i].getTrajectoryFileKind());
              system_checkpoint_kinds.push_back(sysvec[i].getCheckpointFileKind());
              sdf_cache.emplace_back();
              system_count += 1;
            }
          }
        }
        break;
      case CoordinateFileKind::AMBER_NETCDF:
      case CoordinateFileKind::AMBER_NETCDF_RST:
        break;
      case CoordinateFileKind::SDF:
        {
          // If the user has specified the frame range [ 0, -1 ), take that as a special case
          // indicating that "all frames" shall be read and produce no warning or error.
          // Otherwise, check the specified range against the file's actual contents.
          const int fr_init  = sysvec[i].getStartingFrame();
          const int fr_final = sysvec[i].getFinalFrame();
          std::vector<MdlMol> frame_selection;
          try {
            frame_selection = (fr_init == 0 && fr_final == -1) ?
              readStructureDataFile(sysvec[i].getInputCoordinateFileName()) :
              readStructureDataFile(sysvec[i].getInputCoordinateFileName(), fr_init, fr_final);
            coordinates_ok = true;
          }
          catch (std::runtime_error) {
            switch (policy) {
            case ExceptionResponse::DIE:
              rtErr("The SD file " + sysvec[i].getInputCoordinateFileName() + " was not properly "
                    "read.", "SystemCache");
            case ExceptionResponse::WARN:
              rtWarn("The SD file " + sysvec[i].getInputCoordinateFileName() + " was not properly "
                     "read.  All frames will be skipped.", "SystemCache");
              break;
            case ExceptionResponse::SILENT:
              break;
            }
          }
          if (coordinates_ok) {
            const int jlim = frame_selection.size();
            for (int j = 0; j < jlim; j++) {

              // Apply a simple check to help ensure that the data in the SDF applies to the stated
              // topology.  SD files can contain different molecules, but if it passes this test
              // then there is good reason to trust the user.
              const int top_atom_count = topology_cache[top_idx].getAtomCount();
              if (frame_selection[j].getAtomCount() != top_atom_count) {
                const std::string base_err("The atom count in frame " + std::to_string(j) +
                                           " of SD file " +
                                           getBaseName(sysvec[i].getInputCoordinateFileName()) +
                                           " (" +
                                           std::to_string(frame_selection[j].getAtomCount()) +
                                           ") does not agree with the atom count in topology " +
                                           getBaseName(sysvec[i].getTopologyFileName()) + " (" +
                                           std::to_string(top_atom_count) + ").");
                switch (policy) {
                case ExceptionResponse::DIE:
                  rtErr(base_err, "SystemCache");
                case ExceptionResponse::WARN:
                  rtErr(base_err + "  This frame will be skipped.", "SystemCache");
                  break;
                case ExceptionResponse::SILENT:
                  break;
                }
              }
              for (int k = 0; k < sysvec[i].getReplicaCount(); k++) {
                coordinates_cache.push_back(frame_selection[j].exportPhaseSpace());
                topology_indices.push_back(top_idx);
                system_input_coordinate_names.push_back(sysvec[i].getInputCoordinateFileName());
                system_trajectory_names.push_back(sysvec[i].getTrajectoryFileName());
                system_checkpoint_names.push_back(sysvec[i].getCheckpointFileName());
                system_labels.push_back(sysvec[i].getLabel());
                system_input_coordinate_kinds.push_back(icrd_kind);
                system_trajectory_kinds.push_back(sysvec[i].getTrajectoryFileKind());
                system_checkpoint_kinds.push_back(sysvec[i].getCheckpointFileKind());
                sdf_cache.push_back(frame_selection[j]);
                system_count += 1;
              }
            }
          }
        }
        break;
      case CoordinateFileKind::UNKNOWN:
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The format of coordinate file " + sysvec[i].getTopologyFileName() +
                " for system " + std::to_string(i) + " could not be understood.",
                "SystemCache");
        case ExceptionResponse::WARN:
          rtWarn("The format of coordinate file " + sysvec[i].getTopologyFileName() +
                 " could not be understood.  System " + std::to_string(i + 1) +
                 " will be skipped.", "SystemCache");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        break;
      }
    }
  }
  
  // Collect examples of all systems, taking the first system in the list to use each topology
  // as the example of coordinates for that topology.
  topology_count = topology_cache.size();
  example_indices.resize(topology_count, -1);
  for (int i = 0; i < system_count; i++) {
    const int top_idx = topology_indices[i];
    if (example_indices[top_idx] == -1) {
      example_indices[top_idx] = i;
    }
  }

  // Make a list of all the systems making use of each topology, plus a bounds list to navigate it
  topology_cases.resize(system_count);
  topology_case_bounds.resize(topology_count + 1, 0);
  for (int i = 0; i < system_count; i++) {
    topology_case_bounds[topology_indices[i]] += 1;
  }
  prefixSumInPlace(&topology_case_bounds, PrefixSumType::EXCLUSIVE, "SystemCache");
  for (int i = 0; i < system_count; i++) {
    const int case_idx = topology_case_bounds[topology_indices[i]];
    topology_cases[case_idx] = i;
    topology_case_bounds[topology_indices[i]] = case_idx + 1;
  }
  for (int i = topology_count; i > 0; i--) {
    topology_case_bounds[i] = topology_case_bounds[i - 1];
  }
  topology_case_bounds[0] = 0;

  // With all labels now known, condense a list of all unique labels and mark each system with an
  // index into this list, based on its own label.
  std::vector<bool> label_coverage(system_count, false);
  label_indices.resize(system_count);
  for (int i = 0; i < system_count; i++) {
    if (label_coverage[i]) {
      continue;
    }
    label_coverage[i] = true;
    label_indices[i] = label_count;
    for (int j = i + 1; j < system_count; j++) {
      if (label_coverage[j]) {
        continue;
      }
      if (system_labels[j] == system_labels[i]) {
        label_coverage[j] = true;
        label_indices[j] = label_count;
      }
    }
    label_cache.push_back(system_labels[i]);
    label_count += 1;
  }

  // Count the degeneracy of multiple systems under each label
  label_degeneracy.resize(label_count, 0);
  for (int i = 0; i < system_count; i++) {
    label_degeneracy[label_indices[i]] += 1;
  }
  label_cases.resize(system_count);
  label_case_bounds.resize(label_count + 1);
  for (int i = 0; i < label_count; i++) {
    label_case_bounds[i] = label_degeneracy[i];
  }
  prefixSumInPlace(&label_case_bounds, PrefixSumType::EXCLUSIVE, "Lay out the bounds array for "
                   "label cases in the SystemCache object");
  std::vector<int> label_case_counters = label_case_bounds;
  for (int i = 0; i < system_count; i++) {
    const int ilbl_idx = label_indices[i];
    const int ilcc = label_case_counters[ilbl_idx]; 
    label_cases[ilcc] = i;
    label_case_counters[ilbl_idx] = ilcc + 1;
  }

  // Loop over all systems, computing sub-indices for each degenerate trajectory and checkpoint
  // name should the user require that different systems feed into individual outputs.
  trajectory_subindices.resize(system_count);
  trajectory_degeneracy.resize(system_count);
  checkpoint_subindices.resize(system_count);
  checkpoint_degeneracy.resize(system_count);
  std::vector<bool> traj_coverage(system_count, false);
  std::vector<bool> chkp_coverage(system_count, false);
  std::vector<int> common_name_indices;
  for (int i = 0; i < system_count; i++) {
    if (traj_coverage[i] == false) {
      trajectory_subindices[i] = 0;
      int traj_counter = 1;
      common_name_indices.resize(1);
      common_name_indices[0] = i;
      for (int j = i + 1; j < system_count; j++) {
        if (system_trajectory_names[j] == system_trajectory_names[i]) {
          trajectory_subindices[j] = traj_counter;
          traj_counter++;
          traj_coverage[j] = true;
          common_name_indices.push_back(j);
        }
      }
      for (int j = 0; j < traj_counter; j++) {
        trajectory_degeneracy[common_name_indices[j]] = traj_counter;
      }
    }
    if (chkp_coverage[i] == false) {
      checkpoint_subindices[i] = 0;
      int chkp_counter = 1;
      common_name_indices.resize(1);
      common_name_indices[0] = i;
      for (int j = i + 1; j < system_count; j++) {
        if (system_checkpoint_names[j] == system_checkpoint_names[i]) {
          checkpoint_subindices[j] = chkp_counter;
          chkp_counter++;
          chkp_coverage[j] = true;
          common_name_indices.push_back(j);
        }
      }
      for (int j = 0; j < chkp_counter; j++) {
        checkpoint_degeneracy[common_name_indices[j]] = chkp_counter;
      }
    }
  }
  
  // Create ChemicalFeatures objects to pair with each system--the majority of the work in
  // computing the ChemicalFeatures lies in ring detection and drawing the Lewis structure to
  // determine formal charges and bond orders, aspects which are invariant to coordinates.
  // Compute the ChemicalFeatures for each topology and then reassess the chiralities for each
  // system.
  features_cache.reserve(system_count);
  for (int i = 0; i < system_count; i++) {
    
    // Determine whether chemical features have been computed for this topology already.
    // If not, create them anew.  If so, copy the features and recompute chiral orientations.
    const int top_idx = topology_indices[i];
    if (i > example_indices[top_idx]) {
      features_cache.emplace_back(features_cache[example_indices[top_idx]]);
    }
    else {
      features_cache.emplace_back(&topology_cache[top_idx],
                                  CoordinateFrameReader(coordinates_cache[i]), map_chemfe_rotators,
                                  300.0, timer_in);
    }
    features_cache[i].findChiralOrientations(CoordinateFrameReader(coordinates_cache[i]));
  }

  // Create Exclusion masks of the appropriate kind
  std::vector<bool> need_static_mask(topology_count, false);
  std::vector<bool> need_forward_mask(topology_count, false);
  for (int i = 0; i < system_count; i++) {
    const int top_idx = topology_indices[i];
    switch (coordinates_cache[i].getUnitCellType()) {
    case UnitCellType::NONE:
      need_static_mask[top_idx] = true;
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      need_forward_mask[top_idx] = true;
      break;
    }
  }
  static_masks_cache.reserve(topology_count);
  forward_masks_cache.reserve(topology_count);
  for (int i = 0; i < topology_count; i++) {
    if (need_static_mask[i]) {
      static_masks_cache.emplace_back(&topology_cache[i]);
    }
    else {
      static_masks_cache.emplace_back(nullptr);
    }
    if (need_forward_mask[i]) {
      forward_masks_cache.emplace_back(&topology_cache[i]);
    }
    else {
      forward_masks_cache.emplace_back(nullptr);
    }
  }

  // Use the chemical features objects to make restraint apparatuses associated with the various
  // labels found in &restraint namelists.
  int nrst_labels = 0;
  const int nrst_nml = rstcon.size();
  std::vector<std::string> rst_group_labels;
  std::vector<std::vector<int>> rst_group_namelists;
  std::vector<bool> rst_label_group_assigned(nrst_nml, false);
  std::vector<std::string> all_nml_labels(nrst_nml);
  for (int i = 0; i < nrst_nml; i++) {
    all_nml_labels[i] = rstcon[i].getSystemLabel();
  }
  for (int i = 0; i < nrst_nml; i++) {
    if (rst_label_group_assigned[i]) {
      continue;
    }
    rst_label_group_assigned[i] = true;
    const std::string ilabel = rstcon[i].getSystemLabel();
    rst_group_labels.push_back(ilabel);
    std::vector<int> tmp_rst_list(1, i);
    for (int j = i + 1; j < nrst_nml; j++) {
      if (ilabel == all_nml_labels[j]) {
        tmp_rst_list.push_back(j);
        rst_label_group_assigned[j] = true;
      }
    }
    rst_group_namelists.push_back(tmp_rst_list);
  }
  const int nrst_groups = rst_group_labels.size();
  
  // Loop over all systems and apply restraints as stated in the labeled groups
  std::vector<int> blank_ra(topology_cache.size(), -1);
  restraint_indices.resize(system_count);
  for (int i = 0; i < system_count; i++) {
    RestraintApparatus ra(&topology_cache[topology_indices[i]]);
    const AtomGraph *iag_ptr = &topology_cache[topology_indices[i]];
    const CoordinateFrameReader cfr(coordinates_cache[i]);
    for (int j = 0; j < nrst_groups; j++) {
      if (rst_group_labels[j] != system_labels[i]) {
        continue;
      }
      const int nrst_applicable = rst_group_namelists[j].size();
      const int m = topology_indices[i];
      for (int k = 0; k < nrst_applicable; k++) {
        ra.addRestraints(rstcon[rst_group_namelists[j][k]].getRestraint(iag_ptr, features_cache[m],
                                                                        cfr));
      }
    }

    // If there are no restraints, log this apparatus as a new case only if it is the first such
    // blank apparatus for the system topology at hand.
    if (ra.getTotalRestraintCount() == 0) {
      if (blank_ra[topology_indices[i]] < 0) {
        restraints_cache.push_back(ra);
        restraint_count += 1;
        blank_ra[topology_indices[i]] = static_cast<int>(restraints_cache.size()) - 1;
      }
      restraint_indices[i] = blank_ra[topology_indices[i]];
    }
    else {
      restraints_cache.push_back(ra);
      restraint_indices[i] = static_cast<int>(restraints_cache.size()) - 1;
    }
  }
}

//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache(const FilesControls &fcon, const ExceptionResponse policy_in,
                         const MapRotatableGroups map_chemfe_rotators,
                         const PrintSituation expectation_in, StopWatch *timer_in) :
    SystemCache(fcon, {}, policy_in, map_chemfe_rotators, expectation_in, timer_in)
{}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCount() const {
  return topology_count;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getLabelCount() const {
  return label_count;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getRestraintCount() const {
  return restraint_count;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemTopologyIndex(const int coord_index) const {
  checkSystemBounds(coord_index, "getSystemTopologyIndex");
  return topology_indices[coord_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getSystemTopologyIndex() const {
  std::vector<int> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = topology_indices[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemExampleIndex(const int topology_index) const {
  return example_indices[topology_index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* SystemCache::getTopologyPointer(const int topology_index) const {
  checkTopologyBounds(topology_index, "getTopologyPointer");
  return &topology_cache[topology_index];
}

//-------------------------------------------------------------------------------------------------
AtomGraph* SystemCache::getTopologyPointer(const int topology_index) {
  checkTopologyBounds(topology_index, "getTopologyPointer");
  return &topology_cache[topology_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> SystemCache::getTopologyPointer() const {
  const size_t ntop = topology_cache.size();
  std::vector<AtomGraph*> result(ntop);
  for (size_t i = 0; i < ntop; i++) {
    result[i] = const_cast<AtomGraph*>(&topology_cache[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph*> SystemCache::getTopologyPointer() {
  const size_t ntop = topology_cache.size();
  std::vector<AtomGraph*> result(ntop);
  for (size_t i = 0; i < ntop; i++) {
    result[i] = &topology_cache[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& SystemCache::getTopology(const int topology_index) const {
  checkTopologyBounds(topology_index, "getTopologyPointer");
  return topology_cache[topology_index];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCacheIndex(const AtomGraph *ag) const {

  // Try matching the topology to an object in memory by its pointer.
  for (int i = 0; i < topology_count; i++) {
    if (ag == &topology_cache[i]) {
      return i;
    }
  }

  // Try matching the topology by its original file name.
  const std::string query_topname = ag->getFileName();
  for (int i = 0; i < topology_count; i++) {
    if (topology_cache[i].getFileName() == query_topname) {
      return i;
    }
  }

  // Return the failed result.
  return topology_count;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCacheIndex(const AtomGraph &ag) const {
  return getTopologyCacheIndex(ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*>
SystemCache::getTopologiesMatchingLabel(const std::string &query_label) const {
  const std::vector<int> matches = getMatchingSystemIndices(query_label);
  const size_t nmatch = matches.size();
  std::vector<bool> is_included(topology_count, false);
  for (size_t i = 0; i < nmatch; i++) {
    is_included[topology_indices[matches[i]]] = true;
  }
  std::vector<AtomGraph*> result;
  result.reserve(sum<int>(is_included));
  for (int i = 0; i < topology_count; i++) {
    if (is_included[i]) {
      result.push_back(const_cast<AtomGraph*>(&topology_cache[i]));
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph*> SystemCache::getTopologiesMatchingLabel(const std::string &query_label) {
  const std::vector<int> matches = getMatchingSystemIndices(query_label);
  const size_t nmatch = matches.size();
  std::vector<bool> is_included(topology_count, false);
  for (size_t i = 0; i < nmatch; i++) {
    is_included[topology_indices[matches[i]]] = true;
  }
  std::vector<AtomGraph*> result;
  result.reserve(sum<int>(is_included));
  for (int i = 0; i < topology_count; i++) {
    if (is_included[i]) {
      result.push_back(&topology_cache[i]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getLabel(const int label_index) const {
  checkLabelBounds(label_index, "getLabel");
  return label_cache[label_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> SystemCache::getLabelsMatchingTopology(const AtomGraph *query_ag) const {
  const std::vector<int> matches = getMatchingSystemIndices(query_ag);
  const size_t nmatch = matches.size();
  std::vector<bool> is_included(label_count, false);
  for (size_t i = 0; i < nmatch; i++) {
    is_included[label_indices[matches[i]]] = true;
  }
  std::vector<std::string> result;
  result.reserve(sum<int>(is_included));
  for (int i = 0; i < label_count; i++) {
    if (is_included[i]) {
      result.push_back(label_cache[i]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> SystemCache::getLabelsMatchingTopology(const AtomGraph &query_ag) const {
  return getLabelsMatchingTopology(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getLabelCacheIndex(const std::string &query) const {
  for (int i = 0; i < label_count; i++) {
    if (label_cache[i] == query) {
      return i;
    }
  }
  return label_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getMatchingSystemIndices(const AtomGraph* query_ag,
                                                       const std::string &query_label) const {

  // Match the topology by its pointer
  const int top_idx = getTopologyCacheIndex(query_ag);

  // Raise an exception if the topology still could not be matched.
  if (top_idx == topology_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A topology originating from file " + getBaseName(query_ag->getFileName()) +
            " could not be matched to any in the cache of " + std::to_string(topology_count) +
            " topologies.", "SystemCache", "getMatchingSystemIndices");
    case ExceptionResponse::WARN:
      rtWarn("A topology originating from file " + getBaseName(query_ag->getFileName()) +
            " could not be matched to any in the cache of " + std::to_string(topology_count) +
            " topologies.  No result will be returned.", "SystemCache",
             "getMatchingSystemIndices");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  
  // Find all matching cases.
  std::vector<int> result;
  for (int i = topology_case_bounds[top_idx]; i < topology_case_bounds[top_idx + 1]; i++) {
    if (system_labels[topology_cases[i]] == query_label) {
      result.push_back(topology_cases[i]);
    }
  }
  if (result.size() == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No topology originating in file " + getBaseName(query_ag->getFileName()) + " could "
            "be matched to label " + query_label + ".", "SystemCache", "getMatchingSystemIndices");
    case ExceptionResponse::WARN:
      rtWarn("No topology originating in file " + getBaseName(query_ag->getFileName()) + " could "
             "be matched to label " + query_label + ".", "SystemCache",
             "getMatchingSystemIndices");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getMatchingSystemIndices(const AtomGraph& query_ag,
                                                       const std::string &query_label) const {
  return getMatchingSystemIndices(query_ag.getSelfPointer(), query_label);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getMatchingSystemIndices(const AtomGraph* query_ag) const {
  const int top_idx = getTopologyCacheIndex(query_ag);
  const int llim = topology_case_bounds[top_idx];
  const int hlim = topology_case_bounds[top_idx + 1];
  std::vector<int> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = topology_cases[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getMatchingSystemIndices(const AtomGraph& query_ag) const {
  return getMatchingSystemIndices(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getMatchingSystemIndices(const std::string &query_label) const {
  const int label_idx = getLabelCacheIndex(query_label);
  const int llim = label_case_bounds[label_idx];
  const int hlim = label_case_bounds[label_idx + 1];
  std::vector<int> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = label_cases[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getFirstMatchingSystemIndex(const AtomGraph *query_ag) const {
  const std::vector<int> sys_idx = getMatchingSystemIndices(query_ag);
  if (sys_idx.size() == 0) {
    rtErr("No system was found matching the topology originating in file " +
          getBaseName(query_ag->getFileName()) + ".", "SystemCache", "getInputCoordinatesKind");
  }
  return sys_idx[0];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getFirstMatchingSystemIndex(const AtomGraph *query_ag,
                                             const std::string &query_label) const {
  const std::vector<int> sys_idx = getMatchingSystemIndices(query_ag, query_label);
  if (sys_idx.size() == 0) {
    rtErr("No system was found matching the topology originating in file " +
          getBaseName(query_ag->getFileName()) + " and label \"" + query_label + "\".",
          "SystemCache", "getInputCoordinatesKind");
  }
  return sys_idx[0];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getFirstMatchingSystemIndex(const std::string &query_label) const {
  const std::vector<int> sys_idx = getMatchingSystemIndices(query_label);
  if (sys_idx.size() == 0) {
    rtErr("No system was found matching label \"" + query_label + "\".", "SystemCache",
          "getInputCoordinatesKind");
  }
  return sys_idx[0];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* SystemCache::getSystemTopologyPointer(const int index) const {
  checkSystemBounds(index, "getSystemTopologyPointer");
  const AtomGraph* tp_cache_ptr = topology_cache.data();
  return &tp_cache_ptr[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
AtomGraph* SystemCache::getSystemTopologyPointer(const int index) {
  checkSystemBounds(index, "getSystemTopologyPointer");
  AtomGraph* tp_cache_ptr = topology_cache.data();
  return &tp_cache_ptr[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> SystemCache::getSystemTopologyPointer() const {
  const size_t nsys = coordinates_cache.size();
  std::vector<AtomGraph*> result(nsys);
  for (size_t i = 0; i < nsys; i++) {
    result[i] = const_cast<AtomGraph*>(&topology_cache[topology_indices[i]]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph*> SystemCache::getSystemTopologyPointer() {
  const size_t nsys = coordinates_cache.size();
  std::vector<AtomGraph*> result(nsys);
  AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &tp_data[topology_indices[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& SystemCache::getSystemTopology(const int index) const {
  checkSystemBounds(index, "getSystemTopology");
  return topology_cache[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis SystemCache::exportTopologySynthesis(const GpuDetails &gpu,
                                                        const ExceptionResponse policy) const {
  return AtomGraphSynthesis(this->getTopologyPointer(), this->getRestraintPointer(),
                            topology_indices, incrementingSeries<int>(0, system_count), policy,
                            gpu);
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace* SystemCache::getCoordinatePointer(const int index) const {
  checkSystemBounds(index, "getCoordinatePointer");
  const PhaseSpace* crd_cache_ptr = coordinates_cache.data();
  return &crd_cache_ptr[index];
}

//-------------------------------------------------------------------------------------------------
PhaseSpace* SystemCache::getCoordinatePointer(const int index) {
  checkSystemBounds(index, "getCoordinatePointer");
  PhaseSpace* crd_cache_ptr = coordinates_cache.data();
  return &crd_cache_ptr[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<PhaseSpace*> SystemCache::getCoordinatePointer() const {
  const PhaseSpace* ps_ptr = coordinates_cache.data();
  const int nsys = coordinates_cache.size();
  std::vector<PhaseSpace*> result(nsys);
  for (size_t i = 0; i < nsys; i++) {
    result[i] = const_cast<PhaseSpace*>(&ps_ptr[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<PhaseSpace*> SystemCache::getCoordinatePointer() {
  PhaseSpace* ps_ptr = coordinates_cache.data();
  const int nsys = coordinates_cache.size();
  std::vector<PhaseSpace*> result(nsys);
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &ps_ptr[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace& SystemCache::getCoordinates(const int index) const {
  checkSystemBounds(index, "getCoordinates");
  return coordinates_cache[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<PhaseSpace>& SystemCache::getCoordinates() const {
  return coordinates_cache;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis SystemCache::exportCoordinateSynthesis(const int globalpos_scale_bits,
                                                           const int velocity_scale_bits,
                                                           const int force_scale_bits) const {
  return PhaseSpaceSynthesis(coordinates_cache, this->getSystemTopologyPointer(),
                             globalpos_scale_bits, 26, velocity_scale_bits, force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
const MdlMol* SystemCache::getStructureDataEntryPointer(int index) const {
  checkSystemBounds(index, "getStructureDataPointer");
  return &sdf_cache[index];
}

//-------------------------------------------------------------------------------------------------
MdlMol* SystemCache::getStructureDataEntryPointer(int index) {
  checkSystemBounds(index, "getStructureDataPointer");
  return &sdf_cache[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<MdlMol*> SystemCache::getStructureDataEntryPointer() const {
  std::vector<MdlMol*> result;
  result.reserve(system_count);
  for (int i = 0; i < system_count; i++) {
    result.emplace_back(const_cast<MdlMol*>(&sdf_cache[i]));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMol*> SystemCache::getStructureDataEntryPointer() {
  std::vector<MdlMol*> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = &sdf_cache[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const MdlMol& SystemCache::getStructureDataEntry(int index) const {
  checkSystemBounds(index, "getStructureDataEntry");
  return sdf_cache[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<MdlMol>& SystemCache::getStructureDataEntry() const {
  return sdf_cache;
}

//-------------------------------------------------------------------------------------------------
const ChemicalFeatures* SystemCache::getFeaturesPointer(const int index) const {
  checkSystemBounds(index, "getFeaturesPointer");
  return &features_cache[index];
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures* SystemCache::getFeaturesPointer(const int index) {
  checkSystemBounds(index, "getFeaturesPointer");
  return &features_cache[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<ChemicalFeatures*> SystemCache::getFeaturesPointer() const {
  std::vector<ChemicalFeatures*> result(system_count);
  for (size_t i = 0; i < system_count; i++) {
    result[i] = const_cast<ChemicalFeatures*>(&features_cache[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ChemicalFeatures*> SystemCache::getFeaturesPointer() {
  std::vector<ChemicalFeatures*> result(system_count);
  for (size_t i = 0; i < system_count; i++) {
    result[i] = &features_cache[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const ChemicalFeatures& SystemCache::getFeatures(const int index) const {
  checkSystemBounds(index, "getFeaturesReference");
  return features_cache[index];
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus* SystemCache::getRestraintPointer(const int index) const {
  checkSystemBounds(index, "getRestraintPointer");
  return &restraints_cache[restraint_indices[index]];
}

//-------------------------------------------------------------------------------------------------
RestraintApparatus* SystemCache::getRestraintPointer(const int index) {
  checkSystemBounds(index, "getRestraintPointer");
  return &restraints_cache[restraint_indices[index]];
}

//-------------------------------------------------------------------------------------------------
const std::vector<RestraintApparatus*> SystemCache::getRestraintPointer() const {
  std::vector<RestraintApparatus*> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = const_cast<RestraintApparatus*>(&restraints_cache[restraint_indices[i]]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<RestraintApparatus*> SystemCache::getRestraintPointer() {
  std::vector<RestraintApparatus*> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = &restraints_cache[restraint_indices[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus& SystemCache::getRestraints(const int index) const {
  checkSystemBounds(index, "getRestraints");
  return restraints_cache[restraint_indices[index]];
}

//-------------------------------------------------------------------------------------------------
const std::vector<StaticExclusionMask*> SystemCache::getUniqueStaticMaskPointers() const {
  std::vector<StaticExclusionMask*> result(topology_count);
  for (int i = 0; i < topology_count; i++) {
    result[i] = const_cast<StaticExclusionMask*>(&static_masks_cache[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<StaticExclusionMask*> SystemCache::getUniqueStaticMaskPointers() {
  std::vector<StaticExclusionMask*> result(topology_count);
  for (int i = 0; i < topology_count; i++) {
    result[i] = &static_masks_cache[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const StaticExclusionMask* SystemCache::getSystemStaticMaskPointer(const int index) const {
  checkSystemBounds(index, "getSystemStaticMaskPointer");
  const int top_idx = topology_indices[index];
  return &static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMask* SystemCache::getSystemStaticMaskPointer(const int index) {
  checkSystemBounds(index, "getSystemStaticMaskPointer");
  const int top_idx = topology_indices[index];
  return &static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
const StaticExclusionMask& SystemCache::getSystemStaticMask(const int index) const {
  checkSystemBounds(index, "getSystemStaticMask");
  const int top_idx = topology_indices[index];
  return static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
const ForwardExclusionMask* SystemCache::getSystemForwardMaskPointer(const int index) const {
  checkSystemBounds(index, "getSystemForwardMaskPointer");
  const int top_idx = topology_indices[index];
  return &forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
ForwardExclusionMask* SystemCache::getSystemForwardMaskPointer(const int index) {
  checkSystemBounds(index, "getSystemForwardMaskPointer");
  const int top_idx = topology_indices[index];
  return &forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
const ForwardExclusionMask& SystemCache::getSystemForwardMask(const int index) const {
  checkSystemBounds(index, "getSystemForwardMask");
  const int top_idx = topology_indices[index];
  return forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCaseCount(const int topology_index) const {
  return topology_case_bounds[topology_index + 1] - topology_case_bounds[topology_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getTopologicalCases(const int topology_index) const {
  const int llim = topology_case_bounds[topology_index];
  const int hlim = topology_case_bounds[topology_index + 1];
  std::vector<int> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = topology_cases[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getInputCoordinatesName(const int system_index) const {
  checkSystemBounds(system_index, "getInputCoordinatesName");
  return system_input_coordinate_names[system_index];
}

//-------------------------------------------------------------------------------------------------
std::string SystemCache::getTrajectoryName(const int system_index) const {
  checkSystemBounds(system_index, "getTrajectoryName");
  return nondegenerateName(system_trajectory_names[system_index], CoordinateFileRole::TRAJECTORY,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getTrajectoryName(const AtomGraph *query_ag) const {
  return system_trajectory_names[getFirstMatchingSystemIndex(query_ag)];
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getTrajectoryName(const AtomGraph &query_ag) const {
  return getTrajectoryName(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getTrajectoryName(const AtomGraph *query_ag,
                                                  const std::string &query_label) const {
  return system_trajectory_names[getFirstMatchingSystemIndex(query_ag, query_label)];
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getTrajectoryName(const AtomGraph &query_ag,
                                                  const std::string &query_label) const {
  return getTrajectoryName(query_ag.getSelfPointer(), query_label);
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getTrajectoryName(const std::string &query_label) const {
  return system_trajectory_names[getFirstMatchingSystemIndex(query_label)];
}

//-------------------------------------------------------------------------------------------------
std::string SystemCache::getCheckpointName(const int system_index) const {
  checkSystemBounds(system_index, "getCheckpointName");
  return nondegenerateName(system_checkpoint_names[system_index], CoordinateFileRole::CHECKPOINT,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getCheckpointName(const AtomGraph *query_ag) const {
  return system_checkpoint_names[getFirstMatchingSystemIndex(query_ag)];
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getCheckpointName(const AtomGraph &query_ag) const {
  return getCheckpointName(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getCheckpointName(const AtomGraph *query_ag,
                                                  const std::string &query_label) const {
  return system_checkpoint_names[getFirstMatchingSystemIndex(query_ag, query_label)];
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getCheckpointName(const AtomGraph &query_ag,
                                                  const std::string &query_label) const {
  return getCheckpointName(query_ag.getSelfPointer(), query_label);
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getCheckpointName(const std::string &query_label) const {
  return system_checkpoint_names[getFirstMatchingSystemIndex(query_label)];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemLabelIndex(int system_index) const {
  return label_indices[system_index];
}

//-------------------------------------------------------------------------------------------------
const std::string& SystemCache::getSystemLabel(const int system_index) const {
  checkSystemBounds(system_index, "getSystemLabel");
  return system_labels[system_index];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getInputCoordinatesKind(const int system_index) const {
  checkSystemBounds(system_index, "getInputCoordinatesKind");
  return system_input_coordinate_kinds[system_index];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getInputCoordinatesKind(const AtomGraph *query_ag) const {
  return system_input_coordinate_kinds[getFirstMatchingSystemIndex(query_ag)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getInputCoordinatesKind(const AtomGraph &query_ag) const {
  return getInputCoordinatesKind(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getInputCoordinatesKind(const AtomGraph *query_ag,
                                                        const std::string &query_label) const {
  return system_input_coordinate_kinds[getFirstMatchingSystemIndex(query_ag, query_label)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getInputCoordinatesKind(const AtomGraph &query_ag,
                                                        const std::string &query_label) const {
  return getInputCoordinatesKind(query_ag.getSelfPointer(), query_label);
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getInputCoordinatesKind(const std::string &query_label) const {
  return system_input_coordinate_kinds[getFirstMatchingSystemIndex(query_label)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getTrajectoryKind(const int system_index) const {
  checkSystemBounds(system_index, "getTrajectoryKind");
  return system_trajectory_kinds[system_index];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getTrajectoryKind(const AtomGraph *query_ag) const {
  return system_trajectory_kinds[getFirstMatchingSystemIndex(query_ag)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getTrajectoryKind(const AtomGraph &query_ag) const {
  return getTrajectoryKind(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getTrajectoryKind(const AtomGraph *query_ag,
                                                  const std::string &query_label) const {
  return system_trajectory_kinds[getFirstMatchingSystemIndex(query_ag, query_label)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getTrajectoryKind(const AtomGraph &query_ag,
                                                  const std::string &query_label) const {
  return getTrajectoryKind(query_ag.getSelfPointer(), query_label);
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getTrajectoryKind(const std::string &query_label) const {
  return system_trajectory_kinds[getFirstMatchingSystemIndex(query_label)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getCheckpointKind(const int system_index) const {
  checkSystemBounds(system_index, "getCheckpointKind");
  return system_checkpoint_kinds[system_index];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getCheckpointKind(const AtomGraph *query_ag) const {
  return system_checkpoint_kinds[getFirstMatchingSystemIndex(query_ag)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getCheckpointKind(const AtomGraph &query_ag) const {
  return getCheckpointKind(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getCheckpointKind(const AtomGraph *query_ag,
                                                  const std::string &query_label) const {
  return system_checkpoint_kinds[getFirstMatchingSystemIndex(query_ag, query_label)];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getCheckpointKind(const AtomGraph &query_ag,
                                                  const std::string &query_label) const {
  return getCheckpointKind(query_ag.getSelfPointer(), query_label);
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getCheckpointKind(const std::string &query_label) const {
  return system_checkpoint_kinds[getFirstMatchingSystemIndex(query_label)];
}

//-------------------------------------------------------------------------------------------------
PrintSituation SystemCache::getPrintingProtocol() const {
  return expectation;
}

//-------------------------------------------------------------------------------------------------
PrintSituation SystemCache::getPrintingProtocol(const CoordinateFileRole purpose,
                                                const int system_index) const {
  if (determineFileMerger(purpose, system_index)) {

    // If the system is the first of those grouped under a particular label, return the original
    // expectation.  Otherwise, set the printing to append to the growing file.
    switch (purpose) {
    case CoordinateFileRole::INITIATE:
      rtErr("No printing protocol is available for input coordinates.", "SystemCache",
            "getPrintingProtocol");
    case CoordinateFileRole::TRAJECTORY:
      return (trajectory_subindices[system_index] > 0) ? PrintSituation::APPEND : expectation;
    case CoordinateFileRole::CHECKPOINT:
      return (checkpoint_subindices[system_index] > 0) ? PrintSituation::APPEND : expectation;
    }
  }
  else {
    return expectation;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SystemCache* SystemCache::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
void SystemCache::setPrintingProtocol(const PrintSituation expectation_in) {
  expectation = expectation_in;
}

//-------------------------------------------------------------------------------------------------
TrajectoryFusion SystemCache::getTrajectoryFusionProtocol() const {
  return file_merger_protocol;
}

//-------------------------------------------------------------------------------------------------
void SystemCache::setTrajectoryFusionProtocol(const TrajectoryFusion file_merger_protocol_in) {
  file_merger_protocol = file_merger_protocol_in;
}

//-------------------------------------------------------------------------------------------------
void SystemCache::checkSystemBounds(const int index, const char* caller) const {
  if (index < 0 || index >= system_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(system_count) + ".", "SystemCache", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void SystemCache::checkTopologyBounds(const int index, const char* caller) const {
  if (index < 0 || index >= topology_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(topology_count) + ".", "SystemCache", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void SystemCache::checkLabelBounds(const int index, const char* caller) const {
  if (index < 0 || index >= label_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(label_count) + ".", "SystemCache", caller);
  }
}

//-------------------------------------------------------------------------------------------------
bool SystemCache::determineFileMerger(const CoordinateFileRole purpose,
                                      const int system_index) const {
  switch (file_merger_protocol) {
  case TrajectoryFusion::ON:
    switch (purpose) {
    case CoordinateFileRole::INITIATE:
      return false;
    case CoordinateFileRole::TRAJECTORY:

      // File formats that can be trajectories are fit to accept multiple frames per file.
      return true;
    case CoordinateFileRole::CHECKPOINT:
      switch (system_checkpoint_kinds[system_index]) {
      case CoordinateFileKind::AMBER_CRD:
      case CoordinateFileKind::AMBER_NETCDF:
      case CoordinateFileKind::SDF:
        return true;
      case CoordinateFileKind::AMBER_INPCRD:
      case CoordinateFileKind::AMBER_ASCII_RST:
      case CoordinateFileKind::AMBER_NETCDF_RST:    

        // Files will still not be concatenated if the format can only support a single frame.
        return false;
      case CoordinateFileKind::UNKNOWN:
        rtErr("The format of an output file in the role of " + getEnumerationName(purpose) +
              " must be known.", "SystemCache", "determineFileMerger");
      }
      break;
    }
    break;
  case TrajectoryFusion::OFF:
    return false;
  case TrajectoryFusion::AUTO:
    switch (purpose) {
    case CoordinateFileRole::INITIATE:
      return false;
    case CoordinateFileRole::TRAJECTORY:
      return false;
    case CoordinateFileRole::CHECKPOINT:
      switch (system_checkpoint_kinds[system_index]) {
      case CoordinateFileKind::AMBER_CRD:
      case CoordinateFileKind::AMBER_NETCDF:
      case CoordinateFileKind::SDF:
        return true;
      case CoordinateFileKind::AMBER_INPCRD:
      case CoordinateFileKind::AMBER_ASCII_RST:
      case CoordinateFileKind::AMBER_NETCDF_RST:

        // Files will not be concatenated if the format can only support a single frame.
        return false;
      case CoordinateFileKind::UNKNOWN:
        rtErr("The format of an output file in the role of " + getEnumerationName(purpose) +
              " must be known.", "SystemCache", "determineFileMerger");
      }
      break;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string SystemCache::nondegenerateName(const std::string &fname_in,
                                           const CoordinateFileRole purpose,
                                           const int system_index) const {

  // The system does not need a non-degenerate name if there is only one system under the label,
  // if the systems use the same topology and the file has a trajectory format, or if the output
  // is a checkpoint and the output format is suitable for multiple frames.
  int outfile_idx, outfile_deg;
  switch (purpose) {
  case CoordinateFileRole::INITIATE:
    return fname_in;
  case CoordinateFileRole::TRAJECTORY:
    outfile_idx = trajectory_subindices[system_index];
    outfile_deg = trajectory_degeneracy[system_index];
    break;
  case CoordinateFileRole::CHECKPOINT:
    outfile_idx = checkpoint_subindices[system_index];
    outfile_deg = checkpoint_degeneracy[system_index];
    break;
  }
  std::string result;
  if (outfile_deg > 1 && determineFileMerger(purpose, system_index) == false) {
    int last_slash = 0;
    int last_dot = 0;
    const char osc = osSeparator();
    for (int i = 0; i < fname_in.size(); i++) {

      // This logic assumes that the OS directory separator is not a '.', but no file system works
      // in that way.
      if (fname_in[i] == osc) {
        last_slash = i;
      }
      else if (fname_in[i] == '.') {
        last_dot = i;
      }
    }
    if (last_dot > last_slash) {
      result = fname_in.substr(0, last_dot) + "_" + std::to_string(outfile_idx) +
               fname_in.substr(last_dot, fname_in.size());
    }
    else {
      result = fname_in + "_" + std::to_string(outfile_idx);
    }
  }
  else {
    result = fname_in;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMaskSynthesis createMaskSynthesis(const SystemCache &sc) {

  // Loop over all unique topologies in the cache.  Check for isolated boundary conditions.
  const int ntop = sc.getTopologyCount();
  bool has_isolated_bc = false;
  bool has_periodic_bc = false;
  for (int i = 0; i < ntop; i++) {
    switch (sc.getTopologyPointer(i)->getUnitCellType()) {
    case UnitCellType::NONE:
      has_isolated_bc = true;
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      has_periodic_bc = true;
      break;
    }
  }
  if (has_isolated_bc && has_periodic_bc == false) {
    const int nsys = sc.getSystemCount();
    std::vector<StaticExclusionMask*> sev(ntop);
    std::vector<int> tp_idx(nsys);
    for (int i = 0; i < sc.getTopologyCount(); i++) {
      const StaticExclusionMask* sei = sc.getSystemStaticMaskPointer(sc.getSystemExampleIndex(i));
      sev[i] = const_cast<StaticExclusionMask*>(sei);
    }
    for (int i = 0; i < nsys; i++) {
      tp_idx[i] = sc.getSystemTopologyIndex(i);
    }
    StaticExclusionMaskSynthesis result(sev, tp_idx);
    return result;
  }
  else if (has_isolated_bc && has_periodic_bc) {
    rtErr("The systems provided contain a mixture of periodic and non-periodic boundary "
          "conditions.  These systems cannot all be compiled together into a synthesis for "
          "simulations.", "createMaskSynthesis");
  }
  else {
    StaticExclusionMaskSynthesis result;
    return result;
  }
  __builtin_unreachable();
}

} // namespace synthesis
} // namespace stormm
