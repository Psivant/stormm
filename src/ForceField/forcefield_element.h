// -*-c++-*-
#ifndef STORMM_FORCEFIELD_ELEMENT_H
#define STORMM_FORCEFIELD_ELEMENT_H

#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "forcefield_enumerators.h"

namespace stormm {
namespace modeling {

using constants::ExceptionResponse;
using topology::AtomGraph;
using topology::TorsionKind;
using topology::VirtualSiteKind;

/// \brief A versatile object for collecting the parameters and scope of applicability of any
///        molecular mechanics force field term.  This relies on enumerations to inform whether
///        the term applies to atom types or atom and residue names and the nature of the term.
class ForceFieldElement {
public:

  /// \brief A variety of constructors can load one or more atoms and their properties.
  ///
  /// Overloaded:
  ///   - Create an empty element with only a ParameterKind (default NONE)
  ///   - Load non-bonded terms describing a specific atom
  ///   - Load valence terms describing a bond, angle, dihedral, Urey-Bradley, or CHARMM improper
  ///   - Load an entire CMAP surface
  ///   - Load a virtual site with frame specifications
  ///
  /// \param kind_in        The kind of force field parameter, obligatory for every constructor
  /// \param atom_i_in      Name or type of atom I in the term
  /// \param atom_j_in      Name or type of atom J in the term
  /// \param atom_k_in      Name or type of atom K in the term
  /// \param atom_l_in      Name or type of atom L in the term
  /// \param atom_m_in      Name or type of atom M in the term
  /// \param resi_i_in      Name or type of residue I in the term
  /// \param resi_j_in      Name or type of residue J in the term
  /// \param resi_k_in      Name or type of residue K in the term
  /// \param resi_l_in      Name or type of residue L in the term
  /// \param resi_m_in      Name or type of residue M in the term
  /// \param surface_in     Surface values for an entire CMAP
  /// \param frame_type_in  Frame type of the virtual site, if that is what this object contains
  /// \{
  ForceFieldElement(ParameterKind kind_in = ParameterKind::NONE);
  
  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, char4 atom_k_in);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, char4 atom_k_in,
                    char4 atom_l_in, TorsionKind tkind_in = TorsionKind::PROPER);

  ForceFieldElement(ParameterKind kind_in, VirtualSiteKind frame_type_in, char4 atom_i_in,
                    char4 atom_j_in, char4 atom_k_in, char4 residue_i_in, char4 residue_j_in,
                    char4 residue_k_in);

  ForceFieldElement(ParameterKind kind_in, VirtualSiteKind frame_type_in, char4 atom_i_in,
                    char4 atom_j_in, char4 atom_k_in, char4 atom_l_in, char4 residue_i_in,
                    char4 residue_j_in, char4 residue_k_in, char4 residue_l_in);

  ForceFieldElement(ParameterKind kind_in, VirtualSiteKind frame_type_in, char4 atom_i_in,
                    char4 atom_j_in, char4 atom_k_in, char4 atom_l_in, char4 atom_m_in,
                    char4 residue_i_in, char4 residue_j_in, char4 residue_k_in, char4 residue_l_in,
                    char4 residue_m_in);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, char4 atom_k_in,
                    char4 atom_l_in, char4 atom_m_in, char4 residue_i_in, char4 residue_j_in,
                    char4 residue_k_in, char4 residue_l_in, char4 residue_m_in,
                    const std::vector<double> &surface_in, const std::vector<int2> &locations_in);
  /// \}

  /// \brief Get the force field parameter kind
  ParameterKind getKind() const;
  
  /// \brief Get the atom name of the I atom in the term
  char4 getNameOfAtom(char atom_rank = 'I') const;
  
  /// \brief Get the atom type of atom I in the term
  char4 getTypeOfAtom(char atom_rank = 'I') const;
    
  /// \brief Get the resiude name of atom I in the term
  char4 getNameOfResidue(char atom_rank = 'I') const;

  /// \brief Get the charge of an atom with a given atom name and residue nam
  double getCharge() const;

  /// \brief Get the Lennard-Jones sigma parameter of an atom
  double getSigma() const;

  /// \brief Get the Lennard-Jones epsilon parameter of an atom
  double getEpsilon() const;

  /// \brief Get the Lennard-Jones rho parameter of an atom (the third parameter, for 12-6-4
  ///        potentials)
  double getRho() const;

  /// \brief Get the stiffness constant of a bond, angle, Urey-Bradley, or CHARMM improper term.
  double getStiffnessConstant() const;

  /// \brief Get the equilibrium constant of a bond, angle, or Urey-Bradley term.
  double getEquilibriumConstant() const;

  /// \brief Get the amplitude of a cosine-based dihedral term.
  double getAmplitude() const;

  /// \brief Get the phase angle of a cosine-based dihedral or CHARMM improper dihedral term.
  double getPhaseAngle() const;

  /// \brief Get the periodicity of a cosine-based dihedral term
  double getPeriodicity() const;

  /// \brief Get the electrostatic scaling factor for an attenuated 1:4 non-bonded interaction.
  double getElectrostaticScaling() const;
  
  /// \brief Get the van-der Waals scaling factor for an attenuated 1:4 non-bonded interaction.
  double getVanDerWaalsScaling() const;

  /// \brief Get the torsion kind, indicating whether these four atoms describe a proper or
  ///        improper torsion interaction.
  TorsionKind getTorsionKind() const;
  
  /// \brief Get the surface value edits associated with a CMAP term.
  std::vector<double> getSurfaceValues() const;
  
  /// \brief Get the locations of point edits on a CMAP term's energy surface.
  std::vector<int2> getSurfaceIndices() const;  
  
  /// \brief Get the virtual site frame type
  VirtualSiteKind getVirtualSiteFrameType() const;

  /// \brief Test whether the sigma value of a van-der Waals parameter is being modified.
  bool testSigmaModification() const;
  
  /// \brief Test whether the epsilon value of a van-der Waals parameter is being modified.
  bool testEpsilonModification() const;

  /// \brief Test whether the rho value of a van-der Waals parameter is being modified.
  bool testRhoModification() const;

  /// \brief Test whether the stiffness constant of a valence parameter is being modified.
  bool testStiffnessModification() const;

  /// \brief Test whether the equilibrium constant of a valence parameter is being modified.
  bool testEquilibriumModification() const;

  /// \brief Test whether the amplitude of a dihedral parameter is being modified.
  bool testAmplitudeModification() const;

  /// \brief Test whether the phase angle of a dihedral or improper dihedral is being modified.
  bool testPhaseAngleModification() const;

  /// \brief Test whether the periodicity of a dihedral parameter is being modified.
  bool testPeriodicityModification() const;
  
  /// \brief Set the stiffness property of one of the valence terms.
  void setStiffness(double stiffness_in);

  /// \brief Set the equilibrium property of one of the valence terms.
  void setEquilibrium(double equilbrium_in);

  /// \brief Set the phase angle property of one of the valence terms.
  void setPhaseAngle(double phase_angle_in);

  /// \brief Set the amplitude of a cosine-based dihedral term.
  void setAmplitude(double amplitude_in);

  /// \brief Set the periodicity of a cosine-based dihedral term.
  void setPeriodicity(double periodicity_in);

  /// \brief Set the electrostatic scaling factor for an attenuated 1:4 interaction.
  void setChargeScaling(double scaling_in);

  /// \brief Set the van-der Waals scaling factor for an attenuated 1:4 interaction.
  void setVanDerWaalsScaling(double scaling_in);

  /// \brief Set the charge parameter of an atom.
  void setCharge(double charge_in);

  /// \brief Set the sigma parameter of an atom.
  void setSigma(double sigma_in);
  
  /// \brief Set the epsilon parameter of an atom.
  void setEpsilon(double epsilon_in);
  
  /// \brief Set the rho parameter of an atom.
  void setRho(double rho_in);

  /// \brief Apply these force field parameters to any terms found in a specific topology.
  ///
  /// \param ag      The topology to modify
  /// \param policy  The way to respond if the topology has no such parameters
  void apply(AtomGraph *ag, ExceptionResponse policy = ExceptionResponse::SILENT) const;
  
private:
  ParameterKind kind;    ///< The type of parameter, i.e. BOND or VIRTUAL_SITE_FRAME
  char4 atom_name_i;     ///< Atom or atom type name for atom I of some valence term's arrangement,
                         ///<   or the name of a virtual site atom
  char4 atom_name_j;     ///< Atom or atom type name for atom J of some valence term's arrangement,
                         ///<   or the name of the parent atom in a virtual site frame
  char4 atom_name_k;     ///< Atom or atom type name for atom K of some valence term's arrangement,
                         ///<   or the name of the second frame atom in a virtual site frame
  char4 atom_name_l;     ///< Atom or atom type name for atom L of some valence term's arrangement,
                         ///<   or the name of the third frame atom in a virtual site frame
  char4 atom_name_m;     ///< Atom or atom type name for atom M of some valence term's arrangement,
                         ///<   or the name of the fourth frame atom in a virtual site frame
  char4 residue_name_i;  ///< Residue name for atom I of some valence term arrangement, or a
                         ///<   virtual site atom
  char4 residue_name_j;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's parent atom
  char4 residue_name_k;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's second atom
  char4 residue_name_l;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's third atom
  char4 residue_name_m;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's fourth atom

  // General-purpose real-values numbers for keeping this parameter's details.  The information in
  // any of these parameters is optional in the namelist, but it may be necessary to specify at
  // least one in order to alter anything about a parameter.  Each of these properties begins with
  // its "activation" set to false, and the flag must be set, post-construction, in order to show
  // that one of the properties was specified by the user.
  double property_a;  ///< Charge, Lennard-Jones sigma, valence parameter stiffness or amplitude,
                      ///<   or virtual site frame dimension 1
  double property_b;  ///< Lennard-Jones epsilon, valence parameter equilibrium or phase angle, or
                      ///<   virtual site frame dimension 2
  double property_c;  ///< Lennard-Jones tertiary parameter (12-6-4 or Buckingham potential),
                      ///<   dihedral periodicity, or virtual site frame dimension 3
  bool activate_a;    ///< Indication that property_a is a quantity actively specified by the user
  bool activate_b;    ///< Indication that property_b is a quantity actively specified by the user
  bool activate_c;    ///< Indication that property_c is a quantity actively specified by the user
  
  
  // Miscellaneous properties
  TorsionKind torsion_kind;     ///< Kind of torsion parameter this describes: PROPER or IMPROPER
  std::vector<double> surface;  ///< Values for up to four selected points on a CMAP surface
  std::vector<int2> locations;  ///< Location of up to four points on the CMAP surface
  VirtualSiteKind frame_type;   ///< Frame type of the virtual site, if that is what this object
                                ///<   contains.
};
  
} // namespace modeling
} // namespace stormm

#endif
