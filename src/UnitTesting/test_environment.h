// -*-c++-*-
#ifndef STORMM_TEST_ENVIRONMENT_H
#define STORMM_TEST_ENVIRONMENT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Namelists/command_line_parser.h"
#include "unit_test_enumerators.h"

namespace stormm {
namespace testing {

using constants::ExceptionResponse;
using namelist::CommandLineParser;
  
/// \brief Object for parsing as well as storing environment and command-line input relevant to
///        unit testing (Test Driven Development) in STORMM.  Environment variables will be sought,
///        but command-line inputs addressing the same information will take precedence.
class TestEnvironment {
public:

  /// \brief The TestEnvironment constructor takes command line arguments but also searches for
  ///        native environment variables.
  ///
  /// \param argc             Number of command line arguments, passed down from the calling
  ///                         program
  /// \param argv             List of command-line arguments
  /// \param tmpdir_required  Confirmation that the temporary directory is required for tests.  If
  ///                         not required and the temporary directory cannot be created or written
  ///                         to, no warnings will be printed.  The directory's status will still
  ///                         be noted.
  /// \param policy           Policy to take in the event that a command line argument is not
  ///                         recognized
  /// \{
  TestEnvironment(int argc, const char* argv[], CommandLineParser *clip,
                  TmpdirStatus tmpdir_required = TmpdirStatus::NOT_REQUIRED,
                  ExceptionResponse policy = ExceptionResponse::WARN);

  TestEnvironment(int argc, const char* argv[],
                  TmpdirStatus tmpdir_required = TmpdirStatus::NOT_REQUIRED,
                  ExceptionResponse policy = ExceptionResponse::WARN);

  TestEnvironment(int argc, const char* argv[], ExceptionResponse policy);
  /// \}

  /// \brief Destructor for the TestEnvironment object must delete the temporary directory if it
  ///        was created by the constructor.  This can only be done after STORMM removes files in
  ///        the temporary directory known to have been created by the program itself.  The actual
  ///        directory removal is accomplished by a syscall to /bin/rmdir, which will not delete a
  ///        directory with contents still in it.  In this way, the temporary directory can be
  ///        removed safely without destroying valuable work left by some other user or program.
  ~TestEnvironment();

  /// \brief Get the verbosity level requested for this test environment.
  TestVerbosity getVerbosity() const;

  /// \brief Get the general tolerance that can be applied in blanket fashion across many tests.
  ///        This is most useful for regression tests or much higher-level unit tests, as low-level
  ///        unit test results are often known to very high precision.
  double getTolerance() const;

  /// \brief Get the random seed from the test environment.  Separate random seeds can be supplied
  ///        to any random number generator object at the developer's discretion.
  int getRandomSeed() const;

  /// \brief Get the STORMM home (test executable and installed library) path.
  std::string getStormmHomePath() const;

  /// \brief Get the STORMM source root path.
  std::string getStormmSourcePath() const;

  /// \brief Get the temporary directory path that the test environment can use to store files.
  std::string getTemporaryDirectoryPath() const;

  /// \brief Check whether the temporary directory is accessible to the program.
  bool getTemporaryDirectoryAccess() const;

  /// \brief Get the list of files that this test program has been given responsibility.  It will
  ///        assume that these files are things it created temporarily and that they should be
  ///        removed when the TestEnvironment object falls out of scope.
  std::vector<std::string> getListOfFilesCreated() const;

  /// \brief Return whether to record a snapshot of the current data into the named file, or
  ///        compare the data to an existing file.
  SnapshotOperation takeSnapshot() const;

  /// \brief Return the directive on displaying timings at the end of a test program.
  bool getDisplayTimingsOrder() const;

  /// \brief Set the file removal flag.
  ///
  /// \param remove_files_in  The new flag value
  void setFileRemoval(const bool remove_files_in);
  
  /// \brief Catalog that a file has been created by the program.  This function will not actually
  ///        check that the file was created by code recently executed; it just serves as a way to
  ///        accumulate a list of files that should be regarded as the property of whatever program
  ///        using the STORMM TestEnvironment object.  The program will, however, detect whether an
  ///        absolute path is in use and make the path absolute if not.
  ///
  /// \param path  Name of the file of interest
  void logFileCreated(const std::string &path);

  /// \brief Catalog that a directory has been created by the program.  As above, this function
  ///        will take the developer at their word and accumulates a list of directories that
  ///        should be regarded as the property of the TestEnvironment object and removed, if
  ///        possible, when the object is destroyed.
  ///
  /// Overloaded:
  ///   - Take a single path name
  ///   - Take multiple path names (this functionality exists in contrast to logFileCreated, as
  ///     making a new directory can involve creating a series of new directories)
  ///
  /// \param path   Name of the directory of interest
  /// \param paths  List of new directories
  /// \{
  void logDirectoryCreated(const std::string &path);
  void logDirectoryCreated(const std::vector<std::string> &paths);
  /// \}

private:

  /// A command line parsing object with an internal namelist
  CommandLineParser user_mods;
  
  /// Reporting level for test summaries
  TestVerbosity verbose_level;

  /// Baseline tolerance for all tests (provided for easy manipulation of test precision, but only
  /// used if the test calls upon this variable, for example as a tolerance to an Approx object)
  double general_tolerance;

  /// Random number generator seed to use (this gives users command-line control over this seed,
  /// if that is desirable)
  int random_seed;

  /// STORMM installation home path
  std::string stormm_home_path;

  /// STORMM source code path
  std::string stormm_source_path;

  /// Path to temporary work directory (this will be created if not currently available when the
  /// object is constructed, then destroyed by the object's destructor)
  std::string tmpdir_path;

  /// Indicator of whether the temporary directory managed by this object was created during the
  /// object construction.  If so, the temporary directory will then be destroyed upon the
  /// object's destruction.
  bool tmpdir_created;

  /// Indicator of whether the temporary directory is actually writeable (tests that depend on the
  /// directory should be skipped and marked as such)
  bool tmpdir_writeable;

  /// List of files created durign testing and catalogged by this object.  Any files placed in
  /// the temporary directory should be catalogged here to ensure its timely destruction, as only
  /// files listed here can be removed by the program prior to a /bin/rmdir call.
  std::vector<std::string> files_created;

  /// List of directories created on the way to making (and including) the temporary directory
  std::vector<std::string> directories_created;

  /// Indicator of whether to remove files created during testing and logged by this object
  bool remove_files;

  /// Indicator of whether to snapshot current results into the saved test files.  Tests can be
  /// run with snapshot deposition turned on, available only via a command-line option to the
  /// specific test program, and this will cause the requested files to be written based on the
  /// available data rather than read for immediate comparison.
  SnapshotOperation snapshot_behavior;

  /// Directive on whether to display timings at the end of a test program
  bool display_time;
};

} // namespace testing
} // namespace stormm

#endif
