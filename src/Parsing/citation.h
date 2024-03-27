// -*-c++-*-
#ifndef STORMM_CITATION_H
#define STORMM_CITATION_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace parse {

/// \brief A means for numbering and tracking citations or other documented sources of information
class Citation {
public:

  /// \brief Constructor can take all attributes, but prioritizes the byline and index
  Citation(const std::string &byline_in, int index_in = 0);
  
  /// \brief Get the citation number
  int getCitationIndex() const;

  /// \brief Get the byline
  std::string getByline() const;

  /// \brief Add one author to the citation
  ///
  /// \param author_in  Name of the author to add (a full name is entered as one string)
  void addAuthor(const std::string &author_in);

  /// \brief Add multiple authors to the citation
  ///
  /// \param authors_in  Names of the authors to add (a full name is entered as one string)
  void addAuthors(const std::vector<std::string> &authors_in);

  /// \brief Add one editor to the citation
  ///
  /// \param editor_in  Name of the editor to add (a full name is entered as one string)
  void addEditor(const std::string &editor_in);

  /// \brief Add multiple editors to the citation
  ///
  /// \param editors_in  Names of the editors to add (a full name is entered as one string)
  void addEditors(const std::vector<std::string> &editors_in);

  /// \brief Add the title of the publication or book
  ///
  /// \param title_in  The title to add
  void addTitle(const std::string &title_in);

  /// \brief Add the name of the journal (if applicable)
  ///
  /// \param journal_in  The journal name to add
  void addJournal(const std::string &journal_in);

  /// \brief Add the name of the publisher (if applicable)
  ///
  /// \param publisher_in  The publisher name to add
  void addPublisher(const std::string &publisher_in);

  /// \brief Add the journal volume
  ///
  /// \param volume_in  The journal volume number
  void addVolume(int volume_in);
  
  /// \brief Add the journal issue
  ///
  /// \param issue_in  The journal issue number
  void addIssue(int issue_in);
  
  /// \brief Add the starting page (for a journal or book)
  ///
  /// \param page_start_in  The starting page number
  void addStartPage(int page_start_in);
  
  /// \brief Add the ending page (for a journal or book)
  ///
  /// \param page_end_in  The ending page number
  void addEndPage(int page_end_in);
  
  /// \brief Add the year of publication
  ///
  /// \param year_in  The year that the piece was published
  void addYear(int year_in);
  
  /// \brief Add the digital object identifier
  ///
  /// \param part_a  First part of the DOI
  /// \param part_b  Second part of the DOI
  void addDoI(const std::string &part_a, const std::string &part_b);

private:
  int index;                               ///< Citation number (for the developer to rapidly track
                                           ///<   this citation)
  std::string byline;                      ///< Brief description of the source material
  std::vector<std::string> authors;        ///< List of authors
  std::vector<std::string> editors;        ///< List of editors
  std::string title;                       ///< The title of the publication
  std::string journal;                     ///< Journal title (if applicable)
  std::string publisher;                   ///< Publishing house (if applicable)
  int volume;                              ///< Journal volume (if applicable)
  int issue;                               ///< Issue within the journal volume
  int page_start;                          ///< Starting page
  int page_end;                            ///< Final page (some journals do not have final page
                                           ///<   numbers)
  int year;                                ///< Year of publication
  std::string digital_object_identifier_a; ///< First part of the Digital Object Identifier
  std::string digital_object_identifier_b; ///< Second part of the Digital Object Identifier
};

} // namespace parse
} // namespace stormm

#endif

