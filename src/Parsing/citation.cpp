#include "copyright.h"
#include "citation.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
Citation::Citation(const std::string &byline_in, const int index_in) :
    index{index_in}, byline{byline_in}, authors{}, editors{}, title{}, journal{}, publisher{},
    volume{-1}, issue{-1}, page_start{-1}, page_end{-1}, year{-1},
    digital_object_identifier_a{}, digital_object_identifier_b{}
{}

//-------------------------------------------------------------------------------------------------
int Citation::getCitationIndex() const {
  return index;
}

//-------------------------------------------------------------------------------------------------
std::string Citation::getByline() const {
  return byline;
}

//-------------------------------------------------------------------------------------------------
void Citation::addAuthor(const std::string &author_in) {
  authors.push_back(author_in);
}

//-------------------------------------------------------------------------------------------------
void Citation::addAuthors(const std::vector<std::string> &authors_in) {
  authors.insert(authors.end(), authors_in.begin(), authors_in.end());
}

//-------------------------------------------------------------------------------------------------
void Citation::addEditor(const std::string &editor_in) {
  editors.push_back(editor_in);
}

//-------------------------------------------------------------------------------------------------
void Citation::addEditors(const std::vector<std::string> &editors_in) {
  editors.insert(editors.end(), editors_in.begin(), editors_in.end());
}

//-------------------------------------------------------------------------------------------------
void Citation::addTitle(const std::string &title_in) {
  title = title_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addJournal(const std::string &journal_in) {
  journal = journal_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addPublisher(const std::string &publisher_in) {
  publisher = publisher_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addVolume(const int volume_in) {
  volume = volume_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addIssue(const int issue_in) {
  issue = issue_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addStartPage(const int page_start_in) {
  page_start = page_start_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addEndPage(const int page_end_in) {
  page_end = page_end_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addYear(const int year_in) {
  year = year_in;
}

//-------------------------------------------------------------------------------------------------
void Citation::addDoI(const std::string &part_a, const std::string &part_b) {
  digital_object_identifier_a = part_a;
  digital_object_identifier_b = part_b;
}

} // namespace parse
} // namespace stormm
