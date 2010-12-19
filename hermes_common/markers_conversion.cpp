
#include "markers_conversion.h"
#include "hermes_logging.h"

using namespace std;


MarkersConversion::MarkersConversion()
{
  conversion_table_for_element_markers = new std::map<int, std::string>;
  conversion_table_for_boundary_markers = new std::map<int, std::string>;
  conversion_table_for_element_markers_inverse = new std::map<std::string, int>;
  conversion_table_for_boundary_markers_inverse = new std::map<std::string, int>;

  min_boundary_marker_unused = 1;
  min_element_marker_unused = 0;
}

MarkersConversion::~MarkersConversion()
{
  delete conversion_table_for_element_markers;
  delete conversion_table_for_boundary_markers;
  delete conversion_table_for_element_markers_inverse;
  delete conversion_table_for_boundary_markers_inverse;
}

void MarkersConversion::insert_element_marker(int internal_marker, std::string user_marker)
{
    // First a check that the string value is not already present.
  if(user_marker != "")
    if(conversion_table_for_element_markers_inverse->find(user_marker) != conversion_table_for_element_markers_inverse->end())
      return;
  if(conversion_table_for_element_markers->size() == 0 || conversion_table_for_element_markers->find(internal_marker) == conversion_table_for_element_markers->end()) {
    conversion_table_for_element_markers->insert(std::pair<int, std::string>(internal_marker, user_marker));
    conversion_table_for_element_markers_inverse->insert(std::pair<std::string, int>(user_marker, internal_marker));
    if(user_marker != "")
      this->min_element_marker_unused++;
  }
  return;
}

void MarkersConversion::insert_boundary_marker(int internal_marker, std::string user_marker)
{
  // First a check that the string value is not already present.
  if(user_marker != "")
    if(conversion_table_for_boundary_markers_inverse->find(user_marker) != conversion_table_for_boundary_markers_inverse->end())
      return;
  if(conversion_table_for_boundary_markers->size() == 0 || conversion_table_for_boundary_markers->find(internal_marker) == conversion_table_for_boundary_markers->end()) {
    conversion_table_for_boundary_markers->insert(std::pair<int, std::string>(internal_marker, user_marker));
    conversion_table_for_boundary_markers_inverse->insert(std::pair<std::string, int>(user_marker, internal_marker));
    if(user_marker != "")
      this->min_boundary_marker_unused++;
  }
  return;
}

std::string MarkersConversion::get_user_element_marker(int internal_marker) 
{
  if(conversion_table_for_element_markers->find(internal_marker) == conversion_table_for_element_markers->end())
    error("MarkersConversions class asked for a non existing marker %d", internal_marker);
  return conversion_table_for_element_markers->find(internal_marker)->second;
}
std::string MarkersConversion::get_user_boundary_marker(int internal_marker) 
{
  if(conversion_table_for_boundary_markers->find(internal_marker) == conversion_table_for_boundary_markers->end())
    error("MarkersConversions class asked for a non existing marker %d", internal_marker);
  return conversion_table_for_boundary_markers->find(internal_marker)->second;
}

int MarkersConversion::get_internal_element_marker(std::string user_marker) 
{
  if(conversion_table_for_element_markers_inverse->find(user_marker) == conversion_table_for_element_markers_inverse->end())
    error("MarkersConversions class asked for a non existing marker %s", user_marker);
  return conversion_table_for_element_markers_inverse->find(user_marker)->second;
}
int MarkersConversion::get_internal_boundary_marker(std::string user_marker) 
{
  if(conversion_table_for_boundary_markers_inverse->find(user_marker) == conversion_table_for_boundary_markers_inverse->end())
    error("MarkersConversions class asked for a non existing marker %s", user_marker);
  return conversion_table_for_boundary_markers_inverse->find(user_marker)->second;
}

void MarkersConversion::check_boundary_marker(int marker)
{
  // This marker is already present.
  if(conversion_table_for_boundary_markers->find(marker) != conversion_table_for_boundary_markers->end())
    // If the present boundary marker is a different one (user-supplied string).
    if(conversion_table_for_boundary_markers->find(marker)->second != "")
    {
      // We have to reassign the markers.
      std::string temp_user_marker = conversion_table_for_boundary_markers->find(marker)->second;
      conversion_table_for_boundary_markers->erase(marker);
      conversion_table_for_boundary_markers_inverse->erase(temp_user_marker);
      insert_boundary_marker(this->min_boundary_marker_unused, temp_user_marker);
    }
  // Now we need to check if the variable min_boundary_marker_unused does not collide with
  // this marker, because we will need to raise the variable if it does as
  // we are about to insert this user-supplied integer marker.
  if(marker == min_boundary_marker_unused)
    min_boundary_marker_unused++;
}

void MarkersConversion::check_element_marker(int marker)
{
  // This marker is already present.
  if(conversion_table_for_element_markers->find(marker) != conversion_table_for_element_markers->end())
    // If the present element marker is a different one (user-supplied string).
    if(conversion_table_for_element_markers->find(marker)->second != "")
    {
      // We have to reassign the markers.
      std::string temp_user_marker = conversion_table_for_element_markers->find(marker)->second;
      conversion_table_for_element_markers->erase(marker);
      conversion_table_for_element_markers_inverse->erase(temp_user_marker);
      insert_element_marker(this->min_element_marker_unused, temp_user_marker);
    }
  // Now we need to check if the variable min_element_marker_unused does not collide with
  // this marker, because we will need to raise the variable if it does as
  // we are about to insert this user-supplied integer marker.
  if(marker == min_element_marker_unused)
    min_element_marker_unused++;
}
