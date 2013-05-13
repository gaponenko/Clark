#include "EventList.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "EventClass.h"


EventList gEventList;

//================================================================
EventList::EventList(const std::string& infile) {
  if(!infile.empty()) {
    std::ifstream file(infile.c_str());
    if(!file) {
      throw std::runtime_error("EventList(): ERROR opening input file "+infile);
    }

    // Expected line format:
    //    run event [anything else]
    std::string line;
    while(std::getline(file, line)) {
      std::istringstream is(line);
      int run(0), event(0);
      if(is>>run>>event) {
        list_.push_back(EventID(run, event));
      }
      else {
        throw std::runtime_error("EventList: ERROR: bad input line: "+line);
      }
    }
  }
}

//================================================================
bool EventList::requested(int run, int event) const {
  return std::binary_search(list_.begin(), list_.end(), EventID(run, event));
}

//================================================================
bool EventList::requested(const EventClass& evt) const {
  return requested(evt.nrun, evt.nevt);
}

//================================================================
