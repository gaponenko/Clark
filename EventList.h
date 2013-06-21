// Andrei Gaponenko, 2013

#ifndef EventList_h
#define EventList_h

#include <vector>
#include <string>

class EventClass;

class EventList {
public:

  EventList(const std::string& infile = "");

  bool requested(int run, int event) const;
  bool requested(const EventClass& evt) const;

  bool empty() const { return list_.empty(); }

private:

  struct EventID {
    int r;
    int e;
    EventID(int run, int event) : r(run), e(event) {}
    bool operator<(const EventID& b) const {
      return (r < b.r) ||((r==b.r) && (e < b.e));
    }
  };

  typedef std::vector<EventID> List;
  List list_;
};

extern EventList gEventList;

#endif/*EventList_h*/
