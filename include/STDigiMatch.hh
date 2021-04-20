#ifndef S1DigiMatch_h
#define S1DigiMatch_h 1

#include <vector>

class STDigiMatch
{
 public:
  STDigiMatch();
  virtual ~STDigiMatch() {}
  void AddMcPointIndex(int index) { fMcPointIndices.push_back(index); }
  int GetNIndices() { return fMcPointIndices.size(); }
  int GetMcPointIndex(int i) { return fMcPointIndices[i]; }

 private:
  std::vector<int> fMcPointIndices;
};

#endif
