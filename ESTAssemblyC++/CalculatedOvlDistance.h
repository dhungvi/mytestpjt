#ifndef ZY_CalculatedOvlDistance_HH
#define ZY_CalculatedOvlDistance_HH

#include <string>
#include <vector>
#include <map>

struct Key {
	int idx1;
	int idx2;
};
struct Dist {
	int ovlDis;
	int ovlLen;
};

class CalculatedOvlDistance {
public:
	void addDistance(int i1, int i2, int ovlDis, int ovlLen);
	std::vector<int> searchDistance(int i1, int i2);

private:
	std::map<Key, Dist> distances;
};

bool operator<(const Key& k1, const Key& k2);

#endif

