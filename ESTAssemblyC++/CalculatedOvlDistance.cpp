#include "CalculatedOvlDistance.h"
#include <limits.h>
using namespace std;

void CalculatedOvlDistance::addDistance(int i1, int i2, int ovlDis, int ovlLen) {
	Key key;
	Dist dist;
	key.idx1 = i1;
	key.idx2 = i2;
	dist.ovlDis = ovlDis;
	dist.ovlLen = ovlLen;
	distances.insert(pair<Key, Dist>(key, dist));
}

/*
 * search in distances to find if there is information between i1 and i2.
 * If does, return int[2] which is same as the return from "getOVLDistance" in OvlDistance.java.
 * int[0] - overlap length; int[1] - overlap distance.
 * If not, return [INT_MAX, INT_MAX].
 */
vector<int> CalculatedOvlDistance::searchDistance(int i1, int i2) {
	vector<int> ret(2);
	Key key1, key2;
	key1.idx1 = i1;
	key1.idx2 = i2;
	key2.idx1 = i2;
	key2.idx2 = i1;
	map<Key,Dist>::iterator it1 = distances.find(key1);
	map<Key,Dist>::iterator it2 = distances.find(key2);

	if (it1 != distances.end()) {
		ret[0] = (it1->second).ovlLen;
		ret[1] = (it1->second).ovlDis;
	} else if (it2 != distances.end()){
		ret[0] = -1 * ((it2->second).ovlLen);
		ret[1] = -1 * ((it2->second).ovlDis);
	} else {
		ret[0] = INT_MAX;
		ret[1] = INT_MAX;
	}
	return ret;

}


bool operator<(const Key& k1, const Key& k2) {
	if (k1.idx1 < k2.idx1)
		return true;
	else if ((k1.idx1 == k2.idx1) && (k1.idx2 < k2.idx2))
		return true;
	else
		return false;
}
