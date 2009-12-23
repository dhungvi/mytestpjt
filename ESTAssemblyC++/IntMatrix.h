#ifndef ZY_IntMatrix_HH
#define ZY_IntMatrix_HH

#include <vector>

class IntMatrix {
private:
	std::vector<std::vector<int> > matrix;
public:
	IntMatrix() {matrix = std::vector<std::vector<int> >(0,0);}

	// construct with dimensions rows x cols
	IntMatrix(int n, int m) {
		matrix = std::vector<std::vector<int> > (n, std::vector<int>(m));
	}

	// return number of rows
	inline int rows() const {
		return (this->matrix).size();
	}

	// return number of cols
	inline int cols() const {
		return (this->matrix)[0].size();
	}

	inline int get(int row, int col) const {
		return (this->matrix)[row][col];
	}

	inline void set(int row, int col, int value) {
		(this->matrix)[row][col] = value;
	}


};


#endif
