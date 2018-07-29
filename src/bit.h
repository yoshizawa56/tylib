#include <Eigen/Core>
#include <vector>
#include <string>
#include <fstream>

namespace tylib{
	
	#ifndef INCLUDE_COMBINATION
	#define INCLUDE_COMBINATION
	long cal_combination(int L, int N);
	#endif
	
	std::vector<std::pair<int, int> > periodic_1d_bond(int L);
	std::vector<std::pair<int, int> > open_1d_bond(int L);
	std::vector<std::pair<int, int> > periodic_square_bond(int H, int W);
	std::vector<std::pair<int, int> > open_square_bond(int H, int W);
	
	class lattice{
	  public:
		lattice() {}
		lattice(int L, std::vector<std::pair<int, int> > nearest_bond);
		lattice(int L, std::string text_name);
		void init(int L, std::vector<std::pair<int, int> > nearest_bond);
		void init(int L, std::string text_name);
		//n次近接ボンドを返す
		std::vector<std::pair<int, int> > ref_bond(int n);
		
		lattice& operator=(const lattice& tmp){
			return *this;
		}
		
	  private:
		int L;
		std::vector<std::vector<std::pair<int, int> > > neighboring_bond;
		
		void cal_bond(std::vector<std::pair<int, int> > &nearest_bond);
	};
	
	class computational_basis{
	  public:
		computational_basis(int L, int N);
		computational_basis() {}
		computational_basis(int L, int N, std::vector<std::pair<int, int> > nearest_bond);
		void init(int L, int N);
		void init(int L, int N, std::vector<std::pair<int, int> > nearest_bond);
		long index(Eigen::VectorXi &bits);
		Eigen::VectorXi bits(long index);
		std::pair<Eigen::VectorXi, Eigen::VectorXi> fermion_bits(long index);
		int ref_L();
		long ref_d();
		int ref_N();
		bool isConserved();
		lattice ref_lattice();
		std::vector<std::pair<int, int> > ref_bond(int n);
		
	  private:
		Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic> C;
		int ones;
		int zeros;
		int L;
		int N;
		long d;
		bool conserved;
		lattice Lattice;
		
	};
}