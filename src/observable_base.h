#include <vector>
#include <Eigen/Sparse>
#include <complex>

#ifndef INCLUDE_BIT
#define INCLUDE_BIT
#include "bit.h"
#endif

namespace tylib{
	
	std::vector<std::pair<int, double> > cal_hopping(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double J);
	std::vector<std::pair<int, double> > cal_hopping(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, std::vector<double> &J);
	std::vector<std::pair<int, double> > cal_HCB_many_body_observable(int l, computational_basis &basis, std::vector<std::vector<std::pair<int, int> > > &obs, std::vector<double> &J);
	std::vector<std::pair<int, double> > cal_Sx(int l, computational_basis &basis, double h);		//一様な横磁場
	std::vector<std::pair<int, double> > cal_Sx_interaction(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double U);
	std::vector<std::pair<int, double> > cal_spin_many_body_observable(int l, computational_basis &basis, std::vector<std::vector<std::pair<int, int> > > &obs, std::vector<double> &J);
	std::vector<std::pair<int, double> > cal_fermion_hopping(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double J);
	std::vector<std::pair<int, double> > cal_fermion_many_body_observable(int l, computational_basis &basis, std::vector<std::vector<std::pair<int, int> > > &obs, std::vector<double> &J);
	
	double cal_Sz(int l, computational_basis &basis, double h);
	double cal_Sz_disorder(int l, computational_basis &basis, std::vector<double> &h);
	double cal_Sz_interaction(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double U);
	double cal_interaction(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double U);
	double cal_interaction_disorder(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, std::vector<double> &U);
	double cal_onsite(int l, computational_basis &basis, std::vector<double> &potential);
	std::vector<std::pair<int, int> > make_bond(int L, int neighbor);
	
	void SparseMatrix_to_Triplets(std::vector<Eigen::Triplet<double> > &elements, Eigen::SparseMatrix<double> &M);
	void SparseMatrix_to_Triplets(std::vector<Eigen::Triplet<std::complex<double> > > &elements, Eigen::SparseMatrix<std::complex<double> > &M);
	
	
	//物理量はobservable_baseを継承して定義すること．コンストラクタおよびcal_indexesの定義のみを書けばよい
	class observable_base{
	  public:
		observable_base(){}
		virtual std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis) {return std::vector<std::pair<int, double> >();}
		void cal_elements_append(int index, computational_basis &basis, std::vector<std::pair<int, double> > &elements);
		void cal_matrix(Eigen::SparseMatrix<double> &H, computational_basis &basis);
		virtual std::pair<bool, double> coeff(int i, int j, computational_basis &basis) {return std::make_pair(false,0.0);}
		
		observable_base& operator=(const observable_base& tmp){
			return *this;
		}
		
	};
	
	class observable_sum : public observable_base{
	  public:
		observable_sum() {}
		observable_sum(observable_base &O1, observable_base &O2);
		observable_sum(observable_base &O1, double coeff1, observable_base &O2, double coeff2);
		observable_sum(std::vector<observable_base*> &O);
		void append_observable(observable_base &O);
		void append_observable(observable_base &O, double coef);
		void reset();
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		std::pair<bool, double> coeff(int i, int j, computational_basis &basis);
		
	  private:
		std::vector<observable_base*> O;
		std::vector<double> coef;
	};
	
	// 物理量の積 O1O2 を定義
	class observable_product : public observable_base{
	  public:
		observable_product() {}
		observable_product(observable_base &O1, observable_base &O2);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		
	  private:
		observable_base O1;
		observable_base O2;
		
	};
	
	class observable_base_complex{
	  public:
		observable_base_complex(){}
		virtual std::vector<std::pair<int, std::complex<double> > > cal_elements(int index, computational_basis &basis) {return std::vector<std::pair<int, std::complex<double> > >();}
		void cal_elements_append(int index, computational_basis &basis, std::vector<std::pair<int, std::complex<double> > > &elements);
		void cal_matrix(Eigen::SparseMatrix<std::complex<double> > &H, computational_basis &basis);
		//virtual std::pair<bool, std::complex<double> > coeff(int i, int j, computational_basis &basis) {return std::make_pair<bool, std::complex<double> >(false,0);}
		
		observable_base_complex& operator=(const observable_base_complex& tmp){
			return *this;
		}
		
	};
	
	
	class observable_sum_complex : public observable_base_complex{
	  public:
		observable_sum_complex(){}
		observable_sum_complex(observable_base_complex &O1, observable_base_complex &O2);
		observable_sum_complex(observable_base &O1, observable_base_complex &O2);
		observable_sum_complex(observable_base_complex &O1, observable_base &O2);
		observable_sum_complex(observable_base &O1, observable_base &O2);
		observable_sum_complex(std::vector<observable_base_complex*> &O);
		void append_observable(observable_base &O);
		void append_observable(observable_base_complex &O);
		void reset();
		std::vector<std::pair<int, std::complex<double> > > cal_elements(int index, computational_basis &basis);
		
	  private:
		std::vector<observable_base_complex*> Ocd;
		std::vector<observable_base*> Od;
	};
	
	
}