#include "observable_base.h"
#include <iostream>

namespace tylib{
	
	std::vector<std::pair<int, double> > cal_hopping(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double J){
		
		if(J == 0) return std::vector<std::pair<int, double> >(0);
		
		int L = basis.ref_L();
		
		std::vector<std::pair<int, double> > elements;
		
		Eigen::VectorXi bits = basis.bits(l);
		int nbonds=bonds.size();
		for(int i = 0; i < nbonds; i++){
			Eigen::VectorXi tmp = bits;
			
			if( bits(bonds[i].first) == 1 && bits(bonds[i].second) == 0 ){
				tmp(bonds[i].first) = 0;
				tmp(bonds[i].second) = 1;
				
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, J));
			}
			
			if( bits(bonds[i].first) == 0 && bits(bonds[i].second) == 1 ){
				tmp(bonds[i].first) = 1;
				tmp(bonds[i].second) = 0;
				
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, J));
			}
		}
		
		return elements;
		
		
	}
	
	std::vector<std::pair<int, double> > cal_hopping_disorderd(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, std::vector<double> &J){
		
		int L = basis.ref_L();
		
		std::vector<std::pair<int, double> > elements;
		
		Eigen::VectorXi bits = basis.bits(l);
		
		int nbonds=bonds.size();
		for(int i = 0; i < nbonds; i++){
			Eigen::VectorXi tmp = bits;
			
			if( bits(bonds[i].first) == 1 && bits(bonds[i].second) == 0 ){
				tmp(bonds[i].first) = 0;
				tmp(bonds[i].second) = 1;
				
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, J[i]));
			}
			
			if( bits(bonds[i].first) == 0 && bits(bonds[i].second) == 1 ){
				tmp(bonds[i].first) = 1;
				tmp(bonds[i].second) = 0;
				
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, J[i]));
			}
		}
		
		return elements;
		
		
	}
	
	std::vector<std::pair<int, double> > cal_HCB_many_body_observable(int l, computational_basis &basis, std::vector<std::vector<std::pair<int, int> > > &obs, std::vector<double> &J){
		std::vector<std::pair<int, double> > elements;
		
		for(int i = 0; i < obs.size(); i++){
			Eigen::VectorXi bits = basis.bits(l);
			
			Eigen::VectorXi tmp = bits;
			bool a = true;
			
			for(int j = obs[i].size() - 1; j >= 0 ; j--){
				if(obs[i][j].second == 0){		//anihilation
					if(tmp(obs[i][j].first) == 0){
						a = false;
						break;
					}
					tmp(obs[i][j].first) = 0;
				}else{							//creation
					if(tmp(obs[i][j].first) == 1){
						a = false;
						break;
					}
					tmp(obs[i][j].first) = 1;
				}
			}
			if(a){
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, J[i]));
			}
		}
		
		return elements;
	}
	
	std::vector<std::pair<int, double> > cal_spin_many_body_observable(int l, computational_basis &basis, std::vector<std::vector<std::pair<int, int> > > &obs, std::vector<double> &J){
		std::vector<std::pair<int, double> > elements;
		
		for(int i = 0; i < obs.size(); i++){
			Eigen::VectorXi bits = basis.bits(l);
			
			Eigen::VectorXi tmp = bits;
			bool a = true;
			double val = J[i];
			
			for(int j = obs[i].size() - 1; j >= 0 ; j--){
				if(obs[i][j].second == 0){			//S_-
					if(tmp(obs[i][j].first) == 0){
						a = false;
						break;
					}
					tmp(obs[i][j].first) = 0;
				}else if(obs[i][j].second == 1){	//S_+
					if(tmp(obs[i][j].first) == 1){
						a = false;
						break;
					}
					tmp(obs[i][j].first) = 1;
				}else if(obs[i][j].second == 2){	//S_z
					if(tmp(obs[i][j].first) == 0){
						val *= -1;
					}
				}else if(obs[i][j].second == 3){	//S_x
					if(basis.isConserved()){
						a = false;
						break;
					}
					tmp(obs[i][j].first) = tmp(obs[i][j].first) ^ 1;
				}
				
			}
			if(a){
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, val));
			}
		}
		
		return elements;
		
	}
	
	std::vector<std::pair<int, double> > cal_Sx(int l, computational_basis &basis, double h){
		if(basis.isConserved()) return std::vector<std::pair<int, double> >();
		std::vector<std::pair<int, double> > elements;
		
		int L = basis.ref_L();
		Eigen::VectorXi bits = basis.bits(l);
		for(int i = 0; i < L; i++){
			Eigen::VectorXi tmp = bits;
			if(bits(i) == 0){
				tmp(i) = 1;
			}else{
				tmp(i) = 0;
			}
			
			int index = basis.index(tmp);
			elements.push_back(std::make_pair(index, h));
		}
		
		return elements;
	}
	
	std::vector<std::pair<int, double> > cal_Sx_interaction(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double U){
		
		if(basis.isConserved()) return std::vector<std::pair<int, double> >();
		
		std::vector<std::pair<int, double> > elements;
		
		int L = basis.ref_L();
		Eigen::VectorXi bits = basis.bits(l);
		for(int i = 0; i < bonds.size(); i++){
			Eigen::VectorXi tmp = bits;
			if(bits(bonds[i].first) == 0){
				tmp(bonds[i].first) = 1;
			}else{
				tmp(bonds[i].first) = 0;
			}
			
			if(bits(bonds[i].second) == 0){
				tmp(bonds[i].second) = 1;
			}else{
				tmp(bonds[i].second) = 0;
			}
			
			int index = basis.index(tmp);
			//std::cout << index << std::endl;
			elements.push_back(std::make_pair(index, U));
		}
		
		return elements;
	}
	
	std::vector<std::pair<int, double> > cal_fermion_hopping(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double J){
		
		if(J == 0) return std::vector<std::pair<int, double> >(0);
		
		int L = basis.ref_L();
		
		std::vector<std::pair<int, double> > elements;
		
		std::pair<Eigen::VectorXi, Eigen::VectorXi> bits = basis.fermion_bits(l);
		
		int nbonds=bonds.size();
		for(int i = 0; i < nbonds; i++){
			Eigen::VectorXi tmp = bits.first;
			int site1, site2;
			if(bonds[i].first < bonds[i].second){
				site1 = bonds[i].first;
				site2 = bonds[i].second;
			}
			
			if( bits.first(site1) == 1 && bits.first(site2) == 0 ){
				tmp(site1) = 0;
				tmp(site2) = 1;
				
				int index = basis.index(tmp);
				
				int sign = 1 ^ bits.second(site1) ^ bits.second(site2);
				
				double val;
				if(sign == 0){
					val = J;
				}else{
					 val = -J;
				}
				
				elements.push_back(std::make_pair(index, val));
			}
			
			if( bits.first(site1) == 0 && bits.first(site2) == 1 ){
				tmp(site1) = 1;
				tmp(site2) = 0;
				
				int index = basis.index(tmp);
				
				int sign = bits.second(site1) ^ bits.second(site2);
				
				double val;
				if(sign == 0){
					val = J;
				}else{
					 val = -J;
				}
				
				elements.push_back(std::make_pair(index, val));
			}
		}
		
		return elements;
		
		
	}
	
	std::vector<std::pair<int, double> > cal_fermion_many_body_observable(int l, computational_basis &basis, std::vector<std::vector<std::pair<int, int> > > &obs, std::vector<double> &J){
		std::vector<std::pair<int, double> > elements;
		
		int L = basis.ref_L();
		
		for(int i = 0; i < obs.size(); i++){
			std::pair<Eigen::VectorXi, Eigen::VectorXi> bits_count = basis.fermion_bits(l);
			Eigen::VectorXi bits = bits_count.first;
			Eigen::VectorXi count = bits_count.second;
			
			Eigen::VectorXi tmp = bits;
			Eigen::VectorXi tmp_count = count;
			bool a = true;
			
			double val = J[i];
			
			for(int j = obs[i].size() - 1; j >= 0 ; j--){
				int site = obs[i][j].first;
				if(obs[i][j].second == 0){		//anihilation
					if(tmp(site) == 0){
						a = false;
						break;
					}
					tmp(site) = 0;
					if(tmp_count[site] == 1){
						val *= -1;
					}
					
					for(int j = i; j < L; j++){
						tmp_count[site] = tmp_count[site] ^ 1;
					}
				}else{							//creation
					if(tmp(site) == 1){
						a = false;
						break;
					}
					tmp(site) = 1;
					
					if(tmp_count[site] == 1){
						val *= -1;
					}
					
					for(int j = i; j < L; j++){
						tmp_count[site] = tmp_count[site] ^ 1;
					}
				}
			}
			if(a){
				int index = basis.index(tmp);
				elements.push_back(std::make_pair(index, val));
			}
		}
		
		return elements;
	}
	
	double cal_Sz(int l, computational_basis &basis, double h){
		int L = basis.ref_L();
		Eigen::VectorXi bits = basis.bits(l);
		double val = 0;
		for(int i = 0; i < L; i++){
			if(bits(i) == 1){
				val += h;
			}else{
				val -= h;
			}
		}
		return val;
	}
	
	double cal_Sz_disorder(int l, computational_basis &basis, std::vector<double> &h){
		int L = basis.ref_L();
		Eigen::VectorXi bits = basis.bits(l);
		double val = 0;
		for(int i = 0; i < L; i++){
			if(bits(i) == 1){
				val += h[i];
			}else{
				val -= h[i];
			}
		}
		return val;
		
	}
	
	double cal_Sz_interaction(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double U){
		if(U == 0) return 0;
		
		int L = basis.ref_L();
		
		Eigen::VectorXi bits = basis.bits(l);
		
		int nbonds = bonds.size();
		
		double val = 0;
		
		for(int i = 0; i < nbonds; i++){
			if(bits(bonds[i].first) == bits(bonds[i].second)){
				val += U;
			}else{
				val -= U;
			}
		}
		
		return val;
	}
	
	double cal_interaction(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, double U){
		if(U == 0) return 0;
		
		int L = basis.ref_L();
		
		Eigen::VectorXi bits = basis.bits(l);
		
		int nbonds = bonds.size();
		
		double val = 0;
		
		for(int i = 0; i < nbonds; i++){
			if(bits(bonds[i].first) * bits(bonds[i].second) == 1){
				val += U;
			}
		}
		
		return val;
		
	}
	
	double cal_interaction_disorder(int l, computational_basis &basis, std::vector<std::pair<int, int> > &bonds, std::vector<double> &U){
		int L = basis.ref_L();
		
		Eigen::VectorXi bits = basis.bits(l);
		
		int nbonds = bonds.size();
		
		double val = 0;
		
		for(int i = 0; i < nbonds; i++){
			if(bits(bonds[i].first) * bits(bonds[i].second) == 1){
				val += U[i];
			}
		}
		
		return val;
		
	}
	
	
	double cal_onsite(int l, computational_basis &basis, std::vector<double> &potential){
		int L = basis.ref_L();
		
		Eigen::VectorXi bits = basis.bits(l);
		
		double val = 0;
		
		for(int i = 0; i < L; i++){
			if(bits(i) == 1){
				val += potential[i];
			}
		}
		
		
		return val;
	}
	
	std::vector<std::pair<int, int> > make_bond(int L, int neighbor){
		std::vector<std::pair<int, int> > bonds;
		
		for(int i = 0; i < L-neighbor; i++){
			bonds.push_back(std::make_pair(i, i+neighbor));
		}
		for(int i = 0; i < neighbor; i++){
			bonds.push_back(std::make_pair(L-neighbor+i, i));
		}
		
		return bonds;
		
	}
	
	void SparseMatrix_to_Triplets(std::vector<Eigen::Triplet<double> > &elements, Eigen::SparseMatrix<double> &M){
		
		int *row = M.outerIndexPtr();
		int *outerStart = M.outerIndexPtr();
		double *value = M.valuePtr();
		int d = M.cols();
		
		for(int i = 0; i < d; i++){
			int nonzero = *(outerStart+1) - *outerStart;
			for(int j = 0; j < nonzero; j++){
				elements.push_back(Eigen::Triplet<double>(*row, i, *value));
				row++;
				value++;
			}
			outerStart++;
		}
	}
	
	void SparseMatrix_to_Triplets(std::vector<Eigen::Triplet<std::complex<double> > > &elements, Eigen::SparseMatrix<std::complex<double> > &M){
		
		int *row = M.outerIndexPtr();
		int *outerStart = M.outerIndexPtr();
		std::complex<double> *value = M.valuePtr();
		int d = M.cols();
		
		for(int i = 0; i < d; i++){
			int nonzero = *(outerStart+1) - *outerStart;
			for(int j = 0; j < nonzero; j++){
				elements.push_back(Eigen::Triplet<std::complex<double> >(*row, i, *value));
				row++;
				value++;
			}
			outerStart++;
		}
	}
	
	void observable_base::cal_elements_append(int index, computational_basis &basis, std::vector<std::pair<int, double> > &elements){
		
		std::vector<std::pair<int, double> > tmp = cal_elements(index, basis);
		for(int i = 0; i < tmp.size(); i++){
			elements.push_back(tmp[i]);
		}
	}
	
	
	void observable_base::cal_matrix(Eigen::SparseMatrix<double> &H, computational_basis &basis){
		int d = basis.ref_d();
		int L = basis.ref_L();
		
		std::vector<Eigen::Triplet<double> > elements;
		
		for(int i = 0; i < d; i++){
			std::vector<std::pair<int, double> > tmp = cal_elements(i, basis);
			for(int j = 0; j < tmp.size(); j++){
				elements.push_back(Eigen::Triplet<double>(i, tmp[j].first, tmp[j].second));
			}
		}
		
		H.setFromTriplets(elements.begin(), elements.end());
	}
	
	observable_sum::observable_sum(std::vector<observable_base*> &O){
		for(int i = 0; i < O.size(); i++){
			this->O.push_back(O[i]);
			coef.push_back(1);
		}
	}
	
	observable_sum::observable_sum(observable_base &O1, observable_base &O2){
		this->O.push_back(&O1);
		this->O.push_back(&O2);
		coef.push_back(1);
		coef.push_back(1);
	}
	
	observable_sum::observable_sum(observable_base &O1, double coeff1, observable_base &O2, double coeff2){
		this->O.push_back(&O1);
		this->O.push_back(&O2);
		coef.push_back(coeff1);
		coef.push_back(coeff2);
	}
	
	void observable_sum::append_observable(observable_base &O){
		this->O.push_back(&O);
		coef.push_back(1);
	}
	
	void observable_sum::append_observable(observable_base &O, double coef){
		this->O.push_back(&O);
		this->coef.push_back(coef);
	}
	
	void observable_sum::reset(){
		O.clear();
	}
	
	std::vector<std::pair<int, double> > observable_sum::cal_elements(int index, computational_basis &basis){
		std::vector<std::pair<int, double> > elements;
		int j = 0;
		int j_after;
		for(int i = 0; i < O.size(); i++){
			double A = coef[i];
			O[i]->cal_elements_append(index, basis, elements);
			j_after = elements.size();
			for(int i = j; j < j_after; i++){
				elements[i].second *= A;
			}
			j = j_after;
		}
		
		return elements;
	}
	
	std::pair<bool, double> observable_sum::coeff(int i, int j, computational_basis &basis){
		bool a = false;
		double element = 0;
		for(int i = 0; i < O.size(); i++){
			double A = coef[i];
			std::pair<bool, double> tmp = O[i]->coeff(i, j, basis);
			if(tmp.first == true){
				a = true;
				element += tmp.second * A;
			}
		}
		
		return std::make_pair(a, element);
	}
	
	
	observable_product::observable_product(observable_base &O1, observable_base &O2){
		this->O1 = O1;
		this->O2 = O2;
	}
	
	std::vector<std::pair<int, double> > observable_product::cal_elements(int index, computational_basis &basis){
		int d = basis.ref_d();
		std::vector<std::pair<int, double> > elements;
		
		std::vector<std::pair<int, double> > elements1 = O1.cal_elements(index, basis);
		
		for(int col = 0; col < d; col++){
			double val = 0;
			bool a = false;
			for(int i = 0; i < elements1.size(); i++){
				std::pair<bool, double> element2 = O2.coeff(elements[i].first, col, basis);
				if(element2.first){
					a = true;
					val += element2.second;
				}
				
			}
			
			if(a) elements.push_back(std::make_pair(col, val));
		}
		
		
		return elements;
	}
	
	
	void observable_base_complex::cal_elements_append(int index, computational_basis &basis, std::vector<std::pair<int, std::complex<double> > > &elements){
		
		std::vector<std::pair<int, std::complex<double> > > tmp = cal_elements(index, basis);
		for(int i = 0; i < tmp.size(); i++){
			elements.push_back(tmp[i]);
		}
	}
	
	void observable_base_complex::cal_matrix(Eigen::SparseMatrix<std::complex<double> > &H, computational_basis &basis){
		int d = basis.ref_d();
		int L = basis.ref_L();
		
		std::vector<Eigen::Triplet<std::complex<double> > > elements;
		
		for(int i = 0; i < d; i++){
			std::vector<std::pair<int, std::complex<double> > > tmp = cal_elements(i, basis);
			for(int j = 0; j < tmp.size(); j++){
				elements.push_back(Eigen::Triplet<std::complex<double> >(i, tmp[j].first, tmp[j].second));
			}
		}
		
		H.setFromTriplets(elements.begin(), elements.end());
	}
	
	
	observable_sum_complex::observable_sum_complex(std::vector<observable_base_complex*> &O){
		for(int i = 0; i < O.size(); i++){
			Ocd.push_back(O[i]);
		}
	}
	
	
	
	observable_sum_complex::observable_sum_complex(observable_base_complex &O1, observable_base_complex &O2){
		Ocd.push_back(&O1);
		Ocd.push_back(&O2);
	}
	
	observable_sum_complex::observable_sum_complex(observable_base &O1, observable_base_complex &O2){
		Od.push_back(&O1);
		Ocd.push_back(&O2);
	}
	
	observable_sum_complex::observable_sum_complex(observable_base_complex &O1, observable_base &O2){
		Ocd.push_back(&O1);
		Od.push_back(&O2);
	}
	
	observable_sum_complex::observable_sum_complex(observable_base &O1, observable_base &O2){
		Od.push_back(&O1);
		Od.push_back(&O2);
	}
	
	void observable_sum_complex::append_observable(observable_base_complex &O){
		Ocd.push_back(&O);
	}
	
	void observable_sum_complex::append_observable(observable_base &O){
		Od.push_back(&O);
	}
	
	void observable_sum_complex::reset(){
		Od.clear();
		Ocd.clear();
	}
	
	std::vector<std::pair<int, std::complex<double> > > observable_sum_complex::cal_elements(int index, computational_basis &basis){
		std::vector<std::pair<int, std::complex<double> > > elements;
		
		for(int i = 0; i < Ocd.size(); i++){
			Ocd[i]->cal_elements_append(index, basis, elements);
		}
		
		for(int i = 0; i < Od.size(); i++){
			std::vector<std::pair<int, double> > tmp = Od[i]->cal_elements(index, basis);
			for(int j = 0; j < tmp.size(); j++){
				elements.push_back(std::make_pair(tmp[j].first, std::complex<double>(tmp[j].second, 0)));
			}
		}
		
		return elements;
	}
	
	
}
