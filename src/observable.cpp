#include "observable.h"
#include <iostream>
namespace tylib{
	
	ising_interaction::ising_interaction(double J, int distance, bool isNormalized){
		this->J = J;
		this->distance = distance;
		this->isNormalized = isNormalized;
	}
	
	std::vector<std::pair<int, double> > ising_interaction::cal_elements(int index, computational_basis &basis){
		
		std::vector<std::pair<int, int> > bond = basis.ref_bond(distance);
		
		std::vector<std::pair<int, double> > tmp;
		if(isNormalized){
			tmp.push_back(std::make_pair(index, cal_Sz_interaction(index, basis, bond, -J/basis.ref_L())));
		}else{
			tmp.push_back(std::make_pair(index, cal_Sz_interaction(index, basis, bond, -J)));
		}
		
		return tmp;
		
	}
	
	void ising_interaction::parameters(double J, int distance, bool isNormalized){
		this->J = J;
		this->distance = distance;
		this->isNormalized = isNormalized;
	}
	
	x_ising_interaction::x_ising_interaction(double J, int distance, bool isNormalized){
		this->J = J;
		this->distance = distance;
		this->isNormalized = isNormalized;
	}
	
	std::vector<std::pair<int, double> > x_ising_interaction::cal_elements(int index, computational_basis &basis){
		
		std::vector<std::pair<int, int> > bond = basis.ref_bond(distance);
		
		if(isNormalized){
			return cal_Sx_interaction(index, basis, bond, -J/basis.ref_L());
		}else{
			return cal_Sx_interaction(index, basis, bond, -J);
		}
		
		
	}
	
	void x_ising_interaction::parameters(double J, int distance, bool isNormalized){
		this->J = J;
		this->distance = distance;
		this->isNormalized = isNormalized;
	}
	
	z_magnetic_field::z_magnetic_field(double h, bool isNormalized){
		this->h = h;
	}
	
	std::vector<std::pair<int, double> > z_magnetic_field::cal_elements(int index, computational_basis &basis){
		std::vector<std::pair<int, double> > tmp;
		if(isNormalized){
			tmp.push_back(std::make_pair(index, cal_Sz(index, basis, h/basis.ref_L())));
		}else{
			tmp.push_back(std::make_pair(index, cal_Sz(index, basis, h)));
		}
		
		return tmp;
	}
	
	void z_magnetic_field::parameters(double h, bool isNormalized){
		this->h = h;
		this->isNormalized = isNormalized;
	}
	
	x_magnetic_field::x_magnetic_field(double h, bool isNormalized){
		this->h = h;
		this->isNormalized = isNormalized;
	}
	
	std::vector<std::pair<int, double> > x_magnetic_field::cal_elements(int index, computational_basis &basis){
		if(isNormalized){
			return cal_Sx(index, basis, h/basis.ref_L());
		}else{
			return cal_Sx(index, basis, h);
		}
		
	}
	
	void x_magnetic_field::parameters(double h, bool isNormalized){
		this->h = h;
		this->isNormalized = isNormalized;
	}
	
	XX_model::XX_model(double J, int distance){
		this->J = J;
		this->distance = distance;
	}
	
	std::vector<std::pair<int, double> > XX_model::cal_elements(int index, computational_basis &basis){
		std::vector<std::pair<int, int> > bond = basis.ref_bond(distance);
		
		return cal_hopping(index, basis, bond, -J);
		
	}
	
	void XX_model::parameters(double J, int distance){
		this->J = J;
		this->distance = distance;
	}
	
	HCB_interaction::HCB_interaction(double U, int distance){
		this->U = U;
		this->distance = distance;
	}
	
	std::vector<std::pair<int, double> > HCB_interaction::cal_elements(int index, computational_basis &basis){
		
		std::vector<std::pair<int, int> > bonds = basis.ref_bond(distance);
		
		std::vector<std::pair<int, double> > tmp;
		tmp.push_back(std::make_pair(index, cal_interaction(index, basis, bonds, U)));
		
		return tmp;
		
	}
	
	void HCB_interaction::parameters(double U, int distance){
		this->U = U;
		this->distance = distance;
	}
	
	spin_XXZ_model::spin_XXZ_model(double J, double delta, int distance){
		XX.parameters(J, distance);
		Z.parameters(J*delta, distance);
		
		append_observable(XX);
		append_observable(Z);
	}
	
	void spin_XXZ_model::parameters(double J, double delta, int distance){
		XX.parameters(J, distance);
		Z.parameters(J*delta, distance);
		
		reset();
		
		append_observable(XX);
		append_observable(Z);
	}
	
	XXZ_model::XXZ_model(double J, double delta, int distance){
		XX.parameters(J, distance);
		Z.parameters(J*delta, distance);
		
		append_observable(XX);
		append_observable(Z);
	}
	
	void XXZ_model::parameters(double J, double delta, int distance){
		XX.parameters(J, distance);
		Z.parameters(J*delta, distance);
		
		reset();
		
		append_observable(XX);
		append_observable(Z);
	}
	
	XXZ_with_NNN::XXZ_with_NNN(double J, double delta, double lambda){
		double J1 = J / (1.0 + lambda);
		double J2 = J * lambda / (1.0 + lambda);
		XXZ1.parameters(J1, delta, 1);
		
		append_observable(XXZ1);
		if(lambda != 0){
			XXZ2.parameters(J2, delta, 2);
			append_observable(XXZ2);
		}
	}
	
	void XXZ_with_NNN::parameters(double J, double delta, double lambda){
		
		reset();
		
		double J1 = J / (1.0 + lambda);
		double J2 = J * lambda / (1.0 + lambda);
		XXZ1.parameters(J1, delta, 1);
		
		reset();
		
		append_observable(XXZ1);
		
		if(lambda != 0){
			XXZ2.parameters(J2, delta, 2);
			append_observable(XXZ2);
		}
		
	}
	
	
	Hubbard_model::Hubbard_model(double J, double U){
		this->J = J;
		this->U = U;
	}
	
	std::vector<std::pair<int, double> > Hubbard_model::cal_elements(int index, computational_basis &basis){
		int L = basis.ref_L();
		
		if(L % 2 == 1) return std::vector<std::pair<int, double> >();
		
		//hopping
		std::vector<std::pair<int, int> > bond = basis.ref_bond(2);
		std::vector<std::pair<int, double> > elements = cal_fermion_hopping(index, basis, bond, J);
		
		//interaction
		std::vector<std::pair<int, int> > bond2;
		for(int i = 0; i < L / 2; i++){
			bond2.push_back(std::make_pair(2*i, 2*i+1));
		}
		double val = cal_interaction(index, basis, bond2, U);
		elements.push_back(std::make_pair(index, val));
		
		return elements;
	}
	
	void Hubbard_model::parameters(double J, double U){
		this->J = J;
		this->U = U;
	}
	
	
	MBL::MBL(double W, double J) : spin_XXZ_model(J,J,1){
		this->W = W;
		
		std::random_device rd;
		mt.seed(rd());
		
		std::uniform_real_distribution<double> rand(-W, W);
		
		h.resize(10000);
		for(int i = 0; i < 10000; i++){
			h[i] = rand(mt);
		}
		
	}
	
	MBL::MBL(double W, double J, int seed) : spin_XXZ_model(J,J,1){
		this->W = W;
		
		mt.seed(seed);
		
		std::uniform_real_distribution<double> rand(-W, W);
		
		h.resize(500);
		for(int i = 0; i < 10000; i++){
			h[i] = rand(mt);
		}
		
	}
	
	std::vector<std::pair<int, double> > MBL::cal_elements(int index, computational_basis &basis){
		std::vector<std::pair<int, double> > tmp = spin_XXZ_model::cal_elements(index, basis);
		
		double val = cal_Sz_disorder(index, basis, h);
		tmp.push_back(std::make_pair(index, val));
		
		return tmp;
		
	}
	
	void MBL::reset_disorder(){
		std::uniform_real_distribution<double> rand(-W, W);
		
		for(int i = 0; i < 500; i++){
			h[i] = rand(mt);
		}
	}
	
	transverse_field_ising::transverse_field_ising(double J, double h){
		SzSz.parameters(J, 1);
		Sx.parameters(h);
		
		append_observable(SzSz);
		append_observable(Sx);
	}
	
	void transverse_field_ising::parameters(double J, double h){
		reset();
		
		SzSz.parameters(J, 1);
		Sx.parameters(h);
		
		append_observable(SzSz);
		append_observable(Sx);
	}
	
	// the definition of class n1n2
	std::vector<std::pair<int, double> > n1n2::cal_elements(int index, computational_basis &basis){
		Eigen::VectorXi bits = basis.bits(index);
		
		std::vector<std::pair<int, double> > indexes;
		
		double val;
		if(bits(0) == 1 && bits(1) == 1){
			val = 1;
			indexes.push_back(std::make_pair(index, val));
		}
		
		return indexes;
	}
	
	// the definition of class hopping
	hopping::hopping(int distance){
		this->distance = distance;
	}
	
	std::vector<std::pair<int, double> > hopping::cal_elements(int index, computational_basis &basis){
		std::vector<std::pair<int,int> > bonds = basis.ref_bond(distance);
		
		return cal_hopping(index, basis, bonds, 1.0/basis.ref_L());
	}
	
	
	// the definition of class momentum
	momentum::momentum(int k){
		this->k = k;
	}
	
	std::vector<std::pair<int, std::complex<double> > > momentum::cal_elements(int index, computational_basis &basis){
		int L = basis.ref_L();
		int N = basis.ref_N();
		
		Eigen::VectorXi bits = basis.bits(index);
		
		std::vector<std::pair<int, std::complex<double> > > elements;
		
		std::vector<int> ones, zeros;
		for(int i = 0; i < L; i++){
			if(bits[i] == 0){
				zeros.push_back(i);
			}else{
				ones.push_back(i);
			}
		}
		
		double dtheta = 2 * M_PI * k / L;
		
		for(int i = 0; i < ones.size(); i++){
			for(int j = 0; j < zeros.size(); j++){
				Eigen::VectorXi tmp = bits;
				tmp[ones[i]] = 0;
				tmp[zeros[j]] = 1;
				
				double theta = dtheta * (-ones[i]+zeros[j]);
				std::complex<double> val = std::complex<double>(cos(theta), sin(theta)) / double(L);
				int l = basis.index(tmp);
				elements.push_back(std::make_pair(l, val));
			}
		}
		elements.push_back(std::make_pair(index, double(N) / double(L)));
		
		return elements;
	}
	
	
	std::vector<std::pair<int,double> > momentum_0::cal_elements(int index, computational_basis &basis){
		int L = basis.ref_L();
		int N = basis.ref_N();
		
		Eigen::VectorXi bits = basis.bits(index);
		
		std::vector<std::pair<int,double> > elements;
		
		std::vector<int> ones, zeros;
		for(int i = 0; i < L; i++){
			if(bits[i] == 0){
				zeros.push_back(i);
			}else{
				ones.push_back(i);
			}
		}
		
		
		for(int i = 0; i < ones.size(); i++){
			for(int j = 0; j < zeros.size(); j++){
				Eigen::VectorXi tmp = bits;
				tmp[ones[i]] = 0;
				tmp[zeros[j]] = 1;
				
				double val = 1.0 / double(L);
				int l = basis.index(tmp);
				elements.push_back(std::make_pair(l, val));
			}
		}
		elements.push_back(std::make_pair(index, double(N) / double(L)));
		
		return elements;
	}
	
	std::pair<bool, double> momentum_0::coeff(int i, int j, computational_basis &basis){
		int L = basis.ref_L();
		int N = basis.ref_N();
		
		if(i == j){
			return std::make_pair(true, double(N)/L);
		}else{
			Eigen::VectorXi bits1 = basis.bits(i);
			Eigen::VectorXi bits2 = basis.bits(j);
			int diff = 0;
			for(int l = 0; l < L; l++){
				if(bits1(l) != bits2(l)){
					diff += 1;
					if(diff > 2) break;
				}
			}
			
			if(diff == 2){
				return  std::make_pair(true, 1.0/L);
			}else{
				return std::make_pair(false, 0.0);
			}
			
		}
	}
	
}
