#include "translation.h"
#include <iostream>

namespace tylib{
	
	
	translational_basis::translational_basis(computational_basis &c_basis){
		
		
		L = c_basis.ref_L();
		N = c_basis.ref_N();
		d = c_basis.ref_d();
		int fill = L / N;
		
		//std::vector<int> period;
		//存在する周期を列挙
		if(c_basis.isConserved()){		//粒子数保存の場合
			for(int i = 1; i <= N; i++){
				if(N%i == 0){
					period.push_back(fill* i);
				}
			}
		}else{							//粒子数非保存の場合
			for(int i = 1; i <= L; i++){
				if(L % i == 0){
					period.push_back(i);
				}
			}
			
		}
		
		p = period.size();
		
		
		//各固有値で固有ベクトルが存在する周期を列挙
		r_period.resize(L);
		for(int i = 0; i < L; i++){
			for(int j = 0; j < p; j++){
				int t = period[j];
				
				if(i % (L/t) == 0){
					r_period[i].push_back(j);
				}
			}
		}
		/*
		for(int i = 0; i < L; i++){
			std::cout << "r=" << i << std::endl;
			for(int j = 0; j < p; j++){
				std::cout << r_period[i][j] << std::endl;
			
			}
		}
		*/
		
		table.resize(d);
		for(int i = 0; i < d; i++){
			table[i].resize(3);
		}
		
		//周期ごとに計算基底のグループを分類
		basis.resize(p);
		D_period = Eigen::VectorXi::Zero(period.size());
		{
			std::vector<int> exclusion(d, 0);
			for(int i = 0; i < p; i++){
				Eigen::VectorXi bits(L);
				int t = period[i];
				int xN;
				if(c_basis.isConserved()){
					xN = N / (L / t);
				}else{
					xN = -1;
				}
				int group_index = 0;
				
				computational_basis xc_basis(t, xN);
				
				int xd = xc_basis.ref_d();
				for(int j = 0; j < xd; j++){
					Eigen::VectorXi xbits = xc_basis.bits(j);
					Eigen::VectorXi bits(L);
					
					for(int l = 0; l < L/t; l++){
						for(int k = 0; k < t; k++){
							bits(l*t+k) = xbits(k);
						}
					}
					int index = c_basis.index(bits);
					
					if(exclusion[index] == 0){
						exclusion[index] = 1;
						
						//tableに情報を記録
						table[index][0] = i;
						table[index][1] = group_index;
						table[index][2] = 0;
						
						std::vector<int> index_tmp(t);
						index_tmp[0] = index;
						
						std::vector<int> tmp;
						tmp.push_back(index);
						exclusion.push_back(index);
						//並進操作で系列を作成
						for(int l = 1; l < t; l++){
							translation(bits);
							index = c_basis.index(bits);
							table[index][0] = i;
							table[index][1] = group_index;
							table[index][2] = l;
							index_tmp[l] = index;
							exclusion[index] = 1;
							tmp.push_back(index);
						}
						//basisに情報を記録
						basis[i].push_back(index_tmp);
						group_index++;
						D_period(i)++;
					}
					
				}
			}
		}
		
		//各固有値に属する固有ベクトルの本数
		D = Eigen::VectorXi::Zero(L);
		for(int i = 0; i < L; i++){
			for(int j = 0; j < p; j++){
				for(int k = 0; k < r_period[i].size(); k++){
					if(r_period[i][k] == j){
						D(i) += D_period[j];
					}
				}
				
			}
		}
		
		t_index_init.resize(L);
		for(int r = 0; r < L; r++){
			t_index_init[r].resize(p, 0);
			for(int i = 0; i < r_period[r].size(); i++){
				for(int j = 0; j < i; j++){
					t_index_init[r][r_period[r][i]] += D_period[r_period[r][j]];
				}
				
			}
		}
		
	}
	
	
	void translational_basis::cal_observable(Eigen::SparseMatrix<std::complex<double> > &obs, observable_base &operation, int r, bool isTranslational){
		computational_basis basis(L, N);
		
		std::vector<Eigen::Triplet<std::complex<double> > > elements;
		
		if(isTranslational){
			
			
			
			for(int i = 0; i < D(r); i++){
				std::vector<int> phi;
				cal_basis(phi, r, i);
				
				int index = phi[0];
				int t_phi = period[table[index][0]];
				
				std::vector<std::pair<int, double> > indexes = operation.cal_elements(index, basis);
				
				for(int j = 0; j < indexes.size(); j++){
					int t_index = cal_t_index(indexes[j].first, r);
					if(t_index == -1) continue;
					int t_psi = period[table[indexes[j].first][0]];
					int delta = table[indexes[j].first][2];		//固有ベクトル中で何番目か
					double theta = -2.0 * M_PI * r * delta / L;
					std::complex<double> val = indexes[j].second * std::complex<double>(cos(theta), sin(theta)) * sqrt(double(t_phi)/double(t_psi));
					
					elements.push_back(Eigen::Triplet<std::complex<double> >(i, t_index, val));
				}
				
			}
			
		}else{
			
			
			for(int i = 0; i < d; i++){
				int t_i = table[i][0];
				int index_i = table[i][1] + t_index_init[r][t_i];
				int theta_i = table[i][2];
				
				
				if((r % (L/period[t_i])) != 0){
					continue;
				}
				
				//計算基底|i>に対して非ゼロの成分を持つ<j|のリストを取得
				std::vector<std::pair<int, double> > tmp = operation.cal_elements(i, basis);
				for(int j = 0; j < tmp.size(); j++){
					int t_j = table[tmp[j].first][0];
					int index_j = table[tmp[j].first][1] + t_index_init[r][t_j];
					int theta_j = table[tmp[j].first][2];
					if((r % (L/period[t_j])) != 0){
						continue;
					}
					double theta = 2 * M_PI * (double(theta_i)*r/L  - double(theta_j)*r/L);
					std::complex<double> val;
					val = tmp[j].second / sqrt(period[t_i] * period[t_j]) * std::complex<double>(cos(theta), sin(theta));
					elements.push_back(Eigen::Triplet<std::complex<double> >(index_i, index_j, val));
				}
			}
			
		}
		
		obs.resize(D(r), D(r));
		obs.setFromTriplets(elements.begin(), elements.end());
		
	}
	
	void translational_basis::cal_observable(Eigen::SparseMatrix<std::complex<double> > &obs, observable_base_complex &operation, int r, bool isTranslational){
		computational_basis basis(L, N);
		
		std::vector<Eigen::Triplet<std::complex<double> > > elements;
		
		if(isTranslational){
			
			
			
			for(int i = 0; i < D(r); i++){
				std::vector<int> phi;
				cal_basis(phi, r, i);
				
				int index = phi[0];
				int t_phi = period[table[index][0]];
				
				std::vector<std::pair<int, std::complex<double> > > indexes = operation.cal_elements(index, basis);
				
				for(int j = 0; j < indexes.size(); j++){
					int t_index = cal_t_index(indexes[j].first, r);
					if(t_index == -1) continue;
					int t_psi = period[table[indexes[j].first][0]];
					int delta = table[indexes[j].first][2];		//固有ベクトル中で何番目か
					double theta = -2.0 * M_PI * r * delta / L;
					std::complex<double> val = indexes[j].second * std::complex<double>(cos(theta), sin(theta)) * sqrt(double(t_phi)/double(t_psi));
					
					elements.push_back(Eigen::Triplet<std::complex<double> >(i, t_index, val));
				}
				
			}
			
		}else{
			
			
			for(int i = 0; i < d; i++){
				int t_i = table[i][0];
				int index_i = table[i][1] + t_index_init[r][t_i];
				int theta_i = table[i][2];
				
				
				if((r % (L/period[t_i])) != 0){
					
					continue;
				}
				
				//計算基底|i>に対して非ゼロの成分を持つ<j|のリストを取得
				std::vector<std::pair<int, std::complex<double> > > tmp = operation.cal_elements(i, basis);
				for(int j = 0; j < tmp.size(); j++){
					int t_j = table[tmp[j].first][0];
					int index_j = table[tmp[j].first][1] + t_index_init[r][t_j];
					int theta_j = table[tmp[j].first][2];
					
					if((r % (L/period[t_j])) != 0){
						continue;
					}
					
					double theta = 2 * M_PI * (double(theta_i)*r/L  - double(theta_j)*r/L);
					std::complex<double> val;
					val = tmp[j].second / sqrt(period[t_i] * period[t_j]) * std::complex<double>(cos(theta), sin(theta));
					elements.push_back(Eigen::Triplet<std::complex<double> >(index_i, index_j, val));
				}
			}
			
		}
		
		obs.resize(D(r), D(r));
		obs.setFromTriplets(elements.begin(), elements.end());
	}
	
	Eigen::VectorXcd translational_basis::eigenvector(int r, int index){
		Eigen::VectorXcd v(D(r));
		
		std::vector<int> c_indexes;
		
		cal_basis(c_indexes, r, index);
		int t = c_indexes.size();
		
		v = Eigen::VectorXcd::Zero(d);
		for(int i = 0; i < t; i++){
			double theta = 2 * M_PI * i * r / L;
			v(c_indexes[i]) = std::complex<double>(cos(theta), sin(theta)) / sqrt(t);
		}
		
		return v;
		
	}
	
	Eigen::VectorXcd translational_basis::c_to_t(Eigen::VectorXd &vc, int r){
		
		Eigen::VectorXcd vt = Eigen::VectorXcd::Zero(D(r));
		
		for(int i = 0; i < d; i++){
			int t_index = cal_t_index(i, r);
			if(t_index != -1){
				int p = period[table[i][0]];
				double theta = 2.0 * M_PI * r * table[i][2] / L;
				std::complex<double> factor = std::complex<double>(cos(theta), sin(theta)) / sqrt(p);
				vt(t_index) += vc(i) * factor;
			}
			
		}
		
		
		return vt;
	}
	
	Eigen::VectorXcd translational_basis::c_to_t(Eigen::VectorXcd &vc, int r){
		
		Eigen::VectorXcd vt = Eigen::VectorXcd::Zero(D(r));
		
		for(int i = 0; i < d; i++){
			int t_index = cal_t_index(i, r);
			if(t_index != -1){
				int p = period[table[i][0]];
				double theta = 2.0 * M_PI * r * table[i][2] / L;
				std::complex<double> factor = std::complex<double>(cos(theta), sin(theta)) / sqrt(p);
				vt(t_index) += vc(i) * factor;
			}
			
		}
		
		
		return vt;
	}
	
	Eigen::VectorXcd translational_basis::t_to_c(Eigen::VectorXcd &vt, int r){
		
		Eigen::VectorXcd vc = Eigen::VectorXcd::Zero(d);
		
		for(int i = 0; i < D(r); i++){
			Eigen::VectorXcd tmp = eigenvector(r, i);
			vc += vt(i) * tmp;
		}
		
		
		return vc;
	}
	
	std::pair<int, std::complex<double> > translational_basis::c_basis_to_t_basis(int c_index, int r){
		int t_index = cal_t_index(c_index, r);
		std::complex<double> factor;
		if(t_index != -1){
			int p = period[table[c_index][0]];
			double theta = 2.0 * M_PI * r * table[c_index][2] / L;
			factor = std::complex<double>(cos(theta), sin(theta)) / sqrt(p);
			//vt(t_index) += vc(i) * factor;
		}else{
			factor = 0;
		}
		
		return std::make_pair(t_index, factor);
	}
	
	std::vector<int> translational_basis::ref_period(){
		return period;
	}
	
	std::vector<int> translational_basis::ref_rperiod(int r){
		return r_period[r];
	}
	
	int translational_basis::ref_d(){
		return d;
	}
	
	int translational_basis::ref_D(int r){
		return D(r);
	}
	
	int translational_basis::ref_D_period(int i){
		return D_period(i);
	}
	
	int translational_basis::cal_t_index(int i, int r){
		
		
		int index = t_index_init[r][table[i][0]];
		
		//固有値rに属さないベクトルの場合は-1を返す
		if(table[i][0] != r_period[r][0] && index == 0){
			return -1;
		}
		
		index += table[i][1];
		
		return index;
		
	}
	
	void translational_basis::cal_basis(std::vector<int> &indexes, int r, int index){
		int tmp = 0;
		
		
		if(index < D(r)){
		
			int t;
			for(int j = 0; j < r_period[r].size(); j++){
				if(tmp + D_period(r_period[r][j]) <= index){
					tmp += D_period(r_period[r][j]);
				}else{
					t = r_period[r][j];
					break;
				}
			}
			tmp = index - tmp;
			
			for(int i = 0; i < period[t]; i++){
				indexes.push_back(basis[t][tmp][i]);
			}
		}
	}
	
	
	void translation(Eigen::VectorXi& bits){
		int L = bits.size();
		int tmp1 = bits(0);
		int tmp2;
		bits(0) = bits(L-1);
		for(int i = 1; i < L; i++){
			tmp2 = bits(i);
			bits(i) = tmp1;
			tmp1 = tmp2;
		}
	}
	
}
