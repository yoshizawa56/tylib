#include "bit.h"

namespace tylib{
	long cal_combination(int L, int N){
		long c = 1;
		for(int i = 1; i <= N; i++){
			c *= L + 1 - i;
			c /= i;
		}
		
		return c;
	}
	
	std::vector<std::pair<int, int> > periodic_1d_bond(int L){
		std::vector<std::pair<int, int> > bond;
		for(int i = 0; i < L-1; i++){
			bond.push_back(std::make_pair(i, i+1));
		}
		bond.push_back(std::make_pair(L-1, 0));
		
		return bond;
		
	}
	
	std::vector<std::pair<int, int> > open_1d_bond(int L){
		std::vector<std::pair<int, int> > bond;
		for(int i = 0; i < L-1; i++){
			bond.push_back(std::make_pair(i, i+1));
		}
		
		return bond;
		
	}
	
	std::vector<std::pair<int, int> > periodic_square_bond(int H, int W){
		std::vector<std::pair<int, int> > bond;
		
		for(int i = 0; i < H; i++){
			for(int j = 0; j < W; j++){
				int index = j + i * W;
				//横方向のボンド
				if(j != W -1){
					bond.push_back(std::make_pair(index, index + 1));
				}
				//縦方向のボンド
				if(i != H -1){
					bond.push_back(std::make_pair(index, index + W));
				}
			}
		}
		
		//周期的境界条件
		if(H > 2){
			for(int i = 0; i < W; i++){
				bond.push_back(std::make_pair(i, (H-1) * W + i));
			}
		}
		if(W > 2){
			for(int i = 0; i < H; i++){
				bond.push_back(std::make_pair(i * W, i * W + W - 1));
			}
		}
		return bond;
		
	}
	
	std::vector<std::pair<int, int> > open_square_bond(int H, int W){
		std::vector<std::pair<int, int> > bond;
		
		for(int i = 0; i < H; i++){
			for(int j = 0; j < W; j++){
				int index = j + i * W;
				//横方向のボンド
				if(j != W -1){
					bond.push_back(std::make_pair(index, index + 1));
				}
				//縦方向のボンド
				if(i != H -1){
					bond.push_back(std::make_pair(index, index + W));
				}
			}
		}
		
		return bond;
		
	}
	
	lattice::lattice(int L, std::vector<std::pair<int, int> > nearest_bond){
		this->L = L;
		
		neighboring_bond.resize(L);
		
		cal_bond(nearest_bond);
		
	}
	
	lattice::lattice(int L, std::string text_name){
		this->L = L;
		
		neighboring_bond.resize(L);
		std::vector<std::pair<int, int> > nearest_bond;
		
		std::ifstream fin(text_name.c_str());
		int i, j;
		while(fin >> i >> j){
			nearest_bond.push_back(std::make_pair(i,j));
		}
		
		cal_bond(nearest_bond);
	}
	
	void lattice::init(int L, std::vector<std::pair<int, int> > nearest_bond){
		this->L = L;
		
		neighboring_bond.clear();
		neighboring_bond.resize(L);
		
		cal_bond(nearest_bond);
	}
	
	void lattice::init(int L, std::string text_name){
		this->L = L;
		
		neighboring_bond.clear();
		neighboring_bond.resize(L);
		
		std::vector<std::pair<int, int> > nearest_bond;
		
		std::ifstream fin(text_name.c_str());
		int i, j;
		while(fin >> i >> j){
			nearest_bond.push_back(std::make_pair(i,j));
		}
		
		cal_bond(nearest_bond);
	}
	
	std::vector<std::pair<int, int> > lattice::ref_bond(int n){
		return neighboring_bond[n-1];
	}
	
	void lattice::cal_bond(std::vector<std::pair<int, int> > &nearest_bond){
		//first < second になるようにソート
		for(int i = 0; i < nearest_bond.size(); i++){
			if(nearest_bond[i].first < nearest_bond[i].second){
				neighboring_bond[0].push_back(std::make_pair(nearest_bond[i].first, nearest_bond[i].second));
			}else{
				neighboring_bond[0].push_back(std::make_pair(nearest_bond[i].second, nearest_bond[i].first));
			}
		}
		
		//重複がでないよう，これまでに出てきたボンドを保存
		std::vector<std::pair<int, int> > all_bond = neighboring_bond[0];
		
		//各サイトがどのサイトと近接しているかのリストを作成
		std::vector<std::vector<int> > nearest_bond_list(L);
		for(int i = 0; i < neighboring_bond[0].size(); i++){
			
			nearest_bond_list[neighboring_bond[0][i].first].push_back(neighboring_bond[0][i].second);
			nearest_bond_list[neighboring_bond[0][i].second].push_back(neighboring_bond[0][i].first);
			
		}
		
		//n+1次近接のサイトを計算
		for(int n = 1; n < L-1; n++){
			//n次近接のボンドリストを作成
			std::vector<std::vector<int> > bond_list(L);
			for(int i = 0; i < neighboring_bond[n-1].size(); i++){
				bond_list[neighboring_bond[n-1][i].first].push_back(neighboring_bond[n-1][i].second);
				bond_list[neighboring_bond[n-1][i].second].push_back(neighboring_bond[n-1][i].first);
				
			}
			
			//n次近接のボンドリストを元にn+1次近接を作成
			for(int i = 0; i < L; i++){
				for(int j = 0; j < bond_list[i].size(); j++){
					int index1 = bond_list[i][j];
					
					for(int k = 0; k < nearest_bond_list[index1].size(); k++){
						int index2 = nearest_bond_list[index1][k];
						
						if(i < index2){
							std::pair<int, int> tmp = std::make_pair(i,index2);
							bool listed = false;
							for(int l = 0; l < all_bond.size(); l++){
								if(all_bond[l] == tmp){
									listed = true;
								}
							}
							if(!listed){
								neighboring_bond[n].push_back(tmp);
								all_bond.push_back(tmp);
							}
							
						}
						
					}
					
				}
			}
		}
	}
	
	computational_basis::computational_basis(int L, int N){
		this->L = L;
		this->N = N;
		
		if(N==-1){			//粒子数非保存
			conserved = false;
			d = pow(2, L);
		}else{
			conserved = true;
			d = cal_combination(L, N);
			ones = N;
			zeros = L - N;
			
			//C.resize(ones+1, std::vector<long>(zeros+1, 0));
			
			C.resize(ones+1, zeros+1);
			
			for(int i = 0; i <= ones; i++){
				for(int j = 0; j <= zeros; j++){
					C(i,j) = cal_combination(i+j, i);
				}
			}
			
		}
		
		Lattice.init(L, periodic_1d_bond(L));
		
	}
	
	computational_basis::computational_basis(int L, int N, std::vector<std::pair<int, int> > nearest_bond){
		this->L = L;
		this->N = N;
		
		if(N==-1){			//粒子数非保存
			conserved = false;
			d = pow(2, L);
		}else{
			conserved = true;
			d = cal_combination(L, N);
			ones = N;
			zeros = L - N;
			
			//C.resize(ones+1, std::vector<long>(zeros+1, 0));
			
			C.resize(ones+1, zeros+1);
			
			for(int i = 0; i <= ones; i++){
				for(int j = 0; j <= zeros; j++){
					C(i,j) = cal_combination(i+j, i);
				}
			}
			
		}
		
		Lattice.init(L, nearest_bond);
		
	}
	
	void computational_basis::init(int L, int N){
		this->L = L;
		this->N = N;
		
		if(N==-1){			//粒子数非保存
			conserved = false;
			d = pow(2, L);
		}else{
			conserved = true;
			d = cal_combination(L, N);
			ones = N;
			zeros = L - N;
			
			//C.resize(ones+1, std::vector<long>(zeros+1, 0));
			C.resize(ones+1, zeros+1);
			
			for(int i = 0; i <= ones; i++){
				for(int j = 0; j <= zeros; j++){
					C(i,j) = cal_combination(i+j, i);
				}
			}
			
		}
		
		Lattice.init(L, periodic_1d_bond(L));
		
	}
	
	void computational_basis::init(int L, int N, std::vector<std::pair<int, int> > nearest_bond){
		this->L = L;
		this->N = N;
		
		if(N==-1){			//粒子数非保存
			conserved = false;
			d = pow(2, L);
		}else{
			conserved = true;
			d = cal_combination(L, N);
			ones = N;
			zeros = L - N;
			
			//C.resize(ones+1, std::vector<long>(zeros+1, 0));
			
			C.resize(ones+1, zeros+1);
			
			for(int i = 0; i <= ones; i++){
				for(int j = 0; j <= zeros; j++){
					C(i,j) = cal_combination(i+j, i);
				}
			}
			
		}
		
		Lattice.init(L, nearest_bond);
		
	}
	
	long computational_basis::index(Eigen::VectorXi& bits){
		long index = 0;
		
		if(!conserved){
			for(int i = 0; i < L; i++){
				if(bits[i] == 1) index += pow(2, i);
			}
			
			return index;
		}
		
		int zeros_tmp, ones_tmp;
		zeros_tmp = zeros;
		ones_tmp = ones;
		
		
		
		for(int i = 0; i < L; i++){
			if(bits(i) == 0){
				index += C(ones_tmp-1,zeros_tmp);
				zeros_tmp--;
			}else{
				ones_tmp--;
			}
			
			if(zeros_tmp == 0){
				break;
			}
			
			if(ones_tmp == 0){
				break;
			}
		}
		return index;
	}
	
	Eigen::VectorXi computational_basis::bits(long index){
		Eigen::VectorXi bits(L);
		
		if(!conserved){
			for(int i = 0; i < L; i++){
				bits[i] = index % 2;
				
				index /= 2;
			}
			return bits;
			
		}
		
		int zeros_tmp, ones_tmp;
		zeros_tmp = zeros;
		ones_tmp = ones;
		
		for(int i = 0; i< L; i++){
			if(zeros_tmp == 0){
				bits(i) = 1;
				continue;
			}
			
			if(ones_tmp == 0){
				bits(i) = 0;
				continue;
			}
			
			if(index >= C(ones_tmp-1,zeros_tmp)){
				index -= C(ones_tmp-1,zeros_tmp);
				bits(i) = 0;
				zeros_tmp--;
			}else{
				bits(i) = 1;
				ones_tmp--;
			}
		}
		
		return bits;
	}
	
	std::pair<Eigen::VectorXi, Eigen::VectorXi> computational_basis::fermion_bits(long index){
		Eigen::VectorXi bits(L);
		Eigen::VectorXi count(L);
		
		count(0) = 0;
		if(!conserved){
			bits(0) = index % 2;
			index /= 2;
			
			for(int i = 1; i < L; i++){
				bits(i) = index % 2;
				count(i) = bits(i-1) ^ count(i-1);
				
				index /= 2;
			}
			return std::make_pair(bits, count);
			
		}
		
		int zeros_tmp, ones_tmp;
		zeros_tmp = zeros;
		ones_tmp = ones;
		
		if(index >= C(ones_tmp-1,zeros_tmp)){
				index -= C(ones_tmp-1,zeros_tmp);
				bits(0) = 0;
				zeros_tmp--;
			}else{
				bits(0) = 1;
				ones_tmp--;
			}
		
		for(int i = 1; i< L; i++){
			if(zeros_tmp == 0){
				bits(i) = 1;
				count(i) = bits(i-1) ^ count(i-1);
				continue;
			}
			
			if(ones_tmp == 0){
				bits(i) = 0;
				count(i) = bits(i-1) ^ count(i-1);
				continue;
			}
			
			if(index >= C(ones_tmp-1,zeros_tmp)){
				index -= C(ones_tmp-1,zeros_tmp);
				bits(i) = 0;
				zeros_tmp--;
			}else{
				bits(i) = 1;
				ones_tmp--;
			}
			
			count(i) = bits(i-1) ^ count(i-1);
		}
		
		return std::make_pair(bits, count);
	}
	
	int computational_basis::ref_L(){
		
		return L;
	}
	
	long computational_basis::ref_d(){
		
		return d;
	}
	
	int computational_basis::ref_N(){
		
		return N;
	}
	
	bool computational_basis::isConserved(){
		
		return conserved;
	}
	
	lattice computational_basis::ref_lattice(){
		
		return Lattice;
	}
	
	std::vector<std::pair<int, int> > computational_basis::ref_bond(int n){
		
		return Lattice.ref_bond(n);
	}

}