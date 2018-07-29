#include "ss.h"
#include <iostream>

namespace tylib{
	
	komega_bicg::komega_bicg(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold){
		
		//define parameters
		d = v.size();
		nz = z.size();
		this->itermax = itermax;
		this->threshold = threshold;
		
		//allocate v12, v2, v14, v4, x and z
		v12 = new std::complex<double>[d];
		v2 = new std::complex<double>[d];
		v14 = new std::complex<double>[d];
		v4 = new std::complex<double>[d];
		x = new std::complex<double>[d*nz];
		this->z = new std::complex<double>[nz];
		
		_status[0] = 0;
		_status[1] = 0;
		_status[2] = 0;
		
		//set z values
		for(int i = 0; i < nz; i++){
			this->z[i] = z(i);
		}
		
		for(int i = 0; i < d; i++){
			v2[i] = v(i);
			v4[i] = v(i);
		}
		
		komega_bicg_init(&d, &d, &nz, x, this->z, &itermax, &threshold);
		
	}
	
	komega_bicg::~komega_bicg(){
		komega_bicg_finalize();
		delete[] x;
		delete[] z;
		delete[] v12;
		delete[] v2;
		delete[] v14;
		delete[] v4;
	}
	
	void komega_bicg::solve(Eigen::SparseMatrix<std::complex<double> > &H){
		
		Eigen::VectorXcd Hr1(d), Hr2(d);
		for(int i = 0; i < d; i++){
			Hr1(i) = v2[i];
		}
		Hr2 = Hr1;
		
		while( !isConverged() && status() == 0){
			Hr1 = H * Hr1;
			Hr2 = H * Hr2;
			update(Hr1, Hr2);
		}
		
		
	}
	
	void komega_bicg::update(Eigen::VectorXcd& Hr1, Eigen::VectorXcd& Hr2){
		for(int i = 0; i < d; i++){
			v12[i] = Hr1(i);
			v14[i] = Hr2(i);
		}
		
		komega_bicg_update(v12, v2, v14, v4, x, v2, _status);
		
		//Hr に残差ベクトルを代入．（次のステップではこのベクトルにハミルトニアンをかければよい）
		for(int i = 0; i < d; i++){
			Hr1(i) = v2[i];
			Hr2(i) = v4[i];
		}
		
	}
	
	void komega_bicg::init(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold){
		
		//define parameters
		d = v.size();
		nz = z.size();
		this->itermax = itermax;
		this->threshold = threshold;
		
		//allocate v12, v2, v14, v4, x and z
		v12 = new std::complex<double>[d];
		v2 = new std::complex<double>[d];
		v14 = new std::complex<double>[d];
		v4 = new std::complex<double>[d];
		x = new std::complex<double>[d*nz];
		this->z = new std::complex<double>[nz];
		
		_status[0] = 0;
		_status[1] = 0;
		_status[2] = 0;
		
		//set z values
		for(int i = 0; i < nz; i++){
			this->z[i] = z(i);
		}
		
		for(int i = 0; i < d; i++){
			v2[i] = v(i);
			v4[i] = v(i);
		}
		
		komega_bicg_init(&d, &d, &nz, x, this->z, &itermax, &threshold);
		
	}
	
	Eigen::VectorXcd komega_bicg::getSolution(int i){
		Eigen::VectorXcd xi(d);
		int j_start = i * d;
		for(int j = 0; j < d; j++){
			xi(j) = x[j+j_start];
		}
		
		return xi;
	}
	
	int komega_bicg::iterations(){
		return abs(_status[0]);
		
	}
	
	bool komega_bicg::isConverged(){
		if(_status[0] < 0 && _status[1] == 0){
			return true;
		}else{
			return false;
		}
		
	}
	
	int komega_bicg::status(){
		return _status[1];
	}
	
	double komega_bicg::residue(){
		Eigen::VectorXcd tmp(d);
		for(int i = 0; i < d; i++){
			tmp(i) = v12[i];
		}
		return tmp.norm();
	}
	
	int komega_bicg::seed(){
		return _status[2];
	}
	
	komega_cocg::komega_cocg(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold){
		
		//define parameters
		d = v.size();
		nz = z.size();
		this->itermax = itermax;
		this->threshold = threshold;
		
		//allocate v12, v2, v14, v4, x and z
		v12 = new std::complex<double>[d];
		v2 = new std::complex<double>[d];
		x = new std::complex<double>[d*nz];
		this->z = new std::complex<double>[nz];
		
		_status[0] = 0;
		_status[1] = 0;
		_status[2] = 0;
		
		//set z values
		for(int i = 0; i < nz; i++){
			this->z[i] = z(i);
		}
		
		for(int i = 0; i < d; i++){
			v2[i] = v(i);
		}
		
		komega_cocg_init(&d, &d, &nz, x, this->z, &itermax, &threshold);
		
	}
	
	komega_cocg::~komega_cocg(){
		komega_cocg_finalize();
		delete[] x;
		delete[] z;
		delete[] v12;
		delete[] v2;
	}
	
	void komega_cocg::solve(Eigen::SparseMatrix<double> &H){
		
		Eigen::VectorXcd Hr(d);
		for(int i = 0; i < d; i++){
			Hr(i) = v2[i];
		}
		
		while( !isConverged() && status() == 0){
			Hr = H * Hr;
			update(Hr);
		}
		
	}
	
	void komega_cocg::update(Eigen::VectorXcd& Hr){
		for(int i = 0; i < d; i++){
			v12[i] = Hr(i);
		}
		
		komega_cocg_update(v12, v2, x, v2, _status);
		
		//Hr に残差ベクトルを代入．（次のステップではこのベクトルにハミルトニアンをかければよい）
		for(int i = 0; i < d; i++){
			Hr(i) = v2[i];
		}
		
	}
	
	void komega_cocg::init(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold){
		
		//define parameters
		d = v.size();
		nz = z.size();
		this->itermax = itermax;
		this->threshold = threshold;
		
		//allocate v12, v2, v14, v4, x and z
		v12 = new std::complex<double>[d];
		v2 = new std::complex<double>[d];
		x = new std::complex<double>[d*nz];
		this->z = new std::complex<double>[nz];
		
		_status[0] = 0;
		_status[1] = 0;
		_status[2] = 0;
		
		//set z values
		for(int i = 0; i < nz; i++){
			this->z[i] = z(i);
		}
		
		for(int i = 0; i < d; i++){
			v2[i] = v(i);
		}
		
		komega_cocg_init(&d, &d, &nz, x, this->z, &itermax, &threshold);
		
	}
	
	Eigen::VectorXcd komega_cocg::getSolution(int i){
		Eigen::VectorXcd xi(d);
		int j_start = i * d;
		for(int j = 0; j < d; j++){
			xi(j) = x[j+j_start];
		}
		
		return xi;
	}
	
	int komega_cocg::iterations(){
		return abs(_status[0]);
		
	}
	
	bool komega_cocg::isConverged(){
		if(_status[0] < 0 && _status[1] == 0){
			return true;
		}else{
			return false;
		}
		
	}
	
	int komega_cocg::status(){
		return _status[1];
	}
	
	double komega_cocg::residue(){
		Eigen::VectorXcd tmp(d);
		for(int i = 0; i < d; i++){
			tmp(i) = v12[i];
		}
		return tmp.norm();
	}
	
	int komega_cocg::seed(){
		return _status[2];
	}
	
	
	random_vector::random_vector(){
		std::random_device rd;
		mt.seed(rd());
	}
	
	Eigen::VectorXcd random_vector::gauss(int d){
		Eigen::VectorXcd v(d);
		
		std::normal_distribution<> gauss(0.0,1.0);
		
		for(int i = 0; i < d; i++){
			v(i) = std::complex<double>(gauss(mt), 0);
			
		}
		
		v.normalize();
		
		return v;
	}
	
	Eigen::VectorXcd random_vector::binary(int d){
		Eigen::VectorXcd v(d);
		
		std::uniform_real_distribution<> uni(0.0,1.0);
		for(int i = 0; i < d; i++){
			if(uni(mt) > 0.5){
				v(i) = 1;
			}else{
				v(i) = -1;
			}
			
		}
		
		return v;
	}
	
	ss_method::ss_method(Eigen::SparseMatrix<double> &H){
		isComplex = false;
		Hd = H;
		d = H.cols();
		n_threads = 1;
		memory_lim = 32;
		itermax = 0;
		threshold = 1e-10;
		contour = 0;
		n = 20;
		contour_para = 0.1;
		m = 4;
		step = 1;
		kappa = 2;
		cutoff = 1e-10;
		
	}
	
	ss_method::ss_method(Eigen::SparseMatrix<std::complex<double> > &H){
		isComplex = true;
		Hcd = H;
		d = H.cols();
		n_threads = 1;
		memory_lim = 32;
		itermax = 0;
		threshold = 1e-10;
		contour = 0;
		n = 20;
		contour_para = 0.1;
		m = 4;
		step = 1;
		kappa = 2;
		cutoff = 1e-10;
		
	}
	
	void ss_method::eigenstate(Eigen::MatrixXcd &evec, Eigen::VectorXd &E, double gamma, double rho){
		
		double begin = omp_get_wtime();
		
		//状態数の見積もり
		int N = kappa * num_states(gamma, rho) / m + 1;
		
		double end = omp_get_wtime();
		
		std::cout << "必要ベクトル数 : " << N * m << std::endl;
		if(n_threads < 20){
			std::cout << "予想計算時間 : " << (end - begin) * N / 20 << "秒" << std::endl;
		}else{
			std::cout << "予想計算時間 : " << (end - begin) * N / n_threads << "秒" << std::endl;
		}
		
		//必要なメモリの確認
		double memory = (d * m * N * 3) * 16 / 1e9;
		std::cout << "予想メモリ使用量 : " << memory << " GB" << std::endl;
		if(memory > memory_lim){
			std::cout << "メモリーの使用が上限を超えます" << std::endl;
			std::cout << "予想使用量 : " << memory << " GB" << std::endl;
			exit(1);
		}
		std::cout << std::endl;
		
		
		Eigen::MatrixXcd S(d, m * N);
		random_vector Rv;
		
		//射影の計算
		#pragma omp parallel for num_threads(n_threads)
		for(int l = 0; l < N; l++){
			Eigen::VectorXcd v = Rv.gauss(d);
			S.block(0, l*m, d, m) = projection(v, gamma, rho, m);
		}
		
		projection_to_eigenstate(evec, E, S);
		
	}
	
	Eigen::MatrixXcd ss_method::projection(Eigen::VectorXcd &v, double gamma, double rho, int m){
		//積分経路の形成
		Eigen::VectorXcd z = make_contour(gamma, rho);
		
		Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(d, m);
		if(isComplex){
			komega_bicg bicg(v, z, itermax, threshold);
			bicg.solve(Hcd);
			
			for(int j = 0; j < n; j++){
				
				std::complex<double> value = prefactor(gamma, rho, j);
				
				for(int i = 0; i < m; i++){
					S.col(i) += value * pow(z(j), i * step) * bicg.getSolution(j);
				}
				
			}
			
		}else{
			komega_cocg cocg(v, z, itermax, threshold);
			cocg.solve(Hd);
			
			for(int j = 0; j < n; j++){
				
				std::complex<double> value = prefactor(gamma, rho, j);
				
				for(int i = 0; i < m; i++){
					S.col(i) += value * pow(z(j), i * step) * cocg.getSolution(j);
				}
				
			}
		}
		
		return S;
		
	}
	
	int ss_method::num_states(double gamma, double rho){
		
		random_vector Rv;
		
		int N;
		if(n_threads < 20){
			N = 20;
		}else{
			N = n_threads;
		}
		
		Eigen::VectorXcd trace(N);
		#pragma omp parallel for num_threads(n_threads)
		for(int l = 0; l < N; l++){
			Eigen::VectorXcd v = Rv.binary(d);
			Eigen::MatrixXcd S = projection(v, gamma, rho, 1);
			trace(l) = v.adjoint() * S.col(0);
			
		}
		
		return int(real(trace.sum()) / N);
		
	}
	
	void ss_method::projection_to_eigenstate(Eigen::MatrixXcd &evec, Eigen::VectorXd &E, Eigen::MatrixXcd &S){
		
		//Eigen::JacobiSVD<Eigen::MatrixXcd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::JacobiSVD<Eigen::MatrixXcd>* svd = new Eigen::JacobiSVD<Eigen::MatrixXcd>(S, Eigen::ComputeThinU);
		S.resize(0,0);
		
		int rank;
		{
			int ntmp = svd->singularValues().size();
			rank = ntmp;
			double sigma0 = svd->singularValues()(0);
			for(int i = 1; i < ntmp; i++){
				if(svd->singularValues()(i)/sigma0 < cutoff){
					rank = i;
					break;
				}
			}
		}
		
		Eigen::MatrixXcd Qt(d,rank);
		for(int j = 0; j < rank; j++){
			for(int i = 0; i < d; i++){
				Qt(i,j) = svd->matrixU()(i,j);
			}
		}
		
		delete svd;
		
		Eigen::MatrixXcd Ht(rank,rank);
		if(isComplex){
			Ht = Qt.adjoint() * (Hcd * Qt);
		}else{
			Ht = Qt.adjoint() * (Hd * Qt);
		}
		
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> est(Ht);
		
		evec.resize(d, rank);
		E.resize(rank);
		
		for(int i = 0; i < rank; i++){
			E(i) = std::real(est.eigenvalues()(i));
			
			evec.col(i) = Qt * est.eigenvectors().col(i);
			evec.col(i).normalize();
		}
		
	}
	
	
	void ss_method::precision_check(Eigen::MatrixXcd &evec, Eigen::VectorXd &E){
		int N = evec.cols();
		
		for(int i = 0; i < N; i++){
			Eigen::VectorXcd tmp;
			if(isComplex){
				tmp = Hcd * evec.col(i);
			}else{
				tmp = Hd * evec.col(i);
			}
			
			double residue = abs(abs(E(i)) - tmp.norm());
			std::cout << E(i) << " 	" << residue << std::endl;
		}
		
		std::cout << std::endl;
	}
	
	void ss_method::enable_openmp(int n_threads){
		this->n_threads = n_threads;
	}
	
	void ss_method::memory_setting(int memory){
		memory_lim = memory;
	}
	
	void ss_method::itermax_setting(int itermax){
		this->itermax = itermax;
	}
	
	void ss_method::threshold_setting(int threshold){
		this->threshold = threshold;
	}
	
	void ss_method::contour_setting(int contour, int n, double contour_para){
		this->contour = contour;
		this->n = n;
		this->contour_para = contour_para;
	}
	
	void ss_method::copy_setting(int m, int step){
		this->m = m;
		this->step = step;
	}
	
	void ss_method::kappa_setting(double kappa){
		this->kappa = kappa;
	}
	
	void ss_method::cutoff_setting(double cutoff){
		this->cutoff = cutoff;
	}
	
	std::complex<double> ss_method::prefactor(double gamma, double rho, int j){
		
		std::complex<double> value;
		if(contour == 0){			//楕円
			double theta = 2.0 * M_PI * (j + 0.5) / n;
			value = rho * std::complex<double>(contour_para * cos(theta), sin(theta)) / double(n);
		}else{						//長方形
			if(j == 0){
				value = std::complex<double>(contour_para, -2.0*rho/(n-2));
			}else if(j == n/2-1){
				value = std::complex<double>(-contour_para,-2.0*rho/(n-2));
			}else if(j == n/2){
				value = std::complex<double>(-contour_para, 2.0*rho/(n-2));
			}else if(j == n-1){
				value = std::complex<double>(contour_para, 2*rho/(n-2));
			}else if(j < n/2){
				value = std::complex<double>(0, -4*rho) / double(n-2);
			}else{
				value = std::complex<double>(0, 4*rho) / double(n-2);
			}
			
			value /= -2*M_PI;
		}
		
		return value;
		
	}
	
	Eigen::VectorXcd ss_method::make_contour(double gamma, double rho){
		Eigen::VectorXcd z(n);
		
		if(contour == 0){			//楕円
			for(int i = 0; i < n; i++){
				double theta = 2.0 * M_PI * (i + 0.5) / n;
				z(i) = gamma + rho * std::complex<double>(cos(theta), contour_para * sin(theta));
			}
		}else{
			for(int i = 0; i < n/2; i++){
				z(i) = std::complex<double>(gamma + rho * (4*i-n+2)/(n-2), contour_para);
				z(n-1-i) = std::complex<double>(gamma + rho * (4*i-n+2)/(n-2), -contour_para);
			}
			
		}
		
		return z;
	}
	
}