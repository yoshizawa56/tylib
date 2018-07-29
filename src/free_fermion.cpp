#include "free_fermion.h"

namespace tylib{
	
	void XX(Eigen::MatrixXd &H, int N, double J, double h, bool isPeriodic){
		int L = H.rows();
		
		H = Eigen::MatrixXd::Zero(L,L);
		for(int i = 0; i < L-1; i++){
			H(i,i) = h;
			H(i, i+1) = -J;
			H(i+1, i) = -J;
		}
		H(L-1, L-1) = h;
		
		if(isPeriodic){
			if(N % 2 == 0){
				H(0, L-1) = J;
				H(L-1, 0) = J;
			}else{
				H(0, L-1) = -J;
				H(L-1, 0) = -J;
			}
		}
	}
	
	
	
	template <typename T> Slater_determinant<T>::Slater_determinant(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &phi){
		L = phi.rows();
		N = phi.cols();
		
		this->phi.resize(L, N);
		this->phi = phi;
		
	}
	
	template <typename T> Slater_determinant<T>::Slater_determinant(computational_basis basis, int index){
		L = basis.ref_L();
		N = basis.ref_N();
		
		phi = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(L,N);
		
		Eigen::VectorXi bits = basis.bits(index);
		int n = 0;
		for(int i = 0; i < L; i++){
			if(bits(i) == 1){
				phi(n, i) = 1;
				n++;
			}
		}
	}
	
	template <typename T> Slater_determinant<T>::Slater_determinant(int L, int N){
		this->L = L;
		this->N = N;
		
		this->phi.resize(L,N);
		this->phi = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(L,N);
		
	}
	
	
	
	template <typename T> int Slater_determinant<T>::ref_L(){
		return L;
	}
	
	template <typename T> int Slater_determinant<T>::ref_N(){
		return N;
	}
	
	template <typename T> T Slater_determinant<T>::green_function(int i, int j){
		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vi = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(L, N+1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vj = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(L, N+1);
		
		
		vi(i, N) = 1.0;
		vj(j, N) = 1.0;
		
		for(int l = 0; l <= i; l++){
			vi.block(l, 0, 1, N) = phi.row(l);
		}
		
		for(int l = i + 1; l < L; l++){
			vi.block(l, 0, 1, N) = -phi.row(l);
		}
		
		for(int l = 0; l <= j; l++){
			vj.block(l, 0, 1, N) = phi.row(l);
		}
		
		for(int l = j + 1; l < L; l++){
			vj.block(l, 0, 1, N) = -phi.row(l);
		}
		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vv = vi.adjoint() * vj;
		
		return vv.determinant();
		
	}
	
	template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Slater_determinant<T>::matrix() const{
		return phi;
	}
	
	template <typename T> template <typename T2> double Slater_determinant<T>::expectation_value(Eigen::Matrix<T2, Eigen::Dynamic, Eigen::Dynamic> O){
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vv(N,N);
		double value = 0.0;
		
		//各行で最も大きな非ゼロ成分をもつ列
		Eigen::VectorXi last = Eigen::VectorXi::Zero(L);
		for(int i = 0; i < L; i++){
			for(int j = L-1; j > i; j--){
				if(abs(O(i,j)) > 1e-15){
					last(i) = j;
					break;
				}
			}
		}
		
		
		for(int i = 0; i < L; i++){
			//diagonal element
			if(abs(O(i,i)) > 1e-15){
				value += std::real(O(i,i) * T(phi.row(i) * phi.row(i).adjoint()));
			}
			
			//non-diagonal elements
			vv = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(N,N);
			for(int j = i+1; j <= last(i); j++){
				vv -= 2.0 * phi.row(j).adjoint() * phi.row(j);
				
				if(abs(O(i,j)) > 1e-15){
					Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp(N+1, N+1);
					tmp.block(0, 0, N, N) = vv;
					tmp.block(N, 0, 1, N) = phi.row(i);
					tmp.block(0, N, N, 1) = -phi.row(j).adjoint();
					tmp(N, N) = 0;
					value += std::real(std::conj(O(i,j)) * tmp.determinant()) * 2;
				}
			}
		}
		
		return value;
	}
	
	template <typename T> T Slater_determinant<T>::operator* (const Slater_determinant<double> &a) const{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp(N,N);
		tmp = phi.adjoint() * a.matrix();
		return tmp.determinant();
	}
	
	template <typename T> std::complex<double> Slater_determinant<T>::operator* (const Slater_determinant<std::complex<double> > &a) const{
		Eigen::MatrixXcd tmp(N,N);
		tmp = phi.adjoint() * a.matrix();
		return tmp.determinant();
	}
	
	template <typename T> Slater_determinant<T> Slater_determinant<T>::operator* (const Eigen::MatrixXd &M) const{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp = phi * M;
		return Slater_determinant<T>(tmp);
		
	}
	
	template <typename T> Slater_determinant<std::complex<double> > Slater_determinant<T>::operator* (const Eigen::MatrixXcd &M) const{
		Eigen::MatrixXcd tmp = phi * M;
		return Slater_determinant<std::complex<double> >(tmp);
	}
	
	template class Slater_determinant<double>;
	template class Slater_determinant<std::complex<double> >;
	template double Slater_determinant<double>::expectation_value<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>);
	template double Slater_determinant<double>::expectation_value<std::complex<double> >(Eigen::MatrixXcd);
	template double Slater_determinant<std::complex<double>>::expectation_value<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>);
	template double Slater_determinant<std::complex<double>>::expectation_value<std::complex<double> >(Eigen::MatrixXcd);
	
	
	free_fermion::free_fermion(Eigen::MatrixXd H, int N){
		L = H.cols();
		this->N = N;
		es.compute(H);
		
		basis.init(L, N);
		d = basis.ref_d();
		
	}
	
	Slater_determinant<double> free_fermion::eigenvector(long index){
		Eigen::MatrixXd P = Eigen::MatrixXd::Zero(L,N);
		
		Eigen::VectorXi bits = basis.bits(index);
		
		int n = 0;
		for(int i = 0; i < L; i++){
			if(bits(i) == 1){
				P.col(n) = es.eigenvectors().col(i);
				n++;
			}
		}
		
		return Slater_determinant<double>(P);
	}
	
	double free_fermion::eigenvalue(long index){
		Eigen::VectorXi bits = basis.bits(index);
		double E = 0;
		
		for(int i = 0; i < L; i++){
			if(bits(i) == 1){
				E += es.eigenvalues()(i);
			}
		}
		
		return E;
	}
	
	Slater_determinant<std::complex<double> > free_fermion::time_evolution(Slater_determinant<double> &phi, double t){
		Eigen::MatrixXcd e = Eigen::MatrixXcd::Zero(L,L);
		for(int i = 0; i < L; i++){
			double theta = es.eigenvalues()(i) * t;
			e(i,i) = std::complex<double>(cos(theta), sin(theta));
		}
		
		Eigen::MatrixXcd Uphi = es.eigenvectors() * e * es.eigenvectors().adjoint() * phi.matrix();
		
		return Slater_determinant<std::complex<double> >(Uphi);
	}
	
	Slater_determinant<std::complex<double> > free_fermion::time_evolution(Slater_determinant<std::complex<double> > &phi, double t){
		Eigen::MatrixXcd e = Eigen::MatrixXcd::Zero(L,L);
		for(int i = 0; i < L; i++){
			double theta = es.eigenvalues()(i) * t;
			e(i,i) = std::complex<double>(cos(theta), sin(theta));
		}
		
		Eigen::MatrixXcd Uphi = es.eigenvectors().adjoint() * e * es.eigenvectors() * phi.matrix();
		
		return Slater_determinant<std::complex<double> >(Uphi);
	}
	
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> free_fermion::ref_es(){
		return es;
	}
	
	long free_fermion::ref_d(){
		return d;
	}
	
	int free_fermion::ref_L(){
		return L;
	}
	
	int free_fermion::ref_N(){
		return N;
	}
}
