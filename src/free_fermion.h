#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>

#ifndef INCLUDE_BIT
#define INCLUDE_BIT
#include "bit.h"
#endif

namespace tylib{
	
	
	
	void XX(Eigen::MatrixXd &H, int N, double J = 1.0, double h = 0, bool isPeriodic = true);
	
	template <typename T> class Slater_determinant{
	  public:
		Slater_determinant(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &phi);
		Slater_determinant(computational_basis basis, int index);
		Slater_determinant(int L, int N);
		int ref_L();
		int ref_N();
		T green_function(int i, int j);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix() const;
		template <typename T2> double expectation_value(Eigen::Matrix<T2, Eigen::Dynamic, Eigen::Dynamic>  O);
		//double expectation_value(Eigen::MatrixXcd O);
		
		T operator* (const Slater_determinant<double> &a) const;
		std::complex<double> operator* (const Slater_determinant<std::complex<double> > &a) const;
		Slater_determinant<T> operator* (const Eigen::MatrixXd &M) const;
		Slater_determinant<std::complex<double> > operator*(const Eigen::MatrixXcd &M) const;
		T & operator() (int i, int j) {return phi(i,j);}
		
		
	  private:
		int L;
		int N;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> phi;
		
	};
	
	
	
	std::ostream & operator<< (std::ostream &os, const Slater_determinant<double> &a){
		os << a.matrix();
		return os;
	}
	
	std::ostream & operator<< (std::ostream &os, const Slater_determinant<std::complex<double> > &a){
		os << a.matrix();
		return os;
	}
	
	class free_fermion{
		
	  public:
		free_fermion(Eigen::MatrixXd H, int N);
		Slater_determinant<std::complex<double> > time_evolution(Slater_determinant<double> &phi, double t);
		Slater_determinant<std::complex<double> > time_evolution(Slater_determinant<std::complex<double> > &phi, double t);
		Slater_determinant<double> eigenvector(long index);
		double eigenvalue(long index);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref_es();
		long ref_d();
		int ref_L();
		int ref_N();
		
	  private:
		computational_basis basis;
		int L;
		int N;
		long d;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
		
	};
	
	
}