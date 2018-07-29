#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <complex>
#include <string>
#include <random>
#include <omp.h>

extern "C"{
	extern void komega_bicg_init(int* ndim,int* nl,int* nz, std::complex<double>* x, std::complex<double>* z, int* itermax, double* threshold);
	extern void komega_bicg_update(std::complex<double>* v12, std::complex<double>* v2, std::complex<double>* v14, std::complex<double>* v4, std::complex<double>* x, std::complex<double>* r_l, int* status);
	extern void komega_bicg_finalize();
	extern void komega_cocg_init(int* ndim,int* nl,int* nz, std::complex<double>* x, std::complex<double>* z, int* itermax, double* threshold);
	extern void komega_cocg_update(std::complex<double>* v12, std::complex<double>* v2, std::complex<double>* x, std::complex<double>* r_l, int* status);
	extern void komega_cocg_finalize();
};


namespace tylib{
	
	class komega_bicg{
	  public:
		komega_bicg() {};
		komega_bicg(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold);
		~komega_bicg();
		void solve(Eigen::SparseMatrix<std::complex<double> > &H);
		void init(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold);
		void update(Eigen::VectorXcd& Hr1, Eigen::VectorXcd& Hr2);
		Eigen::VectorXcd getSolution(int i);
		int iterations();
		bool isConverged();
		int status();
		int seed();
		double residue();
		
	  private:
		int _status[3];
		int itermax;
		double threshold;
		int nz, d;
		
		std::complex<double>* v12;
		std::complex<double>* v2;
		std::complex<double>* v14;
		std::complex<double>* v4;
		std::complex<double>* x;
		std::complex<double>* z;
		
		
	};
	
	class komega_cocg{
	  public:
		komega_cocg() {};
		komega_cocg(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold);
		~komega_cocg();
		void solve(Eigen::SparseMatrix<double> &H);
		void init(Eigen::VectorXcd& v, Eigen::VectorXcd& z, int itermax, double threshold);
		void update(Eigen::VectorXcd& Hr);
		Eigen::VectorXcd getSolution(int i);
		int iterations();
		bool isConverged();
		int status();
		int seed();
		double residue();
		
	  private:
		int _status[3];
		int itermax;
		double threshold;
		int nz, d;
		
		std::complex<double>* v12;
		std::complex<double>* v2;
		std::complex<double>* x;
		std::complex<double>* z;
		
		
	};
	
	class random_vector{
	  public:
		random_vector();
		random_vector(double seed);
		Eigen::VectorXcd gauss(int d);		// uniform distribution
		Eigen::VectorXcd binary(int d);		// 1 or -1
		
		
	  private:
		std::mt19937 mt;
		
	};
	
	class ss_method{
	  public:
		
		ss_method(Eigen::SparseMatrix<double> &H);
		ss_method(Eigen::SparseMatrix<std::complex<double> > &H);
		
		//固有ベクトルの計算
		//Eigen::MatrixXd eigenvector(double gamma, double rho);
		void eigenstate(Eigen::MatrixXcd &evec, Eigen::VectorXd &E, double gamma, double rho);
		
		//射影ベクトルの計算
		Eigen::MatrixXcd projection(Eigen::VectorXcd &v, double gamma, double rho, int m);
		
		//経路内の固有状態数の見積もり
		int num_states(double gamma, double rho);
		
		//固有ベクトルの精度の確認
		void precision_check(Eigen::MatrixXcd &evec, Eigen::VectorXd &E);
		
		//射影ベクトルから固有ベクトルを計算
		void projection_to_eigenstate(Eigen::MatrixXcd &evec, Eigen::VectorXd &E, Eigen::MatrixXcd &S);
		
		//パラメータの設定
		void enable_openmp(int n_threads);
		void memory_setting(int memory);
		void itermax_setting(int itermax);
		void threshold_setting(int threshold);
		void contour_setting(int contour, int n, double contour_para);
		void copy_setting(int m, int step);
		void kappa_setting(double kappa);
		void cutoff_setting(double cutoff);
		
	  private:
		bool isComplex;		//行列が実か複素か
		Eigen::SparseMatrix<double> Hd;
		Eigen::SparseMatrix<std::complex<double> > Hcd;
		int d;
		
		//計算上のパラメータ
		int n_threads;	//openmpの並列数
		int memory_lim;		//使用できるメモリの上限
		int itermax;		//shifted bicg method の反復回数の上限
		double threshold;	//shifted bicg method の残さの上限
		int contour;		//積分経路の形状 "circle" "ellipse" "rectangle"
		int n;				//複素積分の計算に用いる点数
		double contour_para;		//積分経路の特徴づけパラメータ
		int m;		//射影ベクトルの複製を何本まで行うか
		int step;		//射影ベクトルの複製を何本間隔で行うか |s^0>,|s^{m_step}>,|s^{2m_step}>,...,|s^{m_lin*m_step}>
		double kappa;	//積分経路内の固有状態数の見積もりの何倍の本数の射影ベクトルを計算するか
		double cutoff;		//特異値分解で最大特異値との比がcutoffを下回ると特異値0とみなす
		
		//メソッド
		std::complex<double> prefactor(double gamma, double rho, int j);
		Eigen::VectorXcd make_contour(double gamma, double rho);
		
	};
	
}