#ifndef INCLUDE_BIT
#define INCLUDE_BIT
#include "bit.h"
#endif
#ifndef INCLUDE_OBSERVABLE_BASE
#define INCLUDE_OBSERVABLE_BASE
#include "observable_base.h"
#endif
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace tylib{
	
	
	void translation(Eigen::VectorXi& bits);
	
	class translational_basis{
		
	  public:
		translational_basis(computational_basis &c_basis);
		int ref_d();
		int ref_D(int r);
		
		//並進対称でない物理量の計算を行う場合には，3つ目の引数にfalseを入れる
		void cal_observable(Eigen::SparseMatrix<std::complex<double> > &obs, observable_base &operation, int r, bool isTranslational = true);
		void cal_observable(Eigen::SparseMatrix<std::complex<double> > &obs, observable_base_complex &operation, int r, bool isTranslational = true);
		
		
		int ref_D_period(int i);
		std::vector<int> ref_period();
		std::vector<int> ref_rperiod(int r);
		//Eigen::VectorXcd cal_c_basis(Eigen::VectorXcd v, int r);
		Eigen::VectorXcd eigenvector(int r, int index);			//並進操作の固有ベクトルを計算基底で取得
		Eigen::VectorXcd c_to_t(Eigen::VectorXd &vc, int r);
		Eigen::VectorXcd c_to_t(Eigen::VectorXcd &vc, int r);
		Eigen::VectorXcd t_to_c(Eigen::VectorXcd &vt, int r);
		std::pair<int, std::complex<double> > c_basis_to_t_basis(int c_index, int r);
		
	  private:
		int L;
		int N;
		int d;
		
		int cal_t_index(int i, int r);		//i番目の計算基底は固有値rの何番目の固有ベクトルに属するか
		void cal_basis(std::vector<int> &indexes, int r, int index);
		
		Eigen::VectorXi D;
		std::vector<std::vector<int> > r_period;		//各固有値で固有ベクトルが存在する周期
		Eigen::VectorXi D_period;		//各周期の固有ベクトル数
		std::vector<int> period;
		int p;						//ありえる周期の数
		
		std::vector<std::vector<std::vector<int> > > basis;			//固有ベクトルのインデックス[周期][グループのインデックス][グループ内でのインデックス] 
		std::vector<std::vector<int> > table;
		std::vector<std::vector<int> > t_index_init;		//各固有値で周期の始めのベクトルが何番目か
		//std::vector<c_t_table> table;				//計算基底が属する，周期，周期中でのグループのインデックス，グループ内でのインデックス
	};
	
}
