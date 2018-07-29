#include <tylib/include/log.h>
#include <tylib/include/observable.h>		//observable_base.h はinclude しなくてよい
#include <tylib/include/translation.h>
#include <tylib/include/ss.h>

using namespace std;
using namespace Eigen;
using namespace tylib;
typedef complex<double> cd;

int main(int argc, char* argv[]){
	
	const int L = atoi(argv[1]);				//サイト数
	const int N = atoi(argv[2]);				//粒子数
	
	//計算基底を張って，次元を計算
	computational_basis basis(L, N);
	long d = basis.ref_d();
	
	
	//実行情報とコメントをログファイルに出力
	stringstream comment;
	comment << "L=" << L << ", N=" << N << ", d=" << d << endl;
	loger Log(argc, argv, comment.str());
	
	
	//ハミルトニアンを定義
	XX_model hami(1.0, 1.0, lambda);
	
	
	//並進対称な基底を形成
	translational_basis t_basis(basis);
	
	
	//並進操作の固有値ごとに固有ベクトルを計算
	for(int r = 0; r < L; r++){
		//次元を計算
		int D = t_basis.ref_D(r);
		
		//ハミルトニアンの並進基底のもとでの行列成分を計算
		SparseMatrix<cd> H(D,D);
		t_basis.cal_observable(H, hami, r);
		
		MatrixXcd evec;
		VectorXd E;
		
		//SS法で[-0.1,0.1]の固有ベクトルを計算
		ss_method SS(H);
		SS.eigenstate(evec, E, 0, 0.1);
		
		//固有ベクトルの精度を確認
		SS.precision_check(evec, E);
		
	}
	
	
	
	return 0;
}