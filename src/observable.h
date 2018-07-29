#ifndef INCLUDE_OBSERVABLE_BASE
#define INCLUDE_OBSERVABLE_BASE
#include "observable_base.h"
#endif
#include <random>

namespace tylib{
	
	class ising_interaction : public observable_base{
	  public:
		ising_interaction() {J = 1; distance = 1; isNormalized = false;}
		ising_interaction(double J, int distance = 1, bool isNormalized = false);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void parameters(double J, int distance = 1, bool isNormalized = false);
		
	  private:
		double J;
		int distance;
		bool isNormalized;
		
	};
	
	class z_magnetic_field : public observable_base{
	  public:
		z_magnetic_field() {h = 1; isNormalized = false;}
		z_magnetic_field(double h, bool isNormalized = false);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void parameters(double h, bool isNormalized = false);
		
	  private:
		double h;
		bool isNormalized;
		
	};
	
	class x_ising_interaction : public observable_base{
	  public:
		x_ising_interaction() {J = 1; distance = 1; isNormalized = false;}
		x_ising_interaction(double J, int distance = 1, bool isNormalized = false);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void parameters(double J, int distance = 1, bool isNormalized = false);
		
	  private:
		double J;
		int distance;
		bool isNormalized;
		
	};
	
	class x_magnetic_field : public observable_base{
	  public:
		x_magnetic_field() {h = 1; isNormalized = false;}
		x_magnetic_field(double h, bool isNormalized = false);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void parameters(double h, bool isNormalized = false);
		
	  private:
		double h;
		bool isNormalized;
	};
	
	class XX_model : public observable_base{
	  public:
		XX_model() {J = 1; distance = 1;}
		XX_model(double J, int distance = 1);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void parameters(double J, int distance = 1);
		
	  private:
		double J;
		int distance;
	};
	
	class HCB_interaction : public observable_base{
	  public:
		HCB_interaction(){}
		HCB_interaction(double U, int distance = 1);
		void parameters(double U, int distance = 1);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		
	  private:
		double U;
		int distance;
		
	};
	
	class spin_XXZ_model : public observable_sum{
	  public:
		spin_XXZ_model() {}
		spin_XXZ_model(double J, double delta, int distance = 1);
		void parameters(double J, double delta, int distance = 1);
		
	  private:
		XX_model XX;
		ising_interaction Z;
		
	};
	
	class XXZ_model : public observable_sum{
	  public:
		XXZ_model() {}
		XXZ_model(double J, double delta, int distance = 1);
		void parameters(double J, double delta, int distance = 1);
		
	  private:
		XX_model XX;
		HCB_interaction Z;
		
	};
	
	class XXZ_with_NNN : public observable_sum{
	  public:
		XXZ_with_NNN(){}
		XXZ_with_NNN(double J, double delta, double lambda);
		void parameters(double J, double delta, double lambda);
		
	  private:
		XXZ_model XXZ1, XXZ2;
		
	};
	
	class Hubbard_model : public observable_base{
	  public:
		Hubbard_model(double J, double U);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void parameters(double J, double U);
		
	  private:
		double J;
		double U;
		
	};
	
	class MBL : public spin_XXZ_model{
	  public:
		MBL(double W, double J = 1.0);
		MBL(double W, double J, int seed);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		void reset_disorder();
		
	  private:
		double W;
		std::mt19937 mt;
		std::vector<double> h;
		
	};
	
	class transverse_field_ising : public observable_sum{
	  public:
		transverse_field_ising() {}
		transverse_field_ising(double J, double h);
		void parameters(double J, double h);
		
	  private:
		ising_interaction SzSz;
		x_magnetic_field Sx;
		
	};
	
	class n1n2 : public observable_base{
	  public:
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		
	  
	};
	
	class hopping : public observable_base{
	  public:
		hopping(int distance);
		std::vector<std::pair<int, double> > cal_elements(int index, computational_basis &basis);
		
	  private:
		int distance;
		
	};
	
	class momentum : public observable_base_complex{
	  public:
		momentum(int k);
		std::vector<std::pair<int, std::complex<double> > > cal_elements(int index, computational_basis &basis);
		
	  private:
		int k;		//wavenumber
		
	};
	
	class momentum_0 : public observable_base{
	  public:
		std::vector<std::pair<int,double> > cal_elements(int index, computational_basis &basis);
		std::pair<bool, double> coeff(int i, int j, computational_basis &basis);
	};
	
}