#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <time.h>

namespace tylib{
	
	class loger{
	  public:
		loger(int argc, char* argv[], std::string _comment = "", std::string _logname = "log.txt" );
		void check_point();
		void show();
		int ref_time(int i, int j=0);
		void append_info(std::string s);
		void output_file(std::string output_name);
		void input_file(std::string input_name);
		~loger();
		
	  private:
		std::string logname;
		std::string comment;
		int npoints;
		std::vector<clock_t> t;
		std::stringstream begin_time;
		std::stringstream begin_date;
		std::stringstream fileinfo;
		
		std::stringstream appendix;
		
	};
	
	loger::loger(int argc, char* argv[], std::string _comment, std::string _logname){
		logname = _logname;
		comment = _comment;
		t.push_back(clock());
		npoints = 1;
		
		time_t now = time(NULL);
		struct tm *pnow = localtime(&now);
		std::stringstream date, time;
		date << "開始日時 : "<< pnow->tm_year+1900 << "年" << pnow->tm_mon + 1 << "月" << pnow->tm_mday << "日  " << pnow->tm_hour << "時" << pnow->tm_min << "分" << pnow->tm_sec << "秒" << std::endl;
		time << "開始時刻 : " << pnow->tm_hour << "時" << pnow->tm_min << "分" << pnow->tm_sec << "秒" << std::endl;
		begin_time << time.str();
		begin_date << date.str();
		
		fileinfo << "実行ファイル : ";	//実行ファイル名
		for(int i = 0; i < argc; i++){
			fileinfo << argv[i] << " ";
		}
		
		show();
		
	}
	
	void loger::check_point(){
		t.push_back(clock());
		npoints++;
	}
	
	void loger::show(){
		std::cout << begin_date.str();
		//std::cout << begin_time.str();
		std::cout << fileinfo.str() << std::endl;
		
		if( comment.size() != 0 ){
			std::stringstream parameter;
			parameter << "パラメータ : " << comment << std::endl;
			std::cout << parameter.str();
		}
		
		if( appendix.str().size() != 0 ){
			std::cout << std::endl << appendix.str() << std::endl;
		}
	}
	
	int loger::ref_time(int i, int j){
		double diff = double(t[i] - t[j]);
		return int(diff/1000000);
		
		
	}
	
	void loger::append_info(std::string s){
		appendix << s;
	}
	
	void loger::output_file(std::string output_name){
		appendix << "出力ファイル : " << output_name << std::endl;
	}
	
	void loger::input_file(std::string input_name){
		appendix << "入力ファイル : " << input_name << std::endl;
	}
	
	loger::~loger(){
		std::ofstream log(logname.c_str(), std::ios::app);
		
		std::stringstream partition;
		partition << std::endl << "------------------------------------------------------------------------" << std::endl;
		log << partition.str();
		
		//日時の出力
		time_t now = time(NULL);
		struct tm *pnow = localtime(&now);
		std::stringstream date, time;
		date << "終了日時 : "<< pnow->tm_year+1900 << "年" << pnow->tm_mon + 1 << "月" << pnow->tm_mday << "日  " << pnow->tm_hour << "時" << pnow->tm_min << "分" << pnow->tm_sec << "秒" << std::endl;
		time << "終了時刻 : " << pnow->tm_hour << "時" << pnow->tm_min << "分" << pnow->tm_sec << "秒" << std::endl;
		log << begin_date.str();
		//log << begin_time.str();
		log << date.str();
		//log << time.str();
		log << fileinfo.str() << std::endl;
		
		if( comment.size() != 0 ){
			std::stringstream parameter;
			parameter << "パラメータ : " << comment << std::endl;
			log << parameter.str();
		}
		
		if( appendix.str().size() != 0 ){
			log << std::endl << appendix.str() << std::endl;
		}
		
		check_point();
		
		log << "実行時間 : " << ref_time(npoints-1) << " 秒" << std::endl << std::endl;
		
		
	}
}
