#include <cv.h>      
#include <ml.h>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <chrono>       // std::chrono::system_clock
#include <random>       // std::default_random_engine

using std::vector;
using std::string;
using namespace std;
vector<int> good_randVec(int a, int b,int num);
vector<int> randVec(int a,int b, int num);

int main(int argc, char *argv[])
{
	 if( argc < 2 )
	 {
		 std::cerr << "Usage: " << std::endl;
		 std::cerr << argv[0] << "offsetSettings  windowSettings " << std::endl;
		 return EXIT_FAILURE;
	 }
	 ofstream outfile1;
	 outfile1.open(argv[1], std::ofstream::out);  
	 ofstream outfile2;
	 outfile2.open(argv[2], std::ofstream::out);  


	 // Generate a 500 points array chosen from the mask voxels randomly
	 const int num1=1000;
	 vector<int> offsetX;
	 vector<int> offsetY;
	 vector<int> offsetZ;
	 /*offsetX=good_randVec(-80,80,num1);
	 offsetY=good_randVec(-80,80,num1);
	 offsetZ=good_randVec(-40,40,num1);*/
	 offsetX=randVec(-80,80,num1);
	 offsetY=randVec(-80,80,num1);
	 offsetZ=randVec(-15,15,num1);
	 outfile1<<num1<<endl;
	// cout<<"Good"<<endl;
 	 for (int i=0;i<num1;i++)
		 outfile1<<offsetX[i]<<" "<<offsetY[i]<<" "<<offsetZ[i]<<endl;
	 vector<int> WindowX;
	 vector<int> WindowY;
	 vector<int> WindowZ;
	 /*WindowX=good_randVec(1,20,num1);
	 WindowY=good_randVec(1,20,num1);
	 WindowZ=good_randVec(1,10,num1);*/
	 WindowX=randVec(1,20,num1);
	 WindowY=randVec(1,20,num1);
	 WindowZ=randVec(1,3,num1);
	 outfile2<<num1<<endl;
 	 for (int i=0;i<num1;i++)
		 outfile2<<WindowX[i]<<" "<<WindowY[i]<<" "<<WindowZ[i]<<endl;

}

vector<int> good_randVec(int a, int b,int num)
{
	static default_random_engine e;
	static uniform_int_distribution<int> u(a,b);
	vector<int> ret;
	for (int i=0;i<num;i++)
		ret.push_back(u(e));
	return ret;
}
vector<int> randVec(int a,int b, int num)
{
	vector<int> ret(num);
	for (int i =0;i<num;i++)
		ret[i]=rand()%(b-a+1)+a;
	return ret;
}