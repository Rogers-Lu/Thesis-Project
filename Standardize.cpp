// This program is used to standardize the features;
//input should be the feature file, and it will output the standardized feature file;

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h" 
#include "itkImageRegionIterator.h" //write the results of computations to the output image
#include <fstream>
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <string>
#include "itkImageDuplicator.h"   // to save another image which include the extracted features
#include "itkImageFileWriter.h"
#include <vector>
#include <cmath>
using std::vector;
using std::string;
using namespace std;

//int const num_offset=30;
typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::Image< float, 3 >                ScanImageType;
#define NUMBER_OF_TRAINING_SAMPLES 30000
#define OFFSET_PER_SAMPLE 120//include all the number of features 
#define NUMBER_OF_TESTING_SAMPLES 10000
#define NUMBER_OF_CLASSES 2
void featureStandardize(vector<float> v);
int read_data_from_csv(const char* filename, 	vector<vector<float>> vv1, vector<float> vv1Class,
                       int n_samples );
int main(int argc, char *argv[])
{
	if( argc < 4 )
	{
		std::cerr << "Usage: " << std::endl;

		//featureSettings is a file containing the number of offset, offset matrix, the number of window and the size of the window 
		std::cerr << argv[0] << "trainingFeature testFeature trainingFetureStandardized testFeatureStandardized" << std::endl;
		return EXIT_FAILURE;
  }
	ifstream infile1;
    infile1.open (argv[1]);
    ifstream infile2;
    infile2.open (argv[2]);
	ofstream outfile1;
	outfile1.open(argv[3], std::ofstream::out);  
	ofstream outfile2;
	outfile2.open(argv[4], std::ofstream::out);  

    //Define 2 two dimensional vectors used to standardize the generated feature file
	vector<vector<float> >vv1(NUMBER_OF_TRAINING_SAMPLES, vector<float>(OFFSET_PER_SAMPLE));
    vector<float >vv1Class(NUMBER_OF_TRAINING_SAMPLES);
	vector<vector<float> >vv2(NUMBER_OF_TESTING_SAMPLES, vector<float>(OFFSET_PER_SAMPLE));
    vector<float >vv2Class(NUMBER_OF_TESTING_SAMPLES);
	vector<vector<float> >vv(NUMBER_OF_TRAINING_SAMPLES+NUMBER_OF_TESTING_SAMPLES, vector<float>(OFFSET_PER_SAMPLE*4+1));



	//***********************************************************************************//
	//Populate the vectors
	//Be aware that the feature files are separated by comma
	
	if(read_data_from_csv(argv[1], vv1, vv1Class, NUMBER_OF_TRAINING_SAMPLES) &&
            read_data_from_csv(argv[2], vv2, vv2Class, NUMBER_OF_TESTING_SAMPLES))
	{

		cout<<vv1[0][0];
		//Merge training data vector and test data vector into one vector
		for(int i = 0;i<NUMBER_OF_TRAINING_SAMPLES;i++)
		{
			for (int j=0;j<(OFFSET_PER_SAMPLE*4+1);j++)
			{
				if(j!=OFFSET_PER_SAMPLE*4)
				{
					vv[i][j]=vv1[i][j];
					//cout<<vv[i][j];
				}					
				else
				{
					vv[i][j]=vv1Class[i];
				}
			}
		}
		for(int i = 0;i<NUMBER_OF_TESTING_SAMPLES;i++)
		{
			for (int j=0;j<(OFFSET_PER_SAMPLE*4+1);j++)
			{
				if(j!=OFFSET_PER_SAMPLE*4)
				{
					vv[i+NUMBER_OF_TRAINING_SAMPLES][j]=vv2[i][j];
				}					
				else
				{
					vv[i+NUMBER_OF_TRAINING_SAMPLES][j]=vv2Class[i];
				}
			}
		}
     //


	 // Standardize the features we got
	 vector<float> vTemp(NUMBER_OF_TRAINING_SAMPLES+NUMBER_OF_TESTING_SAMPLES);
	 for(int a=0;a<OFFSET_PER_SAMPLE*4;a++)
	 {
		 for(int b=0;b<NUMBER_OF_TRAINING_SAMPLES+NUMBER_OF_TESTING_SAMPLES;b++)
		 {
			 vTemp[b]=vv[b][a];
		 }	
		 featureStandardize(vTemp);
		 for(int c=0;c<OFFSET_PER_SAMPLE*4;c++)
		 {
			 vv[c][a]=vTemp[c];
			 //cout<<vv1[c][a];
		 }
	 }

	 //Write the standardized features to the training file and test file
	 for(int i=0;i<NUMBER_OF_TRAINING_SAMPLES;i++)
	 {
		 for(int j=0;j<(OFFSET_PER_SAMPLE*4+1);j++)
		 {
			 outfile1<<vv[i][j];

		 }
	 }

	 	 for(int i=0;i<NUMBER_OF_TESTING_SAMPLES;i++)
	 {
		 for(int j=0;j<(OFFSET_PER_SAMPLE*4+1);j++)
		 {
			 outfile2<<vv[i+NUMBER_OF_TRAINING_SAMPLES][j];

		 }
	 }
	}


	 infile1.close();
     infile2.close();
     outfile1.close();
     outfile2.close();
     return EXIT_SUCCESS;

}


void featureStandardize(vector<float> v)
{
	int length=v.size();
	float mean=0.0, dev=0.0,sum1=0.0,sum2=0.0;
	for(int i=0;i<length;i++)
	{
		sum1=sum1+v[i];
	}
	mean=sum1/length;
	for(int j=0;j<length;j++)
	{
		sum2 = sum2+(v[j]-mean)*(v[j]-mean);
	}
	dev=sqrt(sum2);
	//standardize
	for(int k=0;k<length;k++)
	{
		v[k]=(v[k]-mean)/dev;
	}
}

int read_data_from_csv(const char* filename, 	vector<vector<float>> vv1, vector<float> vv1Class,
                       int n_samples )
{
    float tmp;

    // if we can't read the input file then return 0
    FILE* f = fopen( filename, "r" );
    if( !f )
    {
        printf("ERROR: cannot read file %s\n",  filename);
        return 0; // all not OK
    }

    // for each sample in the file

    for(int line = 0; line < n_samples; line++)
    {

        // for each attribute on the line in the file

        for(int attribute = 0; attribute < (OFFSET_PER_SAMPLE + 1); attribute++)
        {
            if (attribute < OFFSET_PER_SAMPLE)
            {

                // The last 10 elements  in each line are the attributes

                fscanf(f, "%f,", &tmp);
                vv1[line][attribute] = tmp;
                // printf("%f,", data.at<float>(line, attribute));

            }
            else if (attribute == OFFSET_PER_SAMPLE)
            {

                // The last attibute is the class label {0,1}

                fscanf(f, "%f,", &tmp);
				vv1Class[line] = tmp;
                // printf("%f\n", classes.at<float>(line, 0));

            }
        }
    }

    fclose(f);

    return 1; // all OK
}

