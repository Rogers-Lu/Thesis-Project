// This program is used to update the mask.
//The inputs are classification results, the corresponding index and the the orignial mask;
//The output is the updated mask image.

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

using std::string;
using namespace std;

typedef itk::Image< unsigned short, 3 >       SegmImageType;

int main(int argc, char *argv[])
{
  /*if( argc < 4 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputofImage inputClassificationFile inputIndexFile outputofImage" << std::endl;
    return EXIT_FAILURE;
  }*/

  typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;
  /*
    // Read input segmentation
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(argv[1]);
  segmReader->Update();
  std::cout << "Loaded image " << argv[1] << std::endl;
  */
 /* ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);*/

  //***************************************************************************************//
  //Count the number of line in features txt file

  int number_of_lines = 0;
  string line;
  ifstream infile;
  infile.open (argv[1]);

  while (std::getline(infile, line))
	{
        ++number_of_lines;
	}

    cout << "Number of lines in text file: " << number_of_lines;

  //***************************************************************************************//
  //infile1.close();
  //infile2.close();
    return EXIT_SUCCESS;

}