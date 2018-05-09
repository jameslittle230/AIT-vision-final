// Source file for image class



// Include files 

#include <vector>
#include <cmath>
#include <iostream>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"




////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}

PointCorrespondence createCorrespondence(Feature ft) {
    PointCorrespondence output = {{ft.centerX, ft.centerY}, {ft.x2, ft.y2}};
    return output;
};

/* ================================================
   ============ IMAGE PROCESSING CODE ============
   ================================================
*/

std::vector<Feature> R2Image::
FirstFrameProcessing() {
  // @TODO

  return this->Harris(2.0);

}

void R2Image::
FrameProcessing(R2Image * mainImage,std::vector<Feature> features) {

  // @TODO
  this->blendImages(mainImage,features);
  //this->blendImages(this,features);
}

void R2Image::
SobelX(void)
{
	R2Image *output = new R2Image(width, height);

  double sobel_x[3][3] = {
    {-1, 0, 1},
    {-2, 0, 2},
    {-1, 0, 1}
  };

  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      double sobelVal =
        (sobel_x[0][0] * Pixel(i-1,j-1).Luminance()) + (sobel_x[0][1] * Pixel(i,j-1).Luminance()) + (sobel_x[0][2] * Pixel(i+1,j-1).Luminance()) +
        (sobel_x[1][0] * Pixel(i-1,j).Luminance())   + (sobel_x[1][1] * Pixel(i,j).Luminance())   + (sobel_x[1][2] * Pixel(i+1,j).Luminance()) +
        (sobel_x[2][0] * Pixel(i-1,j+1).Luminance()) + (sobel_x[2][1] * Pixel(i,j+1).Luminance()) + (sobel_x[2][2] * Pixel(i+1,j+1).Luminance());

      R2Pixel *newPixel = new R2Pixel(sobelVal, sobelVal, sobelVal, 1);
      output->SetPixel(i, j, *newPixel);
    }
  }

  this->pixels = output->pixels;
  output->pixels = nullptr;
  delete output;
}

void R2Image::
SobelY(void)
{
	R2Image *output = new R2Image(width, height);

  double sobel_y[3][3] = {
    {-1, -2, -1},
    {0, 0, 0},
    {1, 2, 1}
  };

  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      double sobelVal =
        (sobel_y[0][0] * Pixel(i-1,j-1).Luminance()) + (sobel_y[0][1] * Pixel(i,j-1).Luminance()) + (sobel_y[0][2] * Pixel(i+1,j-1).Luminance()) +
        (sobel_y[1][0] * Pixel(i-1,j).Luminance())   + (sobel_y[1][1] * Pixel(i,j).Luminance())   + (sobel_y[1][2] * Pixel(i+1,j).Luminance()) +
        (sobel_y[2][0] * Pixel(i-1,j+1).Luminance()) + (sobel_y[2][1] * Pixel(i,j+1).Luminance()) + (sobel_y[2][2] * Pixel(i+1,j+1).Luminance());

      R2Pixel *newPixel = new R2Pixel(sobelVal, sobelVal, sobelVal, 1);
      output->SetPixel(i, j, *newPixel);
    }
  }

  this->pixels = output->pixels;
  output->pixels = nullptr;
  delete output;
}

void R2Image::
Blur(double sigma)
{
  R2Image *temp = new R2Image(width, height);

  // Generate 1D kernel
  int k_rad = 3 * sigma;
  int k_len = 6 * sigma + 1;
  double *k = (double *) malloc(k_len * sizeof(double));
  double sum = 0.0;

  double norm = 1.0 / 2.50662 * sigma;
  double coeff = 2.0 * sigma * sigma;

  for(int i=-1 * k_rad; i<=k_rad; i++) {
    double val = norm * std::exp(-i * i / coeff);
    k[i+k_rad] = val;
    sum += k[i+k_rad];
  }

  for(int i=0; i<k_len; i++) {
    k[i] /= sum;
  }

  int x, y, i;

  for(y=0; y<height; y++) {
    for(x=0; x<width; x++) {
      double r = 0.0, g = 0.0, b = 0.0;
      for(i=-1*k_rad; i<=k_rad; i++) {
        r += Pixel(x, y+i).Red() * k[i+k_rad];
        g += Pixel(x, y+i).Green() * k[i+k_rad];
        b += Pixel(x, y+i).Blue() * k[i+k_rad];
      }
      // std::cout << val << std::endl;
      R2Pixel* p = new R2Pixel(r, g, b, 1);
      temp->SetPixel(x, y, *p);
    }
  }

  this->pixels = temp->pixels;
  temp->pixels = (R2Pixel *) malloc(npixels * sizeof(R2Pixel));

  for(y=k_rad; y<height-k_rad; y++) {
    for(x=k_rad; x<width-k_rad; x++) {
      double r = 0.0, g = 0.0, b = 0.0;
      for(i=-1*k_rad; i<=k_rad; i++) {
        r += Pixel(x+i, y).Red() * k[i+k_rad];
        g += Pixel(x+i, y).Green() * k[i+k_rad];
        b += Pixel(x+i, y).Blue() * k[i+k_rad];
      }
      // std::cout << val << std::endl;
      R2Pixel* p = new R2Pixel(r, g, b, 1);
      temp->SetPixel(x, y, *p);
    }
  }

  this->pixels = temp->pixels;
  temp->pixels = nullptr;
  free(k);
  delete temp;
}

void R2Image::
Square()
{
  int x, y;
  double r, g, b;
  for(x=0; x<width; x++) {
    for(y=0; y<height; y++) {
      r = Pixel(x, y).Red() * Pixel(x, y).Red();
      g = Pixel(x, y).Green() * Pixel(x, y).Green();
      b = Pixel(x, y).Blue() * Pixel(x, y).Blue();

      Pixel(x, y).Reset(r, g, b, 1);
    }
  }
}

/**
 * Draws a line from (x1, y1) to (x2, y2) with the RGB color specified by
 * R, G, and B (which are integers from 0 to 255) and draws a box around the
 * point (x2, y2)
 */
void R2Image::
drawLineWithBox(int x1, int y1, int x2, int y2, int r, int g, int b) {
  float rf = float(r)/255.0;
  float gf = float(g)/255.0;
  float bf = float(b)/255.0;

  int m, n, x;

  for(m=-4; m<=4; m++) {
      for(n=-4; n<=4; n++) {
        this->Pixel(x2+m, y2+n).Reset(rf, gf, bf, 1);
      }
    }

    int dx = x2 - x1;
    int dy = y2 - y1;

    if(dx != 0) { // avoid div by zero errors
      for(x=std::min(x1, x2); x<=std::max(x1, x2); x++) {
        int y = int(std::round(y1 + (double(dy * (x-x1)) / double(dx))));
        this->Pixel(x, y).Reset(rf, gf, bf, 1);
      }
    }
}

std::vector<Feature> R2Image::
Harris(double sigma)
{
  const R2Image self = *this;
  R2Image *t1 = new R2Image(self); t1->SobelX(); t1->Square();
  R2Image *t2 = new R2Image(self); t2->SobelY(); t2->Square();
  R2Image *t3 = new R2Image(Width(), Height());
  R2Image *t4 = new R2Image(Width(), Height());
  // Set T3 to product of T1 and T2
  for(int x=0; x<Width(); x++) {
    for(int y=0; y<Height(); y++) {
      double v = t1->Pixel(x, y)[0] * t2->Pixel(x, y).Red();
      t3->Pixel(x, y).Reset(v, v, v, 1);
    }
  }

  t1->Blur(sigma);
  t2->Blur(sigma);
  t3->Blur(sigma);

  for(int x=0; x<Width(); x++) {
    for(int y=0; y<Height(); y++) {
      double t1v = t1->Pixel(x, y)[0];
      double t2v = t2->Pixel(x, y)[0];
      double t3v = t3->Pixel(x, y)[0];
      double v = t1v * t2v - t3v * t3v - 0.04 * ((t1v + t2v) * (t1v + t2v));
      v += 0.5;
      t4->Pixel(x, y).Reset(v, v, v, 1);
    }
  }

  std::vector<Feature> features;
  std::vector<Feature> featuresOut;

  for(int x=0; x<Width(); x++) {
    for(int y=0; y<Height(); y++) {
      R2Pixel p = t4->Pixel(x, y);
      double v = p[0];

      double sensitivity = 0.50;

      if(v > sensitivity) {
        features.push_back(Feature(x, y, p));
      }
    }
  }

  std::sort(features.begin(), features.end());
  std::reverse(features.begin(), features.end());

  int featuresCount = 150;

  int ct=0, index=0;
  while(ct < featuresCount && index < features.size()) {
    bool skip = false;
    Feature ft = features.at(index);

    for(int i=0; i<index; i++) {
      if(ft.closeTo(features.at(i))) {
        skip = true;
        break; // goes to end of for loop
      }
    }

    if(ft.centerX < 12 || ft.centerX > (this->Width() - 12) || ft.centerY < 12 || ft.centerY > (this->Height() - 12)) {
      skip = true;
    }

    if(!skip) {
      featuresOut.push_back(features.at(index));
      ct++;
    }

    index++;
  }

  featuresOut.resize(std::min(int(featuresOut.size()), featuresCount));

  std::cout << "Harris is done" << std::endl;
  return featuresOut;
}

void computeInverseMatrix(double *i, double *output) {
  /** Algorithm based on publically available algorithm found here:
   *  https://forgetcode.com/C/1781-Inverse-Matrix-of-3x3#
   */
  double a[3][3] = {{i[0],i[1],i[2]},{i[3],i[4],i[5]},{i[6],i[7],i[8]}};
  double determinant = 0;

  for(int i=0;i<3;i++) {
    determinant = determinant + (a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]));
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) {
      output[3*j+i] = (((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/ determinant);
    }
  }
}

void computeHomographyMatrix(std::vector<PointCorrespondence> correspondences, double *k) {
    if(correspondences.size() < 4) {
        return;
    }

    int numberOfRows = correspondences.size() * 2;
    double **a = dmatrix(1, numberOfRows, 1, 9);

    for(int i=0; i<correspondences.size(); i++) {
      double x = correspondences.at(i).before.x;
      double y = correspondences.at(i).before.y;
      double u = correspondences.at(i).after.x;
      double v = correspondences.at(i).after.y;

      a[2*i+1][1] = -1 * x;
      a[2*i+1][2] = -1 * y;
      a[2*i+1][3] = -1;
      a[2*i+1][4] = 0;
      a[2*i+1][5] = 0;
      a[2*i+1][6] = 0;
      a[2*i+1][7] = u * x;
      a[2*i+1][8] = u * y;
      a[2*i+1][9] = u;
      a[2*i+2][1] = 0;
      a[2*i+2][2] = 0;
      a[2*i+2][3] = 0;
      a[2*i+2][4] = -1 * x;
      a[2*i+2][5] = -1 * y;
      a[2*i+2][6] = -1;
      a[2*i+2][7] = v * x;
      a[2*i+2][8] = v * y;
      a[2*i+2][9] = v;
    }

    double w[10]; // 1..9

    double **v = dmatrix(1, 9, 1, 9);

    svdcmp(a, numberOfRows, 9, w, v);

    // find the smallest singular value:
    int mi = 1;
    for(int i=2;i<=9;i++) if(w[i]<w[mi]) mi=i;

    // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
    for(int i=1;i<=9;i++) k[i-1]=v[i][mi];
}

void R2Image::
blendImages(R2Image * otherImage,std::vector<Feature> features)
{
  R2Image *output = new R2Image(*this);
//	std::vector<Feature> features = this->Harris(3); // passed by value
  std::vector<Feature>::iterator it;

  int searchSpaceXDim = this->Width() / 10; // half the search space dimension
  int searchSpaceYDim = this->Height() / 10;
  int windowDimension = 12; // half the window dimension

  for(it=features.begin(); it != features.end(); it++) {
    int i, j, m, n;
    double min_ssd = std::numeric_limits<double>::max();
    int min_ssd_x = 0, min_ssd_y = 0;

    Feature ft = *it;

    // Loop through search space

    for(
      i = std::max(ft.centerX - searchSpaceXDim, windowDimension);
      i <= std::min(ft.centerX + searchSpaceXDim, this->Width() - windowDimension);
      i++
    ) {
      for(
        j = std::max(ft.centerY - searchSpaceYDim, windowDimension);
        j <= std::min(ft.centerY + searchSpaceYDim, this->Height() - windowDimension);
        j++
      ) {

        // For each pixel (i, j) in the search space
        double ssd = 0;

        // Calculate the SSD with the feature assuming (i, j) is the center of the new feature
        for(m=-1*windowDimension; m<windowDimension; m++) {
          for(n=-1*windowDimension; n<windowDimension; n++) {
            double oldLuminance = otherImage->Pixel(ft.centerX + m, ft.centerY + n).Luminance();
            double newLuminance = this->Pixel(i + m, j + n).Luminance();
            double diff = oldLuminance - newLuminance;
            ssd += diff * diff;
          }
        }

        // If the computed SSD is lower than the current minimum, set the current minimum to (i, j)
        if(ssd < min_ssd) {
          min_ssd = ssd;
          min_ssd_x = i;
          min_ssd_y = j;
        }
      }
    }

    ft.x2 = min_ssd_x;
    ft.y2 = min_ssd_y;

    // std::cout << "Found match from feature at " << ft.centerX << ", " << ft.centerY << " at point " << ft.x2 << ", " << ft.y2 << std::endl;
    *it = ft;
  }
  std::cout<<"End of for loop"<<std::endl;
  int numberOfTrials = 3000;
  int maxInliers = 0;
  std::vector<int> inlierIndices;
  double* bestK = nullptr;
  double threshold = 8.0;

  srand(time(NULL));
  for(int i=0; i<numberOfTrials; i++) {

    std::vector<int> tempInlierIndices;

    // Randomly select a single track
    int randomIndices[] = {rand() % features.size(), rand() % features.size(), rand() % features.size(), rand() % features.size()};
    std::vector<PointCorrespondence> correspondences;
    correspondences.push_back(createCorrespondence(features.at(randomIndices[0])));
    correspondences.push_back(createCorrespondence(features.at(randomIndices[1])));
    correspondences.push_back(createCorrespondence(features.at(randomIndices[2])));
    correspondences.push_back(createCorrespondence(features.at(randomIndices[3])));

    double *k = (double *) malloc(sizeof(double) * 9);
    computeHomographyMatrix(correspondences, k);

    // Check all other features, and see if their motion vector is similar
    int inliers = 0;

    for(int i=0; i<features.size(); i++) {
      Feature ft = features.at(i);
      double ha_x, ha_y, ha_z;

      // matrix multiplication
      ha_x = ft.centerX * k[0] + ft.centerY * k[1] + 1*k[2];
      ha_y = ft.centerX * k[3] + ft.centerY * k[4] + 1*k[5];
      ha_z = ft.centerX * k[6] + ft.centerY * k[7] + 1*k[8];

      // normalization
      ha_x /= ha_z;
      ha_y /= ha_z;

      double diffVectorLength = abs(sqrt((ha_x-ft.x2)*(ha_x-ft.x2)+(ha_y-ft.y2)*(ha_y-ft.y2)));

      // Count the number of points whose feature match is within a distance
      // threshold of the original point translated by the translation matrix
      if(diffVectorLength < threshold) {
        tempInlierIndices.push_back(i);
        inliers++;
      }
    }

    // If the number of inliers is less than some threshold repeat the above
    if(inliers > maxInliers) {
      maxInliers = inliers;
      bestK = k;
      inlierIndices = tempInlierIndices;
    }
  }

 //std::cout << inlierIndices.size() << " Do we get here?" << std::endl;

 // output->drawLineWithBox(10, 10, 20, 20, 0, 255, 0);

  for(int i=0; i<inlierIndices.size(); i++) {
    Feature f = features.at(inlierIndices.at(i));
    std::cout << f.centerX << "\t" << f.centerY << "\t" << f.x2 << "\t" << f.y2 << std::endl;
    output->drawLineWithBox(f.centerX, f.centerY, f.x2, f.y2, 255, 0, 0);
  }

  /*

  std::vector<PointCorrespondence> bestCorr;
  for(int i=0; i<inlierIndices.size(); i++) {
    bestCorr.push_back(createCorrespondence(features.at(inlierIndices.at(i))));
  }

  double *k = (double *) malloc(sizeof(double) * 9);
  double *invk = (double *) malloc(sizeof(double) * 9);
  computeHomographyMatrix(bestCorr, k);
  computeInverseMatrix(k, invk);

  for(int i=0; i<Width(); i++) {
    for(int j=0; j<Height(); j++) {
      // matrix multiplication
      double inv_x = i * invk[0] + j * invk[1] + invk[2];
      double inv_y = i * invk[3] + j * invk[4] + invk[5];
      double inv_z = i * invk[6] + j * invk[7] + invk[8];

      // normalization
      inv_x /= inv_z;
      inv_y /= inv_z;

      double floating_x = inv_x - floor(inv_x);
      double floating_y = inv_y - floor(inv_y);

      if(inv_x < 0 || inv_x > Width() || inv_y < 0 || inv_y > Height()) {
        continue;
      }

      double r =
        (Pixel(floor(inv_x), floor(inv_y)).Red() * floating_x +
        Pixel(ceil(inv_x),  floor(inv_y)).Red() * (1.0 - floating_x)) * floating_y +
        (Pixel(floor(inv_x), ceil(inv_y)).Red() * floating_x +
        Pixel(ceil(inv_x),  ceil(inv_y)).Red() * (1.0 - floating_x)) * (1.0 - floating_y);

      double g =
        (Pixel(floor(inv_x), floor(inv_y)).Green() * floating_x +
        Pixel(ceil(inv_x),  floor(inv_y)).Green() * (1.0 - floating_x)) * floating_y +
        (Pixel(floor(inv_x), ceil(inv_y)).Green() * floating_x +
        Pixel(ceil(inv_x),  ceil(inv_y)).Green() * (1.0 - floating_x)) * (1.0 - floating_y);

      double b =
        (Pixel(floor(inv_x), floor(inv_y)).Blue() * floating_x +
        Pixel(ceil(inv_x),  floor(inv_y)).Blue() * (1.0 - floating_x)) * floating_y +
        (Pixel(floor(inv_x), ceil(inv_y)).Blue() * floating_x +
        Pixel(ceil(inv_x),  ceil(inv_y)).Blue() * (1.0 - floating_x)) * (1.0 - floating_y);

      bool debug = true;

      if(debug) {
        r += 0.25;
        g += 0.25;
        b += 0.25;
      }

      output->Pixel(i, j).SetRed(output->Pixel(i, j).Red() * 0.5 + r * 0.5);
      output->Pixel(i, j).SetGreen(output->Pixel(i, j).Green() * 0.5 + g * 0.5);
      output->Pixel(i, j).SetBlue(output->Pixel(i, j).Blue() * 0.5 + b * 0.5);
    }
  }

  */

  this->pixels = output->pixels;
  // output->pixels = nullptr;
  // delete output;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


	

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 100, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
	
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






