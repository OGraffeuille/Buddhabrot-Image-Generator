#include <iostream>
#include <tuple>
#include <algorithm> 
#include <fstream>
#include <complex>
#include <random>
#include <fstream>
#define NOMINMAX // Stops min and max of windows.h interfering
#include <windows.h>
#include <vector>
#include <chrono>
#ifdef _WIN32
#define WINPAUSE system("pause")
#endif

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::seconds;
using namespace std::literals::chrono_literals;
time_point<Clock> start;

// Returns a RGB tuple based on an alpha value between [0 1], and palette chosen
std::tuple<int, int, int> find_rgb(double a, const std::string& palette) {
	a = std::min(std::max(a, 0.), 1.);
	int r, g, b;
	if (palette == "Sunset") {
		r = std::max(std::min((int)(255 - 510 * a), 255), 0);
		g = std::max(std::min((int)(510 * a - 255), 255), 0);
		b = (int)(255 * -a * a * (2 * a - 3));
	} else if (palette == "Gold") {
		r = (int)(128 + 127 * std::cos((a) * 4));
		g = (int)(128 + 127 * std::cos((a + 0.333) * 4));
		b = (int)(128 + 127 * std::cos((a + 0.666) * 4));
	} else if (palette == "Fire") {
		r = std::min((int)(a * (3 - 2 * a) * 255), 255);
		g = (int)(a * a * 250);//std::max(0, (int)(510 * a - 255));
		b = 0;	
	} else { // Greyscale
		r = (int)(255 * a);
		g = (int)(255 * a);
		b = (int)(255 * a);
	}
	return std::make_tuple(r, g, b);
}

void tic() {
	start = Clock::now();
}
int toc() {
	seconds diff = duration_cast<seconds>(Clock::now() - start);
	return diff.count();
}

class Buddhabrot
{
public:
	// Minimum and maximum number of allowable iterations
	int min;
	int max;

	// Number of successful points found and recorded
	int pointsDoneTot;

	// Size of matrix to store counts
	int xDet;
	int yDet;

	// X, Y boundaries of region to store counts (generally centered around 0, -0.4)
	float xMin;
	float xMax;
	float yMin;
	float yMax;

	// Count, histogram matrices
	std::vector<int> count;
	std::vector<int> dist;

	// Class constructor for new matrix
	Buddhabrot(int minIt, int maxIt, const int xDetail, const int yDetail, float xMinLim, float xMaxLim, float yMinLim, float yMaxLim) {
		pointsDoneTot = 0;
		min = std::min(minIt, maxIt);
		max = std::max(minIt, maxIt);
		xDet = xDetail;
		yDet = yDetail;
		xMin = xMinLim;
		xMax = xMaxLim;
		yMin = yMinLim;
		yMax = yMaxLim;
		count.resize(xDet*yDet, 0);
		dist.resize(max + 2, 0);
	}

	// Class constructor when loading old matrix
	Buddhabrot(const std::string& filename) {

		std::ifstream inFile;
		std::ofstream outFile;

		inFile.open(filename);
		if (inFile.fail()) {
			std::cout << "File was not found\n";
			system("PAUSE");
			exit(1);
		}

		std::cout << "Starting loading matrix...\n";
		tic();

		inFile >> pointsDoneTot >> min >> max >> xDet >> yDet >> xMin >> xMax >> yMin >> yMax;
		int val;
		while (!inFile.eof()) {
			inFile >> val;
			count.push_back(val);
		}
		dist.resize(max + 2, 0);

		std::cout << "Finished loading matrix (" << toc() << "s).\n\n";
	}

	// Creates a buddhabrot in the matrix, using calculations with complex numbers (SLOW)
	void CalcComplex(int pointsToDo) {
		std::cout << "Starting calculating matrix...\n";
		tic();

		// Seed, random number generator, distribution on which to apply the generator
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::uniform_int_distribution<long unsigned> distribution(0, 1e9);

		int dispInterval = pow(10, std::max((int)(log10(pointsToDo) - 2), 0));

		int pointsDone = 0;
		while (pointsDone < pointsToDo) {

			// Find sensible start point
			double a = distribution(generator) / 1.e9;
			double b = distribution(generator) / 1.e9;
			std::complex<double> c0(a * 4 - 2, b * 4 - 2);
			double q = pow(real(c0) - 0.25, 2) + pow(imag(c0), 2);

			while ((pow(real(c0) + 1, 2) + pow(imag(c0), 2) < 0.0625)		// Testing if in small circle
				|| (q*(q + real(c0) - 0.25) < 0.25 * pow(imag(c0), 2))	// Testing if in main cardioid
				|| (abs(c0) > 2)) {										// Testing if outside large circle

				a = distribution(generator) / 1.e9;
				b = distribution(generator) / 1.e9;
				c0 = std::complex<double>(a * 4 - 2, b * 4 - 2);
				q = pow(real(c0) - 0.25, 2) + pow(imag(c0), 2);
			}

			// Iterate until max is reached, or it diverges
			std::complex<double> c = c0;
			int iter = 0;

			while (iter <= max) {
				c = c * c + c0;
				if (abs(c) > 2) {
					break;
				}
				iter++;
			}

			// Add point to matrix if it diverges at correct moment
			if (iter >= min && iter <= max) {
				c = c0;
				for (int i = 0; i <= iter; i++) {
					int x = (int)(xDet * (real(c) - xMin) / (xMax - xMin));
					int y = (int)(yDet * (imag(c) - yMin) / (yMax - yMin));
					if (x >= 0 && x < xDet && y >= 0 && y < yDet) {
						count[x*yDet + y]++;
					}
					c = c * c + c0;
				}
				pointsDoneTot++;
				pointsDone++;
				if (pointsDone % dispInterval == 0) {
					std::cout << "Completed: " << pointsDone << " out of " << pointsToDo << " points (" << 100.*pointsDone/pointsToDo << "%, " << pointsDoneTot << " total points)" << std::endl;
				}
			}
		}
		std::cout << "Finished calculating matrix (" << toc() << "s).\n\n";
	}

	// Creates a buddhabrot in the matrix, using calculations with doubles (FAST)
	void CalcDouble(int pointsToDo, bool countDist) {
		std::cout << "Starting calculating matrix...\n";
		tic();

		// Seed, random number generator, distribution on which to apply the generator
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::uniform_int_distribution<long unsigned> distribution(0, 1e9);

		int dispInterval = pow(10, std::max((int)(log10(pointsToDo) - 2), 0));

		double a, b, re, re0, im, im0, q, re_, iter, x, y;
		int pointsDone = 0;
		while (pointsDone < pointsToDo) {

			// Find sensible start point
			a = distribution(generator) / 1.e9;
			b = distribution(generator) / 1.e9;
			re = a * 4 - 2;
			im = b * 4 - 2;
			q = (re-0.25)*(re-0.25) + im*im;
			if (countDist) dist[0]++;

			while (((re+1)*(re+1) + im*im < 0.0625)		// Testing if in small circle
				|| (q*(q + re - 0.25) < 0.25 * im*im)	// Testing if in main cardioid
				|| (re*re + im*im > 4)) {				// Testing if outside large circle
				
				a = distribution(generator) / 1.e9;
				b = distribution(generator) / 1.e9;
				re = a * 4 - 2;
				im = b * 2 - 2;
				q = (re-0.25)*(re-0.25) + im*im;
				if (countDist) dist[0]++;
			}

			// Iterate until max is reached, or it diverges
			re0 = re;
			im0 = im;
			for (iter = -1; (iter <= max) && (re*re + im*im < 4); iter++) {
				re_ = re;
				re *= re;
				re += re0 - im * im;
				im *= 2 * re_;
				im += im0;
			}

			//std::cout << max << " " << iter << "\n";
			if (countDist) dist[iter]++;

			// Add point to matrix if it diverges at correct moment
			if (iter >= min && iter <= max) {
				re = re0;
				im = im0;
				for (int i = 0; i <= iter; i++) {
					x = (int)(xDet * (re - xMin) / (xMax - xMin));
					y = (int)(yDet * (im - yMin) / (yMax - yMin));
					if (x >= 0 && x < xDet && y >= 0 && y < yDet) {
						count[x*yDet + y]++;
					}
					re_ = re;
					re *= re;
					re += re0 - im * im;
					im *= 2 * re_;
					im += im0;
				}
				pointsDoneTot++;
				pointsDone++;
				if (pointsDone % dispInterval == 0) {
					std::cout << "Completed: " << pointsDone << " out of " << pointsToDo << " points (" << 100.*pointsDone / pointsToDo << "%, " << pointsDoneTot << " total points)" << std::endl;
				}
			}
		}
		std::cout << "Finished calculating matrix (" << toc() << "s).\n\n";
	}

	// Fills count matrix with a gradient
	void TestMatrix() {

		for (int i = 0; i < xDet*yDet; i++) {
			count[i] = int(255 * (float(i%xDet) / (xDet - 1) + float(i/yDet) / (yDet - 1)) * 0.5);
		}
	}

	// Prints the matrix contents into the console
	void DispMatrix() {

		int maxChars = int(log10(FindMaxVal()) + 2);
		for (int y = yDet - 1; y >= 0; y--) {
			for (int x = 0; x < xDet; x++) {
				int chars = (count[x] > 0) ? int(log10(count[x]) + 1) : 1;
				for (int k = 0; k < maxChars - chars; k++) {
					std::cout << " ";
				}
				std::cout << int(count[x]);
			}
			std::cout << std::endl;
		}
	}

	// Finds greatest value in matrix
	int FindMaxVal() {
		int maxVal = -1;
		for (int i = 0; i < xDet*yDet; i++) {
			if (count[i] > maxVal) {
				maxVal = count[i];
			}
		}
		//std::cout << "Max val is: " << maxVal << std::endl;
		return maxVal;
	}

	// Find percentile value, for example, if the matrix is 100x100, and percentile = 80, this finds the 2000th largest value
	int FindPercentileVal(float percentile) {

		int maxVal = FindMaxVal();

		int* histogram = new int[maxVal + 1];
		for (int i = 0; i <= maxVal; i++) {
			histogram[i] = 0;
		}
		for (int i = 0; i < xDet*yDet; i++) {
			histogram[count[i]]++;
		}

		int pointsTot = xDet * yDet;
		int pointsLeft = xDet * yDet * percentile / 100;
		for (int i = maxVal; i > 0; i--) {
			pointsTot -= histogram[i];
			if (pointsTot <= pointsLeft) {
				return i;
			}
		}
		return 1;
	}

	// Creates a PPM file of the matrix
	void PrintPPM(float percentile, const std::string& palette, const std::string& filename = "Buddhabrot") {
		std::cout << "Saving PPM image file...\n";
		tic();

		std::ostringstream oss;
		oss << filename << " [" << xDet << "x" << yDet << "] [" << min << "-" << max << "] " << pointsDoneTot << " points.ppm";
		std::string fullFilename = oss.str();

		std::ofstream img(fullFilename);
		img << "P3" << std::endl;
		img << xDet << " " << yDet << std::endl;
		img << "255" << std::endl;

		double a;
		int cutoff = FindPercentileVal(percentile);
		for (int y = yDet - 1; y >= 0; y--) {
			for (int x = 0; x < xDet; x++) {
				a = std::min(count[x*yDet + y], cutoff) / (cutoff*1.);
				std::tuple<int, int, int> rgb = find_rgb(a, palette);
				img << std::get<0>(rgb) << " " << std::get<1>(rgb) << " " << std::get<2>(rgb) << std::endl;
			}
		}
		std::cout << "Finished saving image (" << toc() << "s).\n\n";
	}

	// Creates a BMP file of the matrix
	void PrintBMP(float percentile, const std::string& palette, const std::string& filename = "Buddhabrot") {
		std::cout << "Saving BMP image file...\n";
		tic();

		unsigned int headers[13];
		FILE * outfile;
		int extrabytes;
		int paddedsize;
		int x; int y; int n;
		int red, green, blue;

		extrabytes = 4 - ((xDet * 3) % 4);					// How many bytes of padding to add to each
		if (extrabytes == 4)                                // horizontal line - the size of which must
			extrabytes = 0;									// be a multiple of 4 bytes.

		paddedsize = ((xDet * 3) + extrabytes) * yDet;

		// Headers
		// Note that the "BM" identifier in bytes 0 and 1 is NOT included in these "headers".
		headers[0] = paddedsize + 54;      // bfSize (whole file size)
		headers[1] = 0;                    // bfReserved (both)
		headers[2] = 54;                   // bfOffbits
		headers[3] = 40;                   // biSize
		headers[4] = xDet;                 // biWidth
		headers[5] = yDet;                 // biHeight

										   // Would have biPlanes and biBitCount in position 6, but they're shorts.
										   // It's easier to write them out separately (see below) than pretend
										   // they're a single int, especially with endian issues...

		headers[7] = 0;                    // biCompression
		headers[8] = paddedsize;           // biSizeImage
		headers[9] = 0;                    // biXPelsPerMeter
		headers[10] = 0;                   // biYPelsPerMeter
		headers[11] = 0;                   // biClrUsed
		headers[12] = 0;                   // biClrImportant

		std::ostringstream oss;
		oss << filename << " [" << xDet << "x" << yDet << "] [" << min << "-" << max << "] " << pointsDoneTot << " points.bmp";
		std::string fullFilename = oss.str();
		outfile = fopen(fullFilename.c_str(), "wb");

		// Begin writing headers
		// When printing ints and shorts, we write out 1 character at a time to avoid endian issues
		fprintf(outfile, "BM");

		for (n = 0; n <= 5; n++)
		{
			fprintf(outfile, "%c", headers[n] & 0x000000FF);
			fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
			fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
			fprintf(outfile, "%c", (headers[n] & (unsigned int)0xFF000000) >> 24);
		}

		// These next 4 characters are for the biPlanes and biBitCount fields.
		fprintf(outfile, "%c", 1);
		fprintf(outfile, "%c", 0);
		fprintf(outfile, "%c", 24);
		fprintf(outfile, "%c", 0);

		for (n = 7; n <= 12; n++) {
			fprintf(outfile, "%c", headers[n] & 0x000000FF);
			fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
			fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
			fprintf(outfile, "%c", (headers[n] & (unsigned int)0xFF000000) >> 24);
		}

		// Write data after headers are completed
		double a;
		int cutoff = FindPercentileVal(percentile);
		for (x = xDet - 1; x >= 0; x--) {
			for (y = yDet - 1; y >= 0; y--) {					// BMP image format is written from bottom to top

				a = std::min(count[x*yDet + y], cutoff) / (cutoff*1.);
				std::tuple<int, int, int> rgb = find_rgb(a, palette);

				fprintf(outfile, "%c", std::get<2>(rgb));		// Backwards because written in bgr format
				fprintf(outfile, "%c", std::get<1>(rgb));
				fprintf(outfile, "%c", std::get<0>(rgb));
			}
			if (extrabytes) {									// BMP lines must be of lengths divisible by 4
				for (n = 1; n <= extrabytes; n++) {
					fprintf(outfile, "%c", 0);
				}
			}
		}
		fclose(outfile);
		std::cout << "Finished saving image (" << toc() << "s).\n\n";
	}

	// Creates a BMP file of the matrix (smoothed with a kernel, untested)
	void PrintBMPSmooth(float percentile, const std::string& palette, const std::string& filename = "Buddhabrot") {
		std::cout << "Saving BMP image file...\n";
		tic();

		unsigned int headers[13];
		FILE * outfile;
		int extrabytes;
		int paddedsize;
		int x; int y; int n;
		int red, green, blue;

		extrabytes = 4 - ((xDet * 3) % 4);					// How many bytes of padding to add to each
		if (extrabytes == 4)                                // horizontal line - the size of which must
			extrabytes = 0;									// be a multiple of 4 bytes.

		paddedsize = ((xDet * 3) + extrabytes) * yDet;

		// Headers
		// Note that the "BM" identifier in bytes 0 and 1 is NOT included in these "headers".
		headers[0] = paddedsize + 54;      // bfSize (whole file size)
		headers[1] = 0;                    // bfReserved (both)
		headers[2] = 54;                   // bfOffbits
		headers[3] = 40;                   // biSize
		headers[4] = xDet;                 // biWidth
		headers[5] = yDet;                 // biHeight

										   // Would have biPlanes and biBitCount in position 6, but they're shorts.
										   // It's easier to write them out separately (see below) than pretend
										   // they're a single int, especially with endian issues...

		headers[7] = 0;                    // biCompression
		headers[8] = paddedsize;           // biSizeImage
		headers[9] = 0;                    // biXPelsPerMeter
		headers[10] = 0;                   // biYPelsPerMeter
		headers[11] = 0;                   // biClrUsed
		headers[12] = 0;                   // biClrImportant

		std::ostringstream oss;
		oss << filename << " [" << xDet << "x" << yDet << "] [" << min << "-" << max << "] " << pointsDoneTot << " points.bmp";
		std::string fullFilename = oss.str();
		outfile = fopen(fullFilename.c_str(), "wb");

		// Begin writing headers
		// When printing ints and shorts, we write out 1 character at a time to avoid endian issues
		fprintf(outfile, "BM");

		for (n = 0; n <= 5; n++)
		{
			fprintf(outfile, "%c", headers[n] & 0x000000FF);
			fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
			fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
			fprintf(outfile, "%c", (headers[n] & (unsigned int)0xFF000000) >> 24);
		}

		// These next 4 characters are for the biPlanes and biBitCount fields.
		fprintf(outfile, "%c", 1);
		fprintf(outfile, "%c", 0);
		fprintf(outfile, "%c", 24);
		fprintf(outfile, "%c", 0);

		for (n = 7; n <= 12; n++) {
			fprintf(outfile, "%c", headers[n] & 0x000000FF);
			fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
			fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
			fprintf(outfile, "%c", (headers[n] & (unsigned int)0xFF000000) >> 24);
		}


		std::vector<int> countSmooth;
		countSmooth.resize(xDet*yDet, 0);
		for (int i = xDet + 1; x < xDet*(yDet - 1) - 1; i++) {
			countSmooth[i] = (count[i - xDet - 1] + 2 * count[i - xDet] + count[i - xDet + 1] + 2 * count[i - 1] + 4 * count[i] + 2 * count[i + 1] + count[i + xDet - 1] + 2 * count[i + xDet] + count[i + xDet + 1]) / 16.;
		}

		// Write data after headers are completed
		double a;
		int cutoff = FindPercentileVal(percentile);
		for (x = xDet - 1; x >= 0; x--) {
			for (y = yDet - 1; y >= 0; y--) {					// BMP image format is written from bottom to top

				a = std::min(countSmooth[x*yDet + y], cutoff) / (cutoff*1.);
				std::tuple<int, int, int> rgb = find_rgb(a, palette);

				fprintf(outfile, "%c", std::get<2>(rgb));		// Backwards because written in bgr format
				fprintf(outfile, "%c", std::get<1>(rgb));
				fprintf(outfile, "%c", std::get<0>(rgb));
			}
			if (extrabytes) {									// BMP lines must be of lengths divisible by 4
				for (n = 1; n <= extrabytes; n++) {
					fprintf(outfile, "%c", 0);
				}
			}
		}
		fclose(outfile);
		std::cout << "Finished saving image (" << toc() << "s).\n\n";
	}

	// Creates a txt file of matrix
	void SaveMatrix(const std::string& filename = "Buddhabrot") {
		std::cout << "Saving matrix as text file...\n";
		tic();

		std::ostringstream oss;
		oss << filename << " [" << xDet << "x" << yDet << "] [" << min << "-" << max << "] " << pointsDoneTot << " points.txt";
		std::string fullFilename = oss.str();

		std::ofstream outFile(fullFilename);
		outFile << pointsDoneTot << "\n" << min << "\n" << max << "\n" << xDet << "\n" << yDet << "\n" << xMin << "\n" << xMax << "\n" << yMin << "\n" << yMax << "\n";
		for (int i = 0; i < xDet*yDet; i++) {
			outFile << count[i] << "\n";
		}
		std::cout << "Finished saving matrix text file " << toc() << "s).\n\n";
	}

	// Creates a txt file of distribution
	void SaveDist(const std::string& filename = "Buddhabrot Dist") {
		std::cout << "Saving distribution as text file...\n";
		tic();

		std::ostringstream oss;
		oss << filename << " [1-" << max << "] " << dist[0] << " points.txt";
		std::string fullFilename = oss.str();
		std::ofstream outFile(fullFilename);

		for (int i = 0; i < max+2; i++) {
			outFile << dist[i] << '\n';
		}
		std::cout << "Finished saving distribution text file (" << toc() << "s).\n\n";
	}
};

int main() {

	// Create new buddhabrot, or load existing text file
	Buddhabrot BB(2e4, 6e4, 2500, 2500, -1.6, 1.2, -1.4, 1.4);
	//Buddhabrot BB("Sizetest [2500x2500] [20000-60000] 40 points.txt");

	// Calculate points in Buddhabrot
	BB.CalcDouble(40, false);

	// Save Buddhabrot as text file and as image
	BB.SaveMatrix("Test");
	BB.PrintBMP(99.2, "Gold", "Test");


	std::cout << "Press any key to continue..." << std::endl;
	std::getchar();
	return 0;
}