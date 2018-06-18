#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>



using namespace std;

class Line {
public:
  int k;
  float obj, tv, error, snr, residual;
  Line() {
    k = -1;
    obj = tv = error = residual = 1.0f/0.0f;
    snr = 0.0f;
  }
};


ostream & operator<<(ostream &output, const Line &op) {
  output << setw(4) << op.k << " & ";
  output << setw(8) << fixed << setprecision(1) << op.error << " & ";
  output << setw(4) << op.snr << " & ";
  output << scientific << setw(8) << setprecision(2) << op.residual;
  return output;
}

istream& operator>>(istream &input, Line &op) {
  input >> op.k >> op.obj >> op.tv >> op.error >> op.snr >> op.residual;
  return input;
}

int main(int argc, char const *argv[]) {

  int nPic = 6;
  int nAlgorithm = 3;
  float eps = 1e-3f;

  string algorithm[nAlgorithm] = {"iadmm", "iadmm", "admm(1)"};
  string picName[nPic] = {
    "cameraman256.png",
    "lena256.png",
    "brain-512.png",
    "heart-512.png",
    "5.3.01-1024.png",
    "5.3.02-1024.png"
    // "5.2.08-512.png",
    // "4.2.03-512.png",
  };
  string code[nPic] = {
    "L0",
    "L1",
    "M0",
    "M1",
    "H0",
    "H1"
    // "5.2.08-512.png",
    // "4.2.03-512.png",
  };

  ifstream finp;
  char c_buffer[100];
  Line line, next_line;
  Line best[nAlgorithm];
  Line &iadmm05 = best[0];
  Line &iadmm02 = best[1];
  Line &admm1   = best[2];
  string picFolder;
  for(int k=0; k<nPic; ++k) {
    for(int j=0; j<nAlgorithm; ++j) {
      if (j==1)
        picFolder = "../testData/alpha02/";
      else
        picFolder = "../testData/alpha05/";

      string filename = picFolder + algorithm[j] + "-"
            + picName[k] + ".txt";
      finp.open(filename.c_str());
      // finp.open("../testData/alpha05/admm(05)-cameraman256.png.txt");
      // finp.open("../testData/alpha05/admm(1)-cameraman256.png.txt");
      // finp.open("../testData/alpha05/iadmm-cameraman256.png.txt");
      finp.getline(c_buffer, 100, '\n');

      line = Line();
      finp >> next_line;
      while(!finp.eof() && next_line.residual < line.residual
                  && line.residual > eps) {
        line = next_line;
        // cout << line << endl;
        finp >> next_line;
      }

      best[j] = line;
      finp.close();

    }
    // cout << c_buffer << endl;
    // cout << admm05 << endl << admm1 << endl << iadmm << endl;
    if(k%2)
      cout << "\\rowcolor[gray]{0.9}" << endl;
    cout << code[k] << " & " << iadmm05 << " & " << iadmm05 << " & "
        << fixed << ((float) iadmm02.k)/iadmm05.k << " & "<< admm1 << " & "
        << fixed << ((float) admm1.k)/iadmm05.k << " \\\\" << endl;
  }

  return 0;
}
