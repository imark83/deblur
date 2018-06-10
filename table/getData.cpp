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
  float eps = 5e-4f;

  string algorithm[nAlgorithm] = {"admm(05)", "admm(1)", "iadmm"};
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
  Line &admm05 = best[0];
  Line &admm1  = best[1];
  Line &iadmm  = best[2];
  string picFolder("../testData/alpha05/");
  for(int k=0; k<nPic; ++k) {
    for(int j=0; j<nAlgorithm; ++j) {
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
    cout << code[k] << " & " << iadmm << " & " << admm05 << " & "
        << fixed << ((float) admm05.k)/iadmm.k << " & "<< admm1 << " & "
        << fixed << ((float) admm1.k)/iadmm.k << " \\\\" << endl;
  }

  return 0;
}
