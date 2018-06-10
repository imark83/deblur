#include <iostream>
#include <iomanip>
#include <fstream>

#include <cstdio>

using namespace std;

class Line {
public:
  int k;
  float obj, tv, error, snr, residual;
};


ostream & operator<<(ostream &output, const Line &op) {
  output << setw(4) << op.k << " ";
  output << setw(8) << fixed << setprecision(1)
      << op.obj << " " << setw(8) << op.tv << " ";
  output << setw(4) << op.snr << " ";
  output << scientific << setw(8) << setprecision(2) << op.residual;
  return output;
}

istream& operator>>(istream &input, Line &op) {
  input >> op.k >> op.obj >> op.tv >> op.error >> op.snr >> op.residual;
  return input;
}

int main(int argc, char const *argv[]) {
  ifstream finp;
  char buffer[100];
  Line line;

  finp.open("../testData/alpha05/iadmm-cameraman256.png.txt");
  finp.getline(buffer, 100, '\n');
  printf("%s\n", buffer);

  while(!finp.eof()) {
    finp >> line;
    cout << line << endl;
  }



  finp.close();


  return 0;
}
