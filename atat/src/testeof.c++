#include "parse.h"

int main(void) {
  char c;
  cin >> c;
  cout << c << endl;
  cout << cin.eof() << endl;
  cin >> c;
  cout << c << endl;
  cout << cin.eof() << endl;
  cin.putback(c);
  cout << cin.eof() << endl;
  cin >> c;
  cout << c << endl;
  cout << cin.eof() << endl;
}
