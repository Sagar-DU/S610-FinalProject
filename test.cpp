#include <iostream>

using namespace std;

int main() {
    cout << "Bismillah ArRahmanir Rohim!" << endl;

    int sum = 0;
    for (int i = 1; i <= 100; ++i) {
        sum += i;
    }
    
    cout << "Sum from 1 to 100 is: " << sum << endl;
    return 0;
}