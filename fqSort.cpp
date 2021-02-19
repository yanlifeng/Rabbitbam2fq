#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct read {
    string name;
    string seq;
    string tag;
    string qual;

    bool operator<(const read it) const {
        return name < it.name;
    }
};

vector <read> G;

int main(int argc, char *argv[]) {
    freopen(argv[1], "r", stdin);
    string nam, sq, tg, qul;
    while (cin >> nam >> sq >> tg >> qul) {
        G.push_back(read{nam, sq, tg, qul});
    }
    fclose(stdin);
    sort(G.begin(), G.end());
    freopen(argv[2], "w", stdout);
    for (read it:G) {
        cout << it.name << "\n" << it.seq << "\n" << it.tag << "\n" << it.qual << "\n";
    }
    fclose(stdout);
    return 0;
}