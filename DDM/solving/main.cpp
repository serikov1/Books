#include <iostream>
#include <fstream>
#include <freefem++.hpp>

using namespace std;
using namespace Fem2D;

int main() {
    // Определяем область "сапог"
    Mesh Th(Grid(20, 10, Square(10, 10, 1)));
    Th = BuildMesh(Th, {{"Gamma1", {1, 0}, "left"},
                        {"Gamma2", {2, 3}, "right"},
                        {"Gamma3", {3, 0}, "bottom"},
                        {"Gamma4", {0, 1}, "top"}});

    // Определяем конечно-элементное пространство
    Fespace Vh(Th, P1);
    Vh u, v;

    // Определяем правую часть и граничные условия
    func f = 1;
    func g = 0;
    u = 0;
    v = 0;
    problem Puasso(Th, u, v);
    Puasso = integrate(Th, f*v)
             - integrate(Th, dx(u)*dx(v) + dy(u)*dy(v))
             + on(Gamma1, u = 0)
             + on(Gamma2, dx(u) = g)
             + on(Gamma3, u = 0)
             + on(Gamma4, dy(u) = 0);

    // Решаем задачу
    Puasso.solve();

    // Выводим результаты
    ofstream fout("result.txt");
    fout << u << endl;
    fout.close();

    return 0;
}

