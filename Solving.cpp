#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <cmath>
using namespace std;
using namespace std::chrono;

const double PI = 3.1415926535897932;
const int MAX_IT = 100;
const double EPS = 1e-10;

void print(double value, int size, int space, ofstream &fout){
    fout << scientific << setprecision(size) << setw(space) << left << value;
}

double f1(double x, double y){
    return cos(x) - y;
}
// окружность с центром в (0,1) и радиусом 1
double f2(double x, double y){
    return x * x + (y - 1) * (y - 1) - 1.0;
}
void Jacoby(double x, double y, double J[2][2]){
    J[0][0] = -sin(x);
    J[0][1] = -1.0;
    J[1][0] = 2.0 * x;
    J[1][1] = 2.0 * (y - 1.0);
}

double residual(double x, double y){
    double r1 = f1(x, y);
    double r2 = f2(x, y);
    return r1*r1 + r2*r2;
}

void iteration(double x0, double y0, int right_sol, double *x_sol, double *y_sol, int *it, double *r){
    double x = x0, y = y0, x_new, y_new;
    double temp;
    int i;

    for (i = 0; i < MAX_IT; i ++){
        *r = residual(x, y);
        if (*r < EPS)
            break;
        y_new = cos(x);                                // новое у из первого уравнения
        temp = 1.0 - (y_new - 1.0) * (y_new - 1.0);   // новое x из второго уравнения

        // нельзя брать корень в окружности
        if (temp < 0){
            printf("Выход за область определения в МПИ\n");
            break;
        }
        if (right_sol)
            x_new = sqrt(temp);
        else
            x_new = -sqrt(temp);

        x = x_new;
        y = y_new;
    }

    *x_sol = x;
    *y_sol = y;
    *it = i + 1;
}

void newton(double x0, double y0, double *x_sol, double *y_sol, int *it, double *r){
    double x = x0, y = y0;          
    double J[2][2];                         
    double F1, F2, det, dx, dy;                    
    int i;
    
    for (i = 0; i < MAX_IT; i ++){
        *r = residual(x, y);
        if (*r < EPS) 
            break;
        F1 = f1(x, y);
        F2 = f2(x, y);
        Jacoby(x, y, J);       
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (fabs(det) < 1e-16) {
            printf("Вырожденная матрица Якоби\n");
            break;
        }
        // Метод Крамера
        dx = (F1*J[1][1] - F2*J[0][1]) / det;
        dy = (F2*J[0][0] - F1*J[1][0]) / det;
        x -= dx;
        y -= dy;
    }

    *x_sol = x;
    *y_sol = y;
    *it = i + 1;
}

void modified_newton(double x0, double y0, double *x_sol, double *y_sol, int *it, double *r){
    double x = x0, y = y0;          
    double J[2][2];                         
    double F1, F2, det, dx, dy;                    
    int i;
        
    Jacoby(x, y, J);                                       // Считаем один раз для всех итераций
    det = J[0][0]*J[1][1] - J[0][1]*J[1][0]; 
    if (fabs(det) < 1e-16){
        printf("Ошибка: вырожденная матрица\n");
        return;
    }    
    for (i = 0; i < MAX_IT; i ++){
        *r = residual(x, y);
        if (*r < EPS) 
            break;
        F1 = f1(x, y);
        F2 = f2(x, y);
        // Метод Крамера
        dx = (F1*J[1][1] - F2*J[0][1]) / det;
        dy = (F2*J[0][0] - F1*J[1][0]) / det;
        x -= dx;
        y -= dy;
    }
    
    *x_sol = x;
    *y_sol = y;
    *it = i + 1;
}

void discrete_newton(double x0, double y0, double *x_sol, double *y_sol, int *it, double *r){
    double x = x0, y = y0;
    double F1, F2;  
    double J[2][2];                      
    double det, dx, dy, h = 1e-8;                     // h - шаг для разностей
    int i;
    
    for (i = 0; i < MAX_IT; i ++) {
        *r = residual(x, y);
        if (*r < EPS) 
            break;
        F1 = f1(x, y);
        F2 = f2(x, y);
        J[0][0] = (f1(x+h, y) - F1) / h;  // dF1/dx
        J[0][1] = (f1(x, y+h) - F1) / h;  // dF1/dy
        J[1][0] = (f2(x+h, y) - F2) / h;  // dF2/dx
        J[1][1] = (f2(x, y+h) - F2) / h;  // dF2/dy
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (fabs(det) < 1e-16) {
            printf("Ошибка: вырожденная матрица Якоби\n");
            break;
        }
        // Метод Крамера
        dx = (F1*J[1][1] - F2*J[0][1]) / det;
        dy = (F2*J[0][0] - F1*J[1][0]) / det;
        x -= dx;
        y -= dy;
    }    
    *x_sol = x;
    *y_sol = y;
    *it = i + 1;
}


int main(){
    int it;
    double x_left, y_left, x_right, y_right, t, r_left, r_right;
    // Начальное приближение
    double x0_right = 0.9, y0_right = 0.6; 
    double x0_left = -0.9, y0_left = 0.6;

    ofstream fout("solve.txt");
    ofstream fout1("time.txt");

    fout << "---------------------------\n";
    fout << "cos(x) - y = 0\n";
    fout << "x^2 + (y-1)^2 - 1 = 0\n";
    fout << "---------------------------\n";

    fout << "\nЛевое решение\t\t\t\t" << "\t\t\tПравое решение\n";
    fout1 << "\t\t\tЛевое решение\t\t" << "Правое решение\n";
    // МПИ
    fout << "\n--------------------------- МПИ ----------------------------\n\n"; fout1 << "МПИ\t\t\t\t";
    iteration(x0_left, y0_left, 0, &x_left, &y_left, &it, &r_left);
    fout1 << "it: " << it;
    iteration(x0_right, y0_right, 1, &x_right, &y_right, &it, &r_right);
    fout1 << "\t\t\t\tit: " << it;

    fout << "x = "; print(x_left, 8, 20, fout);
    fout << "\t\t\t\tx = "; print(x_right, 8, 20, fout); fout << endl;   
    fout << "y = "; print(y_left, 8, 20, fout);
    fout << "\t\t\t\ty = "; print(y_right, 8, 20, fout); fout << endl;
    fout << "r = "; print(r_left, 3, 20, fout);
    fout << "\t\t\t\tr = "; print(r_right, 3, 20, fout); fout << endl;
    
    // Ньютон
    fout << "\n------------------------- Ньютон ---------------------------\n\n"; fout1 << "\nНьютон\t\t\t";
    newton(x0_left, y0_left, &x_left, &y_left, &it, &r_left);
    fout1 << "it: " << it;
    newton(x0_right, y0_right, &x_right, &y_right, &it, &r_right);
    fout1 << "\t\t\t\tit: " << it;

    fout << "x = "; print(x_left, 8, 20, fout);
    fout << "\t\t\t\tx = "; print(x_right, 8, 20, fout); fout << endl;   
    fout << "y = "; print(y_left, 8, 20, fout);
    fout << "\t\t\t\ty = "; print(y_right, 8, 20, fout); fout << endl;
    fout << "r = "; print(r_left, 3, 20, fout);
    fout << "\t\t\t\tr = "; print(r_right, 3, 20, fout); fout << endl;
    
    // Модифицированный Ньютон
    fout << "\n---------------- Модифицированный Ньютон -------------------\n\n"; fout1 << "\nМод Ньютон\t\t";
    modified_newton(x0_left, y0_left, &x_left, &y_left, &it, &r_left);
    fout1 << "it: " << it;
    modified_newton(x0_right, y0_right, &x_right, &y_right, &it, &r_right);
    fout1 << "\t\t\t\tit: " << it;

    fout << "x = "; print(x_left, 8, 20, fout);
    fout << "\t\t\t\tx = "; print(x_right, 8, 20, fout); fout << endl;   
    fout << "y = "; print(y_left, 8, 20, fout);
    fout << "\t\t\t\ty = "; print(y_right, 8, 20, fout); fout << endl;
    fout << "r = "; print(r_left, 3, 20, fout);
    fout << "\t\t\t\tr = "; print(r_right, 3, 20, fout); fout << endl;
    
    // Дискретный Ньютон
    fout << "\n------------------- Дискретный Ньютон ----------------------\n\n"; fout1 << "\nДис Ньютон\t\t";
    discrete_newton(x0_left, y0_left, &x_left, &y_left, &it, &r_left);
    fout1 << "it: " << it;
    discrete_newton(x0_right, y0_right, &x_right, &y_right, &it, &r_right);
    fout1 << "\t\t\t\tit: " << it;

    fout << "x = "; print(x_left, 8, 20, fout);
    fout << "\t\t\t\tx = "; print(x_right, 8, 20, fout); fout << endl;   
    fout << "y = "; print(y_left, 8, 20, fout);
    fout << "\t\t\t\ty = "; print(y_right, 8, 20, fout); fout << endl;
    fout << "r = "; print(r_left, 3, 20, fout);
    fout << "\t\t\t\tr = "; print(r_right, 3, 20, fout); fout << endl;

    fout.close();
    return 0;
}
