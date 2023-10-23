#define _USE_MATH_DEFINES
#include <locale.h>
#include <time.h>
#include <ppl.h>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <iostream>
#include "Task1.h"
#define NNN 100



#define _USE_MATH_DEFINES //Simps_Tst.cpp
#include <iostream>
#include <cmath>
#include <time.h>
#include "Concurrent_Simps.h"
#include <locale.h>
using namespace std;

class My_Class {
private:
	double y;
	const int NNN7 = 200000; //Число разбиений отрезка интегрирования
	const double _Inf = 5000.0; //Фактический верхний предел интегрирования
public:
	My_Class(double _y = 0.0) { y = _y; }
	double Sub_Int_Func(double x)
	{
		double Tmp = 15 * log(10.0) + log(abs(x)) - sqrt(x);
		int N = Tmp > 0 ? ceil(Tmp * Tmp + 1) : 1;
		double P = exp(-abs(x));
		double Tmp2 = cos(x + y) * cos(x + y);
		for (int k = 0; k <= N; k++)
			P *= cos(x / (Tmp2 + exp(sqrt((double)k))));
		return P;
	}

	double Quad()
	{
		return MethodCall::Simpson(0.0, _Inf, NNN7, [this](double x) {return Sub_Int_Func(x); });
	}

	double Concurrent_Quad()
	{
		return MethodCall::Concurrent_Simpson(0.0, _Inf, NNN7, [this](double x) {return Sub_Int_Func(x); });
	}

};


int main(void)
{
	setlocale(LC_ALL, ".ACP");
	double y;
	cout << "y="; cin >> y;

	My_Class TObj(y);
	double Tms = clock();
	double F = TObj.Quad();
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout.precision(8);
	cout << "F=" << F << endl << "Время=" << Tms << " с" << endl;
	Tms = clock();
	F = TObj.Concurrent_Quad();
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "F=" << F << endl << "Время=" << Tms << " с" << endl;
}


// conveyor-transformer.cpp 
// compile with: /EHsc
//#include <agents.h>
//#include <cmath>
//#include <iostream>
//#include <locale.h>
//#include <time.h>
//
//using namespace concurrency;
//using namespace std;
//
//int const N = 100; //Размерность рядов
//int const NNNN = 1000; //Длина массива на входе конвейера 
//int const Chunk = 100; //Порция обрабатываемых элементов
//
//double Func1(double x)//Первая ступень конвейера
//{
//	double Tmp = 0, Tmpa = abs(x);
//	for (int n = 0; n <= N; n++) {
//		double Tmpxn = sqrt(abs(x)+n), Tmpn3 = n * n * n * n * n;
//		for (int k = 0; k <= N; k++) {
//			double Tmpk = k * k;
//			for (int j = 0; j <= N; j++)
//				Tmp += Tmpxn / (sqrt(1 + Tmpa + Tmpn3) + pow( k * k + j * j, 3/2));
//		}
//	}
//	return Tmp;
//}
//
//double Func2(double x)//Вторая ступень конвейера
//{
//	double Tmp = 0, Tmpa = abs(x);
//	for (int n = 0; n <= N; n++) {
//		double Tmpn = n * n;
//		for (int k = 0; k <= N; k++) {
//			double Tmpk = k * k * k;
//			for (int j = 0; j <= N; j++)
//				Tmp += 1.0 / (1 + Tmpa + Tmpn + Tmpk + j * j);
//		}
//	}
//	return x * Tmp;
//}
//
//double Func3(double x)//Третья ступень конвейера
//{
//	double Tmp = 0, Tmp2 = x * x;
//	for (int n = 0; n <= N; n++) {
//		double Tmpn = n * n;
//		for (int k = 0; k <= N; k++) {
//			double Tmpk = k * k;
//			for (int j = 0; j <= N; j++)
//				Tmp += (Tmp2 + j) / (1 + Tmp2 + Tmpn + Tmpk + j * j * j);
//		}
//	}
//	return Tmp;
//}
//
//double Func4(double x)//Четвертая ступень конвейера
//{
//	double Tmp = 0, Tmpa = abs(x);
//	for (int n = 0; n <= N; n++) {
//		double Tmpn = n * n * n;
//		for (int k = 0; k <= N; k++) {
//			double Tmpk = Tmpa - k,
//				Tmpk4 = k * k; Tmpk4 *= Tmpk4;
//			for (int j = 0; j <= N; j++)
//				Tmp += Tmpk / (1 + Tmpa + Tmpn + Tmpk4 + j * j * j);
//		}
//	}
//	return Tmp;
//}
//
//
//
//int main()
//{
//	setlocale(LC_ALL, ".ACP");
//	vector<double> X(NNNN), Y(NNNN), Z(NNNN);
//	//Заполнение входного вектора
//	for (int k = 0; k < X.size(); k++)
//		X[k] = 10 * sin((double)k);
//	double Tms = clock();
//	//Последовательные вычисления
//	for (int k = 0; k < X.size(); k++) {
//		Y[k] = Func1(X[k]);
//		Y[k] = Func2(Y[k]);
//		Y[k] = Func3(Y[k]);
//		Y[k] = Func4(Y[k]);
//	}
//	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
//	cout << "Время последовательного алгоритма: " << Tms << " c." << endl;
//	//Ступени конвейера
//	// t1 -> t2 -> t3 -> t4
//	transformer<double, double> t1(Func1), t2(Func2), t3(Func3), t4(Func4);
//	t1.link_target(&t2);
//	t2.link_target(&t3);
//	t3.link_target(&t4);
//
//	Tms = clock();
//	//Параллельные вычисления
//	for (int k = 0; k < Chunk; k++)
//		send(t1, X[k]);
//
//	for (int k = Chunk; k < X.size(); k++) {
//		Z[k - Chunk] = receive(t4);
//		send(t1, X[k]);
//	}
//	for (int k = X.size() - Chunk; k < X.size(); k++) {
//		Z[k] = receive(t4);
//	}
//	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
//	cout << "Время параллельного алгоритма: " << Tms << " c." << endl;
//}